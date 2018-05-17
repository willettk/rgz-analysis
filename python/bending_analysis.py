import logging, urllib2, time, json, os
import pymongo, bson
import numpy as np
import StringIO, gzip
import csv
from astropy.io import fits
from astropy import wcs, coordinates as coord, units as u
from astropy.cosmology import Planck13 as cosmo
from scipy.optimize import brentq, curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy import stats
from skgof import ad_test
from sklearn.metrics import r2_score

from consensus import rgz_path, data_path, db, version
from processing import *
import contour_path_object as cpo
completed_file = '%s/bending_completed%s.txt' % (rgz_path, version)

# For internal use
from pprint import pprint
import itertools
from sklearn.decomposition import PCA
from corner import corner
import matplotlib.pyplot as plt
from matplotlib import rc
plt.ion()

# Connect to Mongo database
version = '_bending'
subjects = db['radio_subjects']
consensus = db['consensus{}'.format(version)]
catalog = db['catalog{}'.format(version)]
whl = db['WHL15']
rm_m = db['redmapper_members']
rm_c = db['redmapper_clusters']
amf = db['AMFDR9']
version = ''
bent_sources = db['bent_sources{}'.format(version)]
bending_15 = db['bending_15{}'.format(version)]
bending_control = db['bending_control{}'.format(version)]

# Final sample cuts
total_cuts = {'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'RGZ.radio_consensus':{'$gte':0.65}, 'using_peaks.bending_corrected':{'$lte':135.}}
params = total_cuts.copy()
params['WHL.population'] = 'outer'
bend = []
for i in bending_15.find(params):
	bend.append(i['using_peaks']['bending_corrected'])
median_bend = np.median(bend)

# Get dictionary for finding the path to FITS files and WCS headers
with open('%s/first_fits.txt' % rgz_path) as f:
	lines = f.readlines()

pathdict = {}
for l in lines:
	spl = l.split(' ')
	pathdict[spl[1].strip()] = '%s/rgz/raw_images/RGZ-full.%i/FIRST-IMGS/%s.fits' % (data_path, int(spl[0]), spl[1].strip())

### Functions ###

def get_data(subject):
	'''
	Returns the radio contours belonging to a single subject field
	'''
	
	link = subject['location']['contours'] # Gets url as Unicode string
	
	# Use local file if available
	
	jsonfile = link.split("/")[-1]
	jsonfile_path = "{0}/rgz/contours/{1}".format(data_path,jsonfile)
	if os.path.exists(jsonfile_path):
		with open(jsonfile_path,'r') as jf:
			data = json.load(jf)
	
	# Otherwise, read from web
	
	else:
		
		# Reform weblink to point to the direct S3 URL, which will work even with older SSLv3
		
		link_s3 = "http://zooniverse-static.s3.amazonaws.com/"+link.split('http://')[-1]
		
		tryCount = 0
		while(True): # In case of error, wait 10 sec and try again; give up after 5 tries
			tryCount += 1
			try:
				compressed = urllib2.urlopen(str(link_s3)).read() #reads contents of url to str
				break
			except (urllib2.URLError, urllib2.HTTPError) as e:
				if tryCount>5:
					output('Unable to connect to Amazon Web Services; trying again in 10 min', logging.exception)
					raise fn.DataAccessError(message)
				logging.exception(e)
				time.sleep(10)
		
		tempfile = StringIO.StringIO(compressed) # Temporarily stores contents as file (emptied after unzipping)
		uncompressed = gzip.GzipFile(fileobj=tempfile, mode='r').read() # Unzips contents to str
		data = json.loads(uncompressed) # Loads JSON object
	
	return data

def get_contours(w, ir_pos, peak_pos, data, peak_count):
	'''
	Returns a list of Path objects corresponding to each outer contour in the data, in RA and dec coordinates
	Removes outer layers until there are two components and a disjoint IR
	Removes any contours that aren't in this source
	'''
	assert (peak_count in [2,3]), 'Not a valid morphology'
	
	# Assemble the contour trees
	contour_trees = []
	for contour in data['contours']:
		tree = cpo.Node(w, contour=contour)
		contour_trees.append(tree)
	
	# Remove each contour that doesn't contain a peak from this source
	contains_peak = []
	for ix, tree in enumerate(contour_trees):
		if any(tree.contains(peak_pos)):
			contains_peak.append(ix)
	
	contour_trees[:] = [tree for ix, tree in enumerate(contour_trees) if ix in contains_peak]
	
	# Combine the entire source into a single tree
	value_at_inf = {'arr':[{'x':-1,'y':-1}, {'x':w._naxis1+1,'y':-1}, {'x':w._naxis1+1,'y':w._naxis2+1}, {'x':-1,'y':w._naxis2+1}, {'x':-1,'y':-1}], 'k':-1}
	source_tree = cpo.Node(w)
	source_tree.insert(cpo.Node(w, value=value_at_inf))
	source_tree.children = contour_trees
	
	# Remove the central source if it's a triple
	if peak_count == 3:
		source_tree.remove_triple_center(ir_pos, peak_pos)
	
	# Increase the contour level until the IR position is outside all the contours
	roots = []
	source_tree.get_equal_disjoint(ir_pos, roots)
	source_tree.children = roots
	
	return source_tree

def get_pos_angle(w, ir, method, method_data):
	'''
	Determines the position angle between the IR position and the given comparison object
	Method:
		Contour: from the IR position to the most distant part of the radio contour (data is contour_list)
		Peak: from the IR position to the peak of each component (data is source['radio']['peaks'])
	'''
	if method == 'contour':
		contour_sky = coord.SkyCoord(w.wcs_pix2world(method_data.vertices,1), unit=(u.deg,u.deg), frame='icrs')
		separation = ir.separation(contour_sky)
		pos_angle = ir.position_angle(contour_sky)[separation==np.max(separation)][0]
	elif method == 'peak':
		pos_angle = ir.position_angle(coord.SkyCoord(method_data['ra'], method_data['dec'], unit=(u.deg,u.deg), frame='icrs'))
	return pos_angle

def get_angles(w, ir, method, method_data):
	'''
	Determines the opening angle of the radio tail and the position angle of the angle bisector
	Method:
		Contour: from the IR position to the most distant part of the radio contour (data is contour_list)
		Peak: from the IR position to the peak of each component (data is source['radio']['peaks'])
	'''
	assert (method in ['contour', 'peak']), 'Not a valid method'
	pos_angle_0 = get_pos_angle(w, ir, method, method_data[0])
	pos_angle_1 = get_pos_angle(w, ir, method, method_data[1])
	opening_angle = np.abs(pos_angle_1-pos_angle_0).wrap_at(2*np.pi*u.rad)
	bending_angle = coord.Angle(np.abs(np.pi*u.rad - opening_angle))
	bisector = (pos_angle_1+pos_angle_0)/2.
	if np.abs(bisector-pos_angle_0) > np.pi/2*u.rad:
		bisector += np.pi*u.rad
	bending_angles = {'pos_angle_0':pos_angle_0, 'pos_angle_1':pos_angle_1, 'bending_angle':bending_angle, 'bisector':bisector.wrap_at(2*np.pi*u.rad)}
	return bending_angles

def get_global_peaks(w, peak_pos, peaks, contour_tree):
	'''
	Determines the position of the global maximum for each component in the contour
	'''
	global_peaks = []
	for child in contour_tree.children:
		global_peak = {'flux':0}
		for peak in [peaks[ix] for ix, elem in enumerate(child.contains(peak_pos)) if elem]:
			if peak['flux'] > global_peak['flux']:
				global_peak = peak
		if global_peak['flux'] > 0:
			global_peaks.append(global_peak)
	return global_peaks

def curve_intersect(fun1, fun2, xmin, xmax):
	'''
	Finds the intersection of two curves, bounded in [xmin, xmax]
	Returns an array of x values
	'''
	diff = lambda x: fun1(x)-fun2(x)
	x_range = np.linspace(xmin, xmax, 100)
	m_sign = np.sign(diff(x_range)).astype(int)
	roots = x_range[np.where(m_sign[1:] - m_sign[:-1] != 0)[0] + 1]
	
	# If they don't cross, return None
	if len(roots) == 0:
		return np.array([])
	
	# If they cross exactly once, find thcone global solution
	elif len(roots) == 1:
		return np.array([brentq(diff, xmin, xmax)])
	
	# If they cross multiple times, find the local solution between each root
	else:
		limits = np.concatenate(([xmin], roots, [xmax]))
		intersections = np.empty(len(limits)-2)
		for ix in range(len(intersections)):
			intersections[ix] = brentq(diff, limits[ix], limits[ix+1])
		return intersections

def get_colinear_separation(w, ir, peak, contour):
	'''
	Finds the distance from the host to the edge of the contour, passing through the peak
	'''
	
	ir_pos = w.wcs_world2pix(np.array([[ir.ra.deg,ir.dec.deg]]), 1)[0]
	peak_pos = w.wcs_world2pix(np.array([[peak['ra'], peak['dec']]]), 1)[0]
	
	# Extrapolate the line connecting the peak to the IR position
	slope = (peak_pos[1]-ir_pos[1])/(peak_pos[0]-ir_pos[0])
	extrap_pos = ir_pos + w._naxis1*np.array([1.,slope])
	extrap_neg = ir_pos - w._naxis1*np.array([1.,slope])
	
	# Split the contours into well-behaved functions
	# Roll the array until the first index is the minimum value
	x, y = contour.vertices.T
	xmin_loc = np.where(x==min(x))[0][0]
	x_rot = np.append(np.roll(x[:-1], len(x)-xmin_loc-1), min(x))
	y_rot = np.append(np.roll(y[:-1], len(x)-xmin_loc-1), y[xmin_loc])
	
	# Find where the contour doubles back on itself along the x-axis
	m_sign = np.sign(x_rot[1:]-x_rot[:-1])
	roots = np.where(m_sign[1:] - m_sign[:-1] != 0)[0] + 1
	limits = np.concatenate(([0], roots, [len(x_rot)-1]))
	
	# Split the contours at the double-back positions
	domains = []
	ranges = []
	for ix in range(len(limits)-1):
		domains.append(x_rot[limits[ix]:limits[ix+1]+1])
		ranges.append(y_rot[limits[ix]:limits[ix+1]+1])
	
	# Interpolate the contour segments
	c_interps = []
	for x_seg, y_seg in zip(domains, ranges):
		c_interp = interp1d(x_seg, y_seg, 'linear')
		c_interps.append(c_interp)
	
	if peak_pos[0] > ir_pos[0]:
		tail = np.vstack((extrap_neg, ir_pos, peak_pos, extrap_pos))
	else:
		tail = np.vstack((extrap_pos, ir_pos, peak_pos, extrap_neg))
	
	tail_interp = interp1d(tail.T[0], tail.T[1], 'linear')
	
	# Find the intersections of the contours and tail
	x_intersects, y_intersects = [], []
	for ix, c_interp in enumerate(c_interps):
		x_intersect = curve_intersect(tail_interp, c_interp, domains[ix][0], domains[ix][-1])
		y_intersect = c_interp(x_intersect)
		x_intersects.append(x_intersect)
		y_intersects.append(y_intersect)
	
	intersects = np.vstack((np.hstack(x_intersects), np.hstack(y_intersects))).T
	
	# Return the maximum separation between host and edge
	intersects_sky = coord.SkyCoord(w.wcs_pix2world(intersects,1), unit=(u.deg,u.deg), frame='icrs')
	return max(ir.separation(intersects_sky))

def get_tail_lengths(w, ir, method, contour_list, peaks=None):
	'''
	Determines angular separation between the IR position and the given comparison object
	Method:
		Contour: from the IR position to the most distant part of the radio contour
		Peak: from the IR position to the peak of the component
	'''
	tail_lengths = []
	if method == 'contour':
		for contour in contour_list:
			contour_sky = coord.SkyCoord(w.wcs_pix2world(contour.vertices,1), unit=(u.deg,u.deg), frame='icrs')
			separation = ir.separation(contour_sky)
			tail_lengths.append(np.max(separation))
	elif method == 'peak':
		assert (peaks is not None), 'No radio peaks provided'
		for contour, peak in zip(contour_list, peaks):
			tail_lengths.append(get_colinear_separation(w, ir, peak, contour))
	return tail_lengths

def peak_edge_ratio(w, ir, peaks, tails):
	'''
	Calculate the ratio of the distance to the peak and to the edge of each tail (measured on the sky)
	'''
	ratios = []
	for peak, tail in zip(peaks, tails):
		peak_pos = coord.SkyCoord(peak['ra'], peak['dec'], unit=(u.deg,u.deg), frame=('icrs'))
		ratios.append(ir.separation(peak_pos).deg/tail.deg)
	return ratios

def get_z(source):
	'''
	Returns the best redshift value and uncertainty for the source
	'''
	if 'SDSS' in source and 'spec_redshift' in source['SDSS']:
		return source['SDSS']['spec_redshift'], source['SDSS']['spec_redshift_err']
	elif 'SDSS' in source and 'photo_redshift' in source['SDSS']:
		return source['SDSS']['photo_redshift'], source['SDSS']['photo_redshift_err']
	elif 'AllWISE' in source and 'photo_redshift' in source['AllWISE']:
		return source['AllWISE']['photo_redshift'], 0

def get_whl(ir, z, z_err, transverse, dz):
	'''
	Find the corresponding galaxy cluster in the WHL15 catalog
	If multiple clusters match, choose the one with least angular separation
	'''
	
	# If the galaxy is too close, physical separations become too great on the sky
	# Restrict redshifts to at least 0.01
	if z < 0.01:
		return None
	
	# Maximum separation
	max_sep = float(transverse * u.Mpc / cosmo.angular_diameter_distance(z) * u.rad / u.deg)
	
	best_sep = np.inf
	cluster = None
	for temp_c in whl.find({'RAdeg':{'$gt':ir.ra.deg-max_sep, '$lt':ir.ra.deg+max_sep}, 'DEdeg':{'$gt':ir.dec.deg-max_sep, '$lt':ir.dec.deg+max_sep}, \
							'$or':[{'zspec':{'$gt':z-dz, '$lt':z+dz}}, {'zphot':{'$gt':z-dz, '$lt':z+dz}}]}):
		current_sep = ir.separation( coord.SkyCoord(temp_c['RAdeg'], temp_c['DEdeg'], unit=(u.deg,u.deg), frame='icrs') )
		if (current_sep < best_sep) and (current_sep < max_sep*u.deg):
			best_sep = current_sep
			cluster = temp_c
	
	return cluster

def get_redmapper(objID, ir, z, z_err, transverse, dz):
	'''
	Find the corresponding galaxy cluster in the redMaPPer catalog
	First check against member catalog, then check radially
	If multiple clusters match, choose the one with least angular separation
	'''
	
	# If the galaxy is too close, physical separations become too great on the sky
	# Restrict redshifts to at least 0.01
	if z < 0.01:
		return None
	
	member = rm_m.find_one({'ObjID':objID})
	if member is not None:
		return rm_c.find_one({'ID':member['ID']})
	
	# Maximum separation
	max_sep = float(transverse * u.Mpc / cosmo.angular_diameter_distance(z) * u.rad / u.deg)
	
	best_sep = np.inf
	cluster = None
	for temp_c in rm_c.find({'RAdeg':{'$gt':ir.ra.deg-max_sep, '$lt':ir.ra.deg+max_sep}, 'DEdeg':{'$gt':ir.dec.deg-max_sep, '$lt':ir.dec.deg+max_sep}, \
							 '$or':[{'zspec':{'$gt':z-dz, '$lt':z+dz}}, {'zlambda':{'$gt':z-dz, '$lt':z+dz}}]}):
		current_sep = ir.separation( coord.SkyCoord(temp_c['RAdeg'], temp_c['DEdeg'], unit=(u.deg,u.deg), frame='icrs') )
		if (current_sep < best_sep) and (current_sep < max_sep*u.deg):
			best_sep = current_sep
			cluster = temp_c
	
	return cluster

def get_amf(ir, z, z_err, transverse, dz):
	'''
	Find the corresponding galaxy cluster in the WHL15 catalog
	If multiple clusters match, choose the one with least angular separation
	'''
	
	# If the galaxy is too close, physical separations become too great on the sky
	# Restrict redshifts to at least 0.01
	if z < 0.01:
		return None
	
	# Maximum separation
	max_sep = float(transverse * u.Mpc / cosmo.angular_diameter_distance(z) * u.rad / u.deg)
	
	best_sep = np.inf
	cluster = None
	for temp_c in amf.find({'ra':{'$gt':ir.ra.deg-max_sep, '$lt':ir.ra.deg+max_sep}, 'dec':{'$gt':ir.dec.deg-max_sep, '$lt':ir.dec.deg+max_sep}, \
							'z':{'$gt':z-dz, '$lt':z+dz}}):
		current_sep = ir.separation( coord.SkyCoord(temp_c['ra'], temp_c['dec'], unit=(u.deg,u.deg), frame='icrs') )
		if (current_sep < best_sep) and (current_sep < max_sep*u.deg):
			best_sep = current_sep
			cluster = temp_c
	
	return cluster

def get_bending(source, peak_count):
	'''
	Calculate all the bending parameters that don't depend on the cluster
	'''
	
	assert (peak_count in [2, 3]), 'Not a valid morphology'
	subject = subjects.find_one({'zooniverse_id':source['zooniverse_id']})
	
	# Get pixel-to-WCS conversion
	fid = subject['metadata']['source']
	fits_loc = pathdict[fid]
	w = wcs.WCS(fits.getheader(fits_loc, 0))

	# Get the location of the source
	ir = coord.SkyCoord(source['SDSS']['ra'], source['SDSS']['dec'], unit=(u.deg,u.deg), frame='icrs') if 'SDSS' in source else coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
	
	ir_pos = w.wcs_world2pix(np.array([[ir.ra.deg,ir.dec.deg]]), 1)
	z, z_err = get_z(source)
	peaks = source['radio']['peaks']
	peak_pos = w.wcs_world2pix(np.array([ [peak['ra'],peak['dec']] for peak in peaks ]), 1)

	# Get image parameters for this source
	data = get_data(subject)
	contour_tree = get_contours(w, ir_pos, peak_pos, data, peak_count)
	peaks = get_global_peaks(w, peak_pos, peaks, contour_tree)
	if len(peaks) != 2:
		output("%s didn't have 2 tails" % source['zooniverse_id'])
		return
	
	contour_list = [child.path for child in contour_tree.children if any(child.contains(peak_pos))]
	
	# Using the 'contour' method
#	bending_angles = get_angles(w, ir, 'contour', contour_list)
#	tail_lengths_apparent = get_tail_lengths(w, ir, 'contour', contour_list)
#	tail_lengths_physical = []
#	for tail in tail_lengths_apparent:
#			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
#	ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
#	asymmetry = ratios[1]/ratios[0]
#	using_contour = {'tail_deg_0':tail_lengths_apparent[0], 'tail_deg_1':tail_lengths_apparent[1], 'size_deg':sum(tail_lengths_apparent), 'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
#	using_contour.update(bending_angles)
#	for key in using_contour.keys():
#		if type(using_contour[key]) is coord.angles.Angle:
#			using_contour[key] = using_contour[key].deg
	
	# Using the 'peak' method
	bending_angles = get_angles(w, ir, 'peak', peaks)
#	tail_lengths_apparent = get_tail_lengths(w, ir, 'peak', contour_list, peaks)
#	tail_lengths_physical = []
#	for tail in tail_lengths_apparent:
#			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
#	ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
#	asymmetry = ratios[1]/ratios[0]
#	using_peaks = {'tail_deg_0':tail_lengths_apparent[0], 'tail_deg_1':tail_lengths_apparent[1], 'size_deg':sum(tail_lengths_apparent), 'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
#	using_peaks.update(bending_angles)
#	for key in using_peaks.keys():
#		if type(using_peaks[key]) is coord.angles.Angle:
#			using_peaks[key] = using_peaks[key].deg
	
	morphology = 'double' if peak_count == 2 else 'triple'
	rgz = {'RGZ_id':source['catalog_id'], 'zooniverse_id':source['zooniverse_id'], 'ir_consensus':source['consensus']['ir_level'], 'radio_consensus':source['consensus']['radio_level'], 'peaks':peaks, 'components':source['radio']['components'], 'morphology':morphology, 'size_arcmin':source['radio']['max_angular_extent']/60., 'size_kpc':float((cosmo.angular_diameter_distance(z)*source['radio']['max_angular_extent']*np.pi/180.)/u.kpc)}
	
#	entry = {'RGZ':rgz, 'using_contour':using_contour, 'using_peaks':using_peaks}
	entry = {'RGZ':rgz, 'using_peaks':{}}
	
	if 'SDSS' in source:
		entry['SDSS'] = source['SDSS']
	if 'AllWISE' in source:
		entry['AllWISE'] = source['AllWISE']
	
	best_ra = entry['SDSS']['ra'] if 'SDSS' in entry else entry['AllWISE']['ra']
	best_dec = entry['SDSS']['dec'] if 'SDSS' in entry else entry['AllWISE']['dec']
	entry['best'] = {'ra':best_ra, 'dec':best_dec, 'redshift':z}
	
	return entry

def make_bent_sources():
	'''
	Generate a collection of radio sources with bending parameters that don't depend on the cluster
	Once generated, various matching schemes to the clusters can be tried efficiently
	'''
	
	# Determine which sources have already been processed
	completed = []
	if os.path.exists(completed_file):
		with open(completed_file, 'r') as f:
			lines = f.readlines()
			for line in lines:
				completed.append(int(line))
	
	z_range = [0.01, 0.8]
	# 'ignore_bending':False
	double_args = {'$and': [{'overedge':0, 'catalog_id':{'$nin':completed}}, \
							{'$or':  [{'radio.number_peaks':2, 'radio.number_components':1}, \
						 			  {'radio.number_components':2}]}, \
							{'$or':  [{'SDSS.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
									  {'SDSS.spec_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
									  {'AllWISE.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}] }]}
	triple_args = {'$and': [{'overedge':0, 'catalog_id':{'$nin':completed}}, \
							{'$or':  [{'radio.number_peaks':3, 'radio.number_components':{'$in':[1,2]}}, \
									  {'radio.number_components':3}]}, \
							{'$or':  [{'SDSS.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
									  {'SDSS.spec_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
									  {'AllWISE.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}] }]}
	
	# Find the bending parameters for each source that matches
	for args,peak_count,morphology in zip([double_args,triple_args], [2,3], ['double','triple']):
		count = bent_sources.find({'RGZ.morphology':morphology}).count()
		with open(completed_file, 'a') as f:
			for source in catalog.find(args).batch_size(50):
				entry = get_bending(source, peak_count)
				if entry is not None:
					count += 1
					output('%i %s' % (count, source['zooniverse_id']))
					bent_sources.insert(entry)
				print >> f, source['catalog_id']

def get_cluster_match(source):
	'''
	Given a source from RGZ, match it to a cluster and calculate the redshift-dependent bending parameters
	'''
	
	ir = coord.SkyCoord(source['SDSS']['ra'], source['SDSS']['dec'], unit=(u.deg,u.deg), frame='icrs') if 'SDSS' in source else \
		 coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
	z, z_err = get_z(source)
	
	# Match to cluster catalogs
	cluster_w = get_whl(ir, z, z_err, 15, 0.04*(1+z))
	whl_prop = {}
	if cluster_w is not None:
		c_pos = coord.SkyCoord(cluster_w['RAdeg'], cluster_w['DEdeg'], unit=(u.deg,u.deg), frame='icrs')
		c_sep_arc = c_pos.separation(ir)
		c_sep_mpc = float(cosmo.angular_diameter_distance(cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot'])/u.Mpc * c_sep_arc.to(u.rad)/u.rad)
		c_pos_angle = c_pos.position_angle(ir)
		if c_sep_mpc/cluster_w['r500'] < 0.01:
			pop = 'central'
		elif c_sep_mpc/cluster_w['r500'] >= 1.5:
			pop = 'outer'
		else:
			pop = 'intermediate'
		whl_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg, 'r/r500':c_sep_mpc/cluster_w['r500'], 'population':pop, 'M500':np.exp(1.08*np.log(cluster_w['RL*500'])-1.37), 'zbest':cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot']}
		for key in ['_id', 'N500', 'N500sp', 'RL*500', 'name', 'r500', 'zphot', 'zspec']:
			if key in cluster_w:
				whl_prop[key] = cluster_w[key]
	
	# Only continue if a cluster was matched
	if cluster_w is None: #and cluster_r is None and cluster_a is None:
		
		output("%s didn't match to a cluster" % source['RGZ']['zooniverse_id'])
		return
	
	else:
		
		z = cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot']
		
		# Using the 'contour' method
#		tail_lengths_apparent = [source['using_contour']['tail_deg_0']*u.deg, source['using_contour']['tail_deg_1']*u.deg]
#		tail_lengths_physical = []
#		for tail in tail_lengths_apparent:
#			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
#		using_contour = {'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc)}
#		source['using_contour'].update(using_contour)
		
		# Using the 'peak' method
#		tail_lengths_apparent = [source['using_peaks']['tail_deg_0']*u.deg, source['using_peaks']['tail_deg_1']*u.deg]
#		tail_lengths_physical = []
#		for tail in tail_lengths_apparent:
#			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
#		using_peaks = {'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc)}
#		source['using_peaks'].update(using_peaks)
		
		# Calculate the orientation angle
#		for cluster,prop in zip([cluster_w,cluster_r,cluster_a], [whl_prop,rm_prop,amf_prop]):
#			if cluster is not None:
		cluster, prop = cluster_w, whl_prop
#		for method,using in zip(['contour','peaks'], [source['using_contour'],source['using_peaks']]):
		for method,using in zip(['peaks'], [source['using_peaks']]):
			orientation = coord.angles.Angle(prop['position_angle'] - using['bisector'], unit=u.deg).wrap_at(360.*u.deg).deg
			if orientation > 180.:
				orientation = 360. - orientation
			prop['orientation_%s' % method] = orientation
			if orientation > 90.:
				orientation = 180. - orientation
			prop['orientation_folded'] = orientation
		
		# Compile all results
		if cluster_w is not None:
			source['WHL'] = whl_prop
#		if cluster_r is not None:
#			source['redMaPPer'] = rm_prop
#		if cluster_a is not None:
#			source['AMFDR9'] = amf_prop
		
		return source

def make_catalog():
	'''
	Find the sources matching the given search parameters and morphology and run the processing pipeline
	'''
	
	# Determine which sources have already been processed
	completed = []
	if os.path.exists(completed_file):
		with open(completed_file, 'r') as f:
			lines = f.readlines()
			for line in lines:
				completed.append(int(line))
	
	# Find the bending and cluster results for each source that matches
	for peak_count,morphology in zip([2,3], ['double','triple']):
		count = bending_15.find({'RGZ.morphology':morphology}).count()
		with open(completed_file, 'a') as f:
			for source in bent_sources.find({'RGZ.RGZ_id': {'$nin':completed}, 'RGZ.morphology':morphology}).batch_size(50):
				entry = get_cluster_match(source)
				if entry is not None:
					count += 1
					output('%i %s' % (count, source['RGZ']['zooniverse_id']))
					bending_15.insert(entry)
				print >> f, source['RGZ']['RGZ_id']
	
	# Apply the bending correction
	bending_correct(plot=False, update=True, methods=['using_peaks'])

def random_control():
	'''
	Generate a control sample by randomizing the positions and redshifts of all the sources meeting the radio morphology requirements
	'''
	
	# Get the original location values
	ras1, ras2, decs1, decs2 = [], [], [], []
	for source in bent_sources.find():
		ra = source['SDSS']['ra'] if 'SDSS' in source else source['AllWISE']['ra']
		dec = source['SDSS']['dec'] if 'SDSS' in source else source['AllWISE']['dec']
		if 90 < ra < 290:
			ras1.append(ra)
			decs1.append(dec)
		else:
			ras2.append(ra)
			decs2.append(dec)
	
	# Shuffle the location values
	np.random.shuffle(ras1)
	np.random.shuffle(ras2)
	np.random.shuffle(decs1)
	np.random.shuffle(decs2)
	loc = np.vstack([ np.append(ras1,ras2), np.append(decs1,decs2) ]).T
	np.random.shuffle(loc)
	
	# Find the bending and cluster results for each source that matches
	for ix, source in enumerate(bent_sources.find().batch_size(50)):
		
		# Assign the randomized values
		if 'SDSS' in source:
			source['SDSS'].update({'ra':loc[ix][0], 'dec':loc[ix][1]})
		else:
			source['SDSS'] = {'ra':loc[ix][0], 'dec':loc[ix][1]}
		
		entry = get_cluster_match(source)
		if entry is not None:
			output('%i %s' % (ix, source['RGZ']['zooniverse_id']))
			bending_control.insert(entry)

def output(string, fn=logging.info):
	'''
	Print a string to screen and the logfile
	'''
	fn(string)
	print string

def to_file(filename, coll, params={}):
	'''
	Print the bending collection to a csv file for analysis
	'''
	rgz_keys = ['RGZ_id', 'zooniverse_id', 'first_id', 'morphology', 'radio_consensus', 'ir_consensus', 'size_arcmin', 'size_kpc', 'solid_angle']
	best_keys = ['ra', 'ra_err', 'dec', 'dec_err', 'redshift']
	sdss_keys = ['ra', 'dec', 'objID', 'photo_redshift', 'photo_redshift_err', 'spec_redshift', 'spec_redshift_err', 'u', 'u_err', 'g', 'g_err', 'r', 'r_err', 'i', 'i_err', 'z', 'z_err', 'number_matches']
	wise_keys = ['ra', 'dec', 'designation', 'photo_redshift', 'w1mpro', 'w1sigmpro', 'w1snr', 'w2mpro', 'w2sigmpro', 'w2snr', 'w3mpro', 'w3sigmpro', 'w3snr', 'w4mpro', 'w4sigmpro', 'w4snr']
	whl_keys = ['name', 'ra', 'dec', 'zphot', 'zspec', 'N500', 'N500sp', 'RL*500', 'M500', 'r500', 'separation_deg', 'separation_Mpc', 'r/r500', 'population', 'position_angle', 'orientation_contour', 'orientation_peaks', 'orientation_folded', 'P', 'P500', 'aligned']
	rm_keys = ['name', 'ra', 'dec', 'zlambda', 'zspec', 'S', 'lambda', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	amf_keys = ['AMF_id', 'ra', 'dec', 'z', 'r200', 'richness', 'core_radius', 'concentration', 'likelihood', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	bending_keys = ['pos_angle_0', 'pos_angle_1', 'bending_angle', 'bending_corrected', 'bending_err', 'bending_excess', 'bisector', 'asymmetry'] #, 'tail_deg_0', 'tail_deg_1', 'size_deg', 'tail_kpc_0', 'tail_kpc_1', 'size_kpc', 'ratio_1', 'ratio_0']
	
	all_keys = [rgz_keys, best_keys, sdss_keys, wise_keys, whl_keys, bending_keys, bending_keys]
	dict_names = ['RGZ', 'best', 'SDSS', 'AllWISE', 'WHL', 'using_contour', 'using_peaks']
	
	success = 0
	with open(filename, 'w') as f:
		
		header = 'final_sample'
		for superkey, key_list in zip(dict_names, all_keys):
			for key in key_list:
				header += ',%s.%s' % (str(superkey), str(key))
		print >> f, header
		
		for entry in coll.find(params):
			try:
				row = str(entry['final_sample'])
				for superkey, key_list in zip(dict_names, all_keys):
					for key in key_list:
						if superkey in entry and key in entry[superkey]:
							if type(entry[superkey][key]) is long or type(entry[superkey][key]) is bson.int64.Int64:
								row += ",'%s'" % str(entry[superkey][key])
							else:
								row += ',%s' % str(entry[superkey][key])
						else:
							row += ',-99'
				print >> f, row
				success += 1
			except BaseException as e:
				output(e, logging.exception)
		
		output('%i/%i successfully printed to %s' % (success, coll.find(params).count(), filename))

def plot_running(x_param, y_param, coll=bending_15, morph=None, pop=None, bin_by=None, bin_count=0, logx=False, logy=False, bent_cut=0, align=None):
	'''
	Plot the running 1st, 2nd, and 3rd quartiles averaged over window data points
	x_param, y_param, and (optional) bin_by need to be 'category.key' to search for coll[category][key]
	morph can be 'double' or 'triple' to select only that morphology
	pop can be 'central', 'affected', or 'outer' to select only that population
	The axes can be plotted on a log scale by specifying logx and logy to True or False
	'''
	
	assert morph in [None, 'double', 'triple'], "morph must be 'double' or 'triple'"
	assert pop in [None, 'central', 'intermediate', 'outer', 'separate', 'non-BCG'], "pop must be 'central', 'intermediate', 'outer', 'separate', or 'non-BCG'"
	assert align in [None, 'radial', 'tangential'], "align must be 'radial' or 'tangential'"
	if bin_by is not None:
		assert type(bin_count) is int and bin_count>0, 'bin_count must be positive int'
	
	# Prepare parameters for search
	params = total_cuts.copy()
	x_param_list = x_param.split('.')
	y_param_list = y_param.split('.')
	for param in [x_param, y_param]:
		if param in params:
			params[param]['$exists'] = True
		else:
			params[param] = {'$exists':True}
	
	if morph is not None:
		params['RGZ.morphology'] = morph
	
	if align is not None:
		params['WHL.aligned'] = align
	
	if pop == 'separate':
		pop_list = ['intermediate', 'outer', 'central']
	elif pop == 'non-BCG':
		pop_list = ['intermediate', 'outer']
	else:
		pop_list = [pop]
	
	params['using_peaks.bending_corrected']['$gte'] = bent_cut
	
	# Open the plotting window
	fig, ax = plt.subplots()
	box = ax.get_position()
	
	# Plot trends for each population of interest
	needs_labels = True
	windows = set()
	if bin_by is not None:
		plot_params = {'RGZ.size_arcmin': {'name':"Size (')", 'fmt':'%.2f\n - %.2f', 'width':0.91}, \
					   'best.redshift':	  {'name':'Redshift', 'fmt':'%.2f\n - %.2f', 'width':0.89}, \
					   'WHL.M500':		  {'name':'Mass', 'fmt':'%.1f\n - %.1f', 'width':0.91}}
		bin_by_list = bin_by.split('.')
		bins = np.arange(bin_count+1) * 100. / bin_count
		vals = []
		for i in coll.find(params):
			vals.append(i[bin_by_list[0]][bin_by_list[1]])
		samples = np.percentile(vals, bins)
		ax.plot(0, 0, c='w', label=plot_params[bin_by]['name'])
		for pop2 in pop_list:
			if pop2 is not None:
				params['WHL.population'] = pop2
			for i in range(len(samples)-1):
				params[bin_by] = {'$gte':samples[i], '$lt':samples[i+1]}
				window_size, run_x_50, run_y_50 = get_trends(params, x_param, y_param, coll, True)
				windows.add(window_size)
				if needs_labels:
					ax.plot(run_x_50, run_y_50, label=plot_params[bin_by]['fmt'] % (samples[i], samples[i+1]), color='C%i'%i)
					ax.set_position([box.x0, box.y0, box.width * plot_params[bin_by]['width'], box.height])
				else:
					ax.plot(run_x_50, run_y_50, color='C%i'%i)
			needs_labels = False
		
	else:
		for pop2 in pop_list:
			if pop2 is not None:
				params['WHL.population'] = pop2
			window_size, run_x_50, run_y_25, run_y_50, run_y_75 = get_trends(params, x_param, y_param, coll, False)
			windows.add(window_size)
			if needs_labels:
				ax.plot(run_x_50, run_y_75, label='75%', color='C0')
				ax.plot(run_x_50, run_y_50, label='50%', color='C1')
				ax.plot(run_x_50, run_y_25, label='25%', color='C2')
				ax.set_position([box.x0, box.y0, box.width * 0.94, box.height])
				needs_labels = False
			else:
				ax.plot(run_x_50, run_y_75, color='C0')
				ax.plot(run_x_50, run_y_50, color='C1')
				ax.plot(run_x_50, run_y_25, color='C2')
	
	# Make the plot pretty
	ax.set_xlabel(x_param)
	ax.set_ylabel(y_param)
	windows_txt = str(min(windows))
	if len(windows) > 1:
		windows_txt += str(-1*max(windows))
	ax.set_title('%s %s%ssources (window size: %s)' % (morph.title() if type(morph) is str else 'All', align+' ' if type(align) is str else '', pop+' ' if type(pop) is str and pop!='separate' else '', windows_txt))
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	if logx:
		ax.set_xscale('log')
	if logy:
		ax.set_yscale('log')
	if y_param == 'using_peaks.bending_corrected':
		params['WHL.population'] = 'outer'
		bend = []
		for i in coll.find(params):
			bend.append(i['using_peaks']['bending_corrected'])
		ax.axhline(np.median(bend), ls='dotted', c='k')
	elif y_param == 'using_peaks.bending_excess':
		ax.axhline(0, ls='dotted', c='k')

def get_trends(params, x_param, y_param, coll, binned):
	x_param_list = x_param.split('.')
	y_param_list = y_param.split('.')
	x, y = [], []
	for i in coll.find(params).sort(x_param, 1):
		x.append(i[x_param_list[0]][x_param_list[1]])
		y.append(i[y_param_list[0]][y_param_list[1]])
	
	x = np.array(x)
	y = np.array(y)
	window_size = min(len(x)/10, 100)
	
	if 'WHL.population' in params and params['WHL.population'] == 'central':
		run_x_50 = [0.01, 0.011]
		run_y_25 = 2*[np.percentile(y, 25)]
		run_y_50 = 2*[np.percentile(y, 50)]
		run_y_75 = 2*[np.percentile(y, 75)]
	else:
		run_x_50 = np.array([np.percentile(x[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
		run_y_25 = np.array([np.percentile(y[ix:ix+window_size], 25) for ix in np.arange(len(y)-window_size+1)])
		run_y_50 = np.array([np.percentile(y[ix:ix+window_size], 50) for ix in np.arange(len(y)-window_size+1)])
		run_y_75 = np.array([np.percentile(y[ix:ix+window_size], 75) for ix in np.arange(len(y)-window_size+1)])
	
	if binned:
		return window_size, run_x_50, run_y_50
	else:
		return window_size, run_x_50, run_y_25, run_y_50, run_y_75

def get_params(param_list, coll=bending_15, morph=None, pop=None):
	'''
	Returns an array of data containing the values of the specified paramters for all sources in the sample
	'''
	
	assert morph in [None, 'double', 'triple'], "morph must be 'double' or 'triple'"
	assert pop in [None, 'central', 'intermediate', 'outer'], "pop must be 'central', 'intermediate', or 'outer'"
	
	params = total_cuts.copy()
	if morph is not None:
		params['RGZ.morphology'] = morph
	if pop is not None:
		params['WHL.population'] = pop
	
	data = np.zeros([len(param_list), coll.find(params).count()])
	for jx, gal in enumerate(coll.find(params)):
		for ix, param in enumerate(param_list):
			datum = gal[param[0]][param[1]]
			if datum=='double':
				datum = 2
			elif datum=='triple':
				datum = 3
			elif datum=='central':
				datum = 0
			elif datum=='intermediate':
				datum = 1
			elif datum=='outer':
				datum = 2
			data[ix,jx] = datum
	return data

def custom_exp(x, a, b):
	return a*np.exp(b*x)

def bending_correct(coll=bending_15, window_size=100, plot=False, update=False, methods=None):
	'''
	Apply a pre-determined correction for the angular size dependence to the bending angle
	'''
	
	if methods is None:
		methods = ['using_peaks', 'using_contour']
	elif type(methods) is str:
		methods = [methods]
	
	# Repeat for both morphologies and angle-measuring methods separately
	for morph in ['double', 'triple']:
		for method in methods:
			
			# Print progress
			print morph, method
			
			# Collect data from outer region
			sizes, angles, separations = [], [], []
			params = {'RGZ.size_arcmin':total_cuts['RGZ.size_arcmin'], 'RGZ.overedge':total_cuts['RGZ.overedge'], 'RGZ.radio_consensus':total_cuts['RGZ.radio_consensus'], method+'.bending_angle':total_cuts['using_peaks.bending_corrected'], 'RGZ.morphology':morph, 'WHL.population':'outer'}
			for gal in bending_15.find(params).sort('RGZ.size_arcmin', 1):
				sizes.append(gal['RGZ']['size_arcmin'])
				angles.append(gal[method]['bending_angle'])
				separations.append(gal['WHL']['r/r500'])
			
			sizes = np.array(sizes)
			angles = np.array(angles)
			separations = np.array(separations)
			
			# Find running trend
			sizes_running = np.array([np.median(sizes[ix:ix+window_size]) for ix in np.arange(len(sizes)-window_size+1)])
			angles_running = np.array([np.median(angles[ix:ix+window_size]) for ix in np.arange(len(angles)-window_size+1)])
			
			# Find best fit
			med_size = np.median(sizes)
			popt, pcov = curve_fit(custom_exp, sizes_running-med_size, angles_running)
			angles_best = custom_exp(sizes-med_size, *popt)
			angles_best_running = np.array([np.median(angles_best[ix:ix+window_size]) for ix in np.arange(len(angles_best)-window_size+1)])
			
			# Save fits for comparison
			if method == 'using_peaks':
				if morph == 'double':
					d_x = sizes
					d_popt = popt
					d_x0 = med_size
				else:
					t_x = sizes
					t_popt = popt
					t_x0 = med_size
			
			# Plot fits
			if plot:
				fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
				ax1.plot(sizes_running, angles_running, label='Median fit')
				ax1.plot(sizes_running, angles_best_running, label='%.1f*e^(%+.1f*(size-%.2f))\nR^2 = %.3f' % (popt[0], popt[1], med_size, r2_score(angles_running, angles_best_running)))
				ax1.legend(loc='upper right')
				ax1.set_ylabel(method+'.bending_angle')
				ax1.set_title('Best fit %s sources' % morph)
				ax2.plot(sizes_running, angles_running-angles_best_running, label='residual')
				ax2.axhline(0, color='k', lw=1)
				ax2.set_xlabel('RGZ.size_arcmin')
				ax2.set_ylabel('residual')
				plt.subplots_adjust(hspace=0.05)
			
			# Insert to Mongo
			if update:
				del params['WHL.population']
				for gal in coll.find():
					best_angle = custom_exp(gal['RGZ']['size_arcmin']-med_size, *popt)
					corr_angle = popt[0] * gal[method]['bending_angle'] / best_angle
					coll.update({'_id':gal['_id']}, {'$set': {method+'.bending_corrected':corr_angle}})
	
	# Plot comparison of fits
	if plot:
		x = np.sort(np.hstack([d_x, t_x]))
		y0 = custom_exp(x-d_x0, *d_popt)
		y1 = custom_exp(x-t_x0, *t_popt)
		fig, ax = plt.subplots(1)
		ax.plot(x, y0, label='Best fit doubles')
		ax.plot(x, y1, label='Best fit triples')
		ax.plot(1, 10, c='w', label='R^2 = %.3f' % r2_score(y0, y1))
		ax.legend(loc='upper right')
		ax.set_xlabel('RGZ.size_arcmin')
		ax.set_ylabel('using_peaks.bending_angle')
		ax.set_title('Best fit comparison')

def pca_analysis(param_space, coll=bending_15, morph=None, pop=None):
	
	assert param_space in ['bending', 'WHL']
	assert morph in [None, 'double', 'triple'], "morph must be 'double' or 'triple'"
	assert pop in [None, 'central', 'intermediate', 'outer', 'non-BCG'], "pop must be 'central', 'intermediate', 'outer', or 'non-BCG'"
	
	if param_space == 'bending':
		param_list = np.array([['using_peaks', 'bending_angle'], ['RGZ', 'size_arcmin'], ['RGZ', 'size_kpc'], ['WHL', 'r/r500'], ['WHL', 'M500'], ['WHL', 'orientation_folded'], ['best', 'redshift']])
		names = np.array(['log(bending)', 'size_arcmin', 'size_kpc', 'r/r500', 'log(M500)', 'orientation', 'redshift'])
	else:
		param_list = np.array([['WHL', 'RL*500'], ['WHL', 'r500'], ['WHL', 'M500'], ['WHL', 'zbest']])
		names = np.array(['log(RL*500)', 'r500', 'log(M500)', 'redshift'])
	
	# Get data
	if pop == 'non-BCG':
		data1 = get_params(param_list, coll=coll, morph=morph, pop='intermediate')
		data2 = get_params(param_list, coll=coll, morph=morph, pop='outer')
		data = np.hstack([data1, data2.T[data2[3]<=10].T])
	else:
		data = get_params(param_list, coll=coll, morph=morph, pop=pop)
	
	# Plot bending, mass and separation on log scales
	if param_space == 'bending':
		for i in [0, 4]:
			data[i] = np.log10(data[i])
		if pop is None:
			data[3] = np.log10(data[3])
	else:
		for i in [0, 2]:
			data[i] = np.log10(data[i])
	
	# Normalize the data
	mean = np.mean(data, 1).reshape(-1,1)*np.ones(data.shape)
	std = np.std(data, 1).reshape(-1,1)*np.ones(data.shape)
	data_normed = (data-mean)/std
	
	# Run PCA
	pca = PCA()
	output = pca.fit_transform(data_normed.T).T
	
	# Plot the features and PCA outputs
	if param_space == 'bending' and pop is None:
		names[3] = 'log(r/r500)'
	corner(data.T, labels=names)
	corner(output.T, labels=['PC%i'%i for i in range(len(names))])
	
	# What the PCs are comprised of, printed as html table
	output = np.hstack([names.reshape(-1,1), pca.components_])
	print '<table style="width: 100%%;" border="0" width="1887" height="127"><tbody>'
	for ix, row in enumerate(output.T):
		if ix == 0:
			out_str = '<tr><td>PC# (contribution to variance)</td><td>'
		else:
			out_str = '<tr><td>PC%i (%.3f)</td><td>' % (ix-1, pca.explained_variance_ratio_[ix-1])
		for val in row:
			try:
				fval = float(val)
				out_str += '%.2f</td><td>' % fval
			except ValueError:
				out_str += val + '</td><td>'
		print out_str[:-4] + '</tr>'
	print '</tbody></table>'

def make_corner_plot():
	# Get data
	param_list = np.array([['using_peaks', 'bending_angle'], ['RGZ', 'size_arcmin'], ['WHL', 'r/r500'], ['WHL', 'M500'], ['WHL', 'orientation_folded'], ['best', 'redshift']])
	data = get_params(param_list, coll=bending_15)

	# Discard BCGs
	#data = data[:, data[5]!=0][:5]

	# Check what outlier cuts need to be made
	labels = np.array(['Bending (deg)', 'Size (arcmin)', 'Separation\n(r/r500)', 'Cluster mass\n(10^14 M_sun)', 'Redshift'])
	corner(data.T, labels=labels)

	# Apply readability cuts
	bending_cut = data[0]<60
	separation_cut = data[2]<10
	mass_cut = data[3]<24
	data = data[:, np.logical_and.reduce((bending_cut, separation_cut, mass_cut))]
	corner(data.T, labels=labels)
	plt.tight_layout()

def sample_numbers():
	print 'Consensus cut (double, triple):', bent_sources.find({'RGZ.morphology':'double', 'RGZ.radio_consensus':{'$gte':.65}}).count(), bent_sources.find({'RGZ.morphology':'triple', 'RGZ.radio_consensus':{'$gte':.65}}).count()
	print 'Redshifts (SDSS, AllWISE):', bent_sources.find({'RGZ.radio_consensus':{'$gte':.65}, '$or':[{'SDSS.spec_redshift':{'$exists':True}}, {'SDSS.photo_redshift':{'$exists':True}}]}).count(), bent_sources.find({'RGZ.radio_consensus':{'$gte':.65}, 'SDSS.spec_redshift':{'$exists':False}, 'SDSS.photo_redshift':{'$exists':False}}).count()
	print 'Size/bending cut (double, triple):', bent_sources.find({'RGZ.morphology':'double', 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_angle':{'$lte':135.}}).count(), bent_sources.find({'RGZ.morphology':'triple', 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_angle':{'$lte':135.}}).count()
	print 'Matched (double, triple):', bending_15.find({'RGZ.morphology':'double', 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_angle':{'$lte':135.}}).count(), bending_15.find({'RGZ.morphology':'triple', 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_angle':{'$lte':135.}}).count()
	print 'Final (double, triple):', bending_15.find({'RGZ.morphology':'double', 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_corrected':{'$lte':135.}}).count(), bending_15.find({'RGZ.morphology':'triple', 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_corrected':{'$lte':135.}}).count()

def contamination():
	z = []
	dz = []
	with open('/home/garon/Documents/RGZdata/bending/dz_whl_with_z.csv', 'r') as f:
		r = csv.reader(f)
		r.next()
		for row in r:
			z.append(float(row[0]))
			dz.append(float(row[1]))
		z = np.array(z)
		dz = np.array(dz)
	
	dz_norm = dz/(1+z)
	mask = np.abs(dz_norm)<0.04
	comp = 0.982 # Completeness fraction using 0.04 cut
	cont = 0.215 # Contamination fraction using 0.04 cut
	
	f, ax = plt.subplots(1)
	rc('text', usetex=True)
	n, bins, patches = ax.hist(dz_norm[np.abs(dz_norm)<=.1], 30)
	ax.axvline(0.04, c='r', label=r'$\Delta z$ threshold')
	ax.axvline(-0.04, c='r')
	ax.axhline(np.median(n[(np.abs(bins)>0.04)[:-1]]), c='k', label='Background')
	ax.plot(0, 0, c='w', label='Completeness: %.3f\nContamination: %.3f' % (comp, cont))
	ax.legend()
	ax.set_xlim(-0.1, 0.1)
	ax.set_xlabel(r'$\Delta z / (1+z)$')
	ax.set_ylabel('Count')
	plt.title('Full RGZ-WHL cross-match')
	rc('text', usetex=False)
	
	sep, n, bins = r500_hist(20)
	mean_cont = cont * sum(n) / (np.pi*np.square(bins[-1])) # Number of contaminating sources per square r500
	bin_area = np.pi * np.array( [np.square(bins[0])] + [np.square(bins[i+1])-np.square(bins[i]) for i in np.arange(len(bins)-1)] )
	cont_per_bin = mean_cont * bin_area
	fractional_cont = cont_per_bin / n
	fc = interp1d(bins, fractional_cont, kind='linear')
	for i in 100*np.logspace(-4, 0, 20):
		xx = brentq(lambda x: fc(x)-i/100., bins[0], bins[-1])
		print '%.1f%% at x=%.2f' % (i,xx)
	
	f, ax = plt.subplots(1)
	density = n/bin_area
	ax.scatter(bins, density, label='Measured sources')
	a, b = np.polyfit(np.log10(bins[1:-3]), np.log10(density[1:-3]), 1)
	ax.plot(bins, pow(10, a*np.log10(bins)+b), c='k', label='Best fit, 0.01<r/r500<10')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.legend()
	ax.set_xlabel('Separation (r/r500)')
	ax.set_ylabel('Count / area of annulus')
	ax.set_title('Count density vs separation')
	ax.set_xlim([0.0075, 85])
	ax.set_ylim([0.00016, 4e6])
	
	f, (ax1, ax2) = plt.subplots(2, sharex=True)
	ax1.plot(bins, n, label='Measured sources')
	ax1.plot(bins, cont_per_bin, label='Contaminating sources')
	ax1.legend()
	ax1.set_xscale('log')
	ax1.set_yscale('log')#, subsy=[2,3,4,5,6,7,8,9,20,30,40,50,60,70,80,90])
	#ax1.set_yticks(10.**np.arange(-4,4))
	ax1.set_ylim(ymin=0.01)
	ax1.set_ylabel('Count')
	ax1.set_title('Contamination')
	ax2.plot(bins, fractional_cont)
	ax2.axhline(1, ls='dotted', c='k')
	ax2.set_yscale('log', subsy=[2,3,4,5,6,7,8,9,20,30,40,50,60,70,80,90])
	ax2.set_yticks(10.**np.arange(-7,3))
	ax2.set_xlabel('WHL.r/r500')
	ax2.set_ylabel('Fractional contamination')
	plt.tight_layout()

def r500_hist(bin_count=20):
	params = total_cuts.copy()
	del params['using_peaks.bending_corrected']
	sep = []
	for i in bending_15.find(params):
		sep.append(i['WHL']['r/r500'])
	
	sep = np.array(sep)
	min_sep = min(np.log(sep))
	sep = np.clip(sep, .01, None) # Combine everything less than 0.01 into one bin
	n, bins, patches = plt.hist(np.log(sep), bins=bin_count)
	bins0 = bins[0]
	bins[0] = min_sep
	n, bins, patches = plt.hist(sep, bins=np.exp(bins))
	plt.xscale('log')
	bins[0] = np.exp(bins0)
	return sep, n, (bins[:-1]+bins[1:])/2.

def orientation(bent_cutoff=None, folded=True, r_min=0.01, r_max=10):
	sep = []
	bend = []
	ori = []
	for i in bending_15.find(total_cuts):
		sep.append(i['WHL']['r/r500'])
		bend.append(i['using_peaks']['bending_corrected'])
		if folded:
			ori.append(i['WHL']['orientation_folded'])
		else:
			ori.append(i['WHL']['orientation_peaks'])
	
	sep = np.array(sep)
	bend = np.array(bend)
	ori = np.array(ori)
	
	near = np.logical_and(0.01<sep, sep<1.5)
	far = np.logical_and(1.5<sep, sep<10)
	seps = np.vstack([near, far])
	print sum(near), 'near sources,', sum(far), 'far sources'
	
	if bent_cutoff is None:
		bent_cutoff = pow(10, np.mean(np.log10(bend)) + np.std(np.log10(bend)))
	straight = bend<bent_cutoff
	bent = bend>=bent_cutoff
	bends = np.vstack([bent, straight])
	print sum(straight), 'straight sources,', sum(bent), 'bent sources'
	
#	f, ax = plt.subplots(2, 2, sharex='col', sharey='row')
#	plt.suptitle('Orientation counts by subsample')
#	for ix in range(2):
#		for jx in range(2):
#			mask = np.logical_and(bends[ix], seps[jx])
#			data = ori[mask]
#			ad = ad_test(data, stats.uniform(0,90))
#			print ad.pvalue, len(data)
#			ax[ix][jx].hist(data, bins=min(20, len(data)/10), label='p=%.3f'%ad.pvalue)
#			ax[ix][jx].legend(loc='upper center')
#	
#	ax[0][0].set_ylabel('Bending > %.0f deg' % bent_cutoff)
#	ax[1][0].set_ylabel('Bending < %.0f deg' % bent_cutoff)
#	ax[1][0].set_xlabel('0.01 < r/r500 < 1.5')
#	ax[1][1].set_xlabel('1.5 < r/r500 < 10.0')
#	plt.tight_layout()
#	plt.subplots_adjust(top=.92)
#	
#	f, ax = plt.subplots(1)
#	ax.set_title('Orientation counts')
#	data = ori[bent]
#	ad = ad_test(data, stats.uniform(0,90))
#	print ad.pvalue, len(data)
#	ax.hist(data, bins=min(20, len(data)/10), label='p=%.3f'%ad.pvalue)
#	ax.legend(loc='upper center')
#	ax.set_ylabel('Bending > %.0f deg' % bent_cutoff)
#	ax.set_xlabel('All separations')
#	plt.tight_layout()
	
	if folded:
		max_ori = 90.
	else:
		max_ori = 180.
	
	f, ax = plt.subplots(1)
	ax.set_title('Orientation distribution (%.2f<r/r500<%.2f; bending>=%.0f deg)' % (r_min, r_max, bent_cutoff))
	data = ori[np.logical_and(bent, np.logical_and(sep>r_min, sep<r_max))]
	ad = ad_test(data, stats.uniform(0,max_ori))
	print ad.pvalue, z_score(ad.pvalue), len(data)
	ax.hist(data, bins=6, label='n=%i; p=%.3f'%(len(data),ad.pvalue) )
	ax.legend(loc='upper center')
	ax.set_ylabel('Count')
	ax.set_xlabel('Orientation angle (deg)')
	plt.tight_layout()
	
	f, ax = plt.subplots(1)
	ax.set_title('Orientation distribution (%.2f<r/r500<%.2f; bending<%.0f deg)' % (r_min, r_max, bent_cutoff))
	data = ori[np.logical_and(straight, np.logical_and(sep>r_min, sep<r_max))]
	ad = ad_test(data, stats.uniform(0,max_ori))
	print ad.pvalue, z_score(ad.pvalue), len(data)
	ax.hist(data, bins=6, label='n=%i; p=%.3f'%(len(data),ad.pvalue) )
	ax.legend(loc='upper center')
	ax.set_ylabel('Count')
	ax.set_xlabel('Orientation angle (deg)')
	plt.tight_layout()
	
	y, yerr = [], []	
	for i in np.arange(int(max_ori/2.)):
		a = sum(data<i)
		b = sum(data>max_ori-i)
		y.append(a-b)
		yerr.append(np.sqrt(a+b))
	
	f, ax = plt.subplots(1)
	ax.errorbar(np.arange(int(max_ori/2.)), y, yerr=yerr)
	ax.axhline(0)
	ax.set_xlabel('Orientation cut')
	ax.set_ylabel('Count')
	ax.set_title('Excess of inward-moving sources')
	plt.tight_layout()
	
	mask = np.logical_and(bent, np.logical_and(sep>0.01, sep<15))
	n = sum(mask)
	bin_count = 6
	bins = int(1.*n/bin_count) * np.arange(bin_count+1)
	bins[-1] = -1
	sep_sort = np.sort(sep[mask])
	ori_sort = ori[mask][np.argsort(sep[mask])]
	x, diff, err = [], [], []
	for i,j in zip(bins[:-1],bins[1:]):
		inward = sum(ori_sort[i:j]<30)
		outward = sum(ori_sort[i:j]>60)
		x.append(np.median(sep_sort[i:j]))
		diff.append(inward-outward)
		err.append(np.sqrt(inward+outward))

	f, ax = plt.subplots(1)
	ax.errorbar(x, diff, err)
	ax.axhline(0, c='k', ls='dotted')
	ax.set_xscale('log')
	ax.set_xlabel('WHL.r/r500')
	ax.set_ylabel('Count')
	ax.set_title('Excess of radially-moving sources')

def size_dependence():
	params0, params1 = total_cuts.copy(), total_cuts.copy()
	params0['WHL.r/r500'] = {'$gte':0.01, '$lt':1.5}
	params1['WHL.r/r500'] = {'$gte':1.5, '$lt':10.}

	# Get trends
	window_size, run_x_50, run_y_25, size0, run_y_75 = get_trends(params0, 'WHL.r/r500', 'RGZ.size_arcmin', bending_15, False)
	sep0 = (run_x_50-run_x_50[0])/(run_x_50[-1]-run_x_50[0])
	window_size, run_x_50, run_y_25, size1, run_y_75 = get_trends(params1, 'WHL.r/r500', 'RGZ.size_arcmin', bending_15, False)
	sep1 = (run_x_50-run_x_50[0])/(run_x_50[-1]-run_x_50[0])

	# Downsample size1 to same length as size0
	mask = (np.linspace(0,1,len(size0)) * (len(size1)-1)).astype(int)

	# Get values
	size0, size1 = [], []
	for i in bending_15.find(params0):
		size0.append(i['RGZ']['size_arcmin'])

	for i in bending_15.find(params1):
		size1.append(i['RGZ']['size_arcmin'])

	# AD test the values
	print 'Different separations:', stats.anderson_ksamp([size0, size1])

	# Repeat for masses
	params0, params1 = total_cuts.copy(), total_cuts.copy()
	params0['WHL.M500'] = {'$lte':10}
	params1['WHL.M500'] = {'$gte':15}
	size0, size1 = [], []
	for i in bending_15.find(params0):
		size0.append(i['RGZ']['size_arcmin'])

	for i in bending_15.find(params1):
		size1.append(i['RGZ']['size_arcmin'])

	# AD test the values
	print 'Different masses:', stats.anderson_ksamp([size0, size1])

def trend_tests():
	
	# Get trends for mass bins
	params0, params1 = total_cuts.copy(), total_cuts.copy()
	params0['WHL.r/r500'] = {'$lte':7.8}
	params1['WHL.r/r500'] = {'$gte':7.8}
	window_size, run_x_50, run_y_25, size0, run_y_75 = get_trends(params0, 'WHL.r/r500', 'RGZ.size_arcmin', bending_15, False)
	sep0 = (run_x_50-run_x_50[0])/(run_x_50[-1]-run_x_50[0])
	window_size, run_x_50, run_y_25, size1, run_y_75 = get_trends(params1, 'WHL.r/r500', 'RGZ.size_arcmin', bending_15, False)
	sep1 = (run_x_50-run_x_50[0])/(run_x_50[-1]-run_x_50[0])
	print stats.mannwhitneyu(size0, size1, alternative='two-sided')

def get_errs(get_errs=False):
	first_err = np.sqrt(0.3**2+0.02**2) # https://arxiv.org/pdf/1501.01555.pdf
	for source in bending_15.find().batch_size(100):
		# Positional errors
		if get_errs:
			if 'SDSS' in source:
				sql = 'select raerr, decerr from photoprimary where objid=%i' % source['SDSS']['objID']
				df = SDSS_select(sql)
				ra_err = df['raerr'][0]
				dec_err = df['decerr'][0]
			else:
				ir_pos = coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
				table = Irsa.query_region(ir_pos, catalog='allwise_p3as_psd', radius=1.*u.arcsec)
				ra_err = table['sigra'][0]
				dec_err = table['sigdec'][0]
		else:
			ra_err = source['best']['ra_err']
			dec_err = source['best']['dec_err']
		pos_err = np.sqrt(ra_err**2 + dec_err**2 + first_err**2)
		
		# Morphology errors
		area = 0
		for comp in source['RGZ']['components']:
			area += comp['solid_angle']
		
		# Total errors
		size = 60.* source['RGZ']['size_arcmin']
		frac_pos_err = pos_err / size / 4. # fractional pos error
		morph_err = area / size**2 * 180. / np.pi # total morph error in deg
		frac_morph_err = morph_err / source['using_peaks']['bending_angle']
		total_err = np.sqrt(frac_pos_err**2 + frac_morph_err**2)
		bend_err = total_err * source['using_peaks']['bending_angle']
		
		bending_15.update({'_id':source['_id']}, {'$set': {'best.ra_err':ra_err, 'best.dec_err':dec_err, 'best.frac_positional_err':frac_pos_err, 'RGZ.solid_angle':area, 'RGZ.frac_morphology_err':frac_morph_err, 'using_peaks.bending_frac_err':total_err, 'using_peaks.bending_err':bend_err}})
	
	window, size, area_25, area_50, area_75 = get_trends(total_cuts, 'RGZ.size_arcmin', 'RGZ.solid_angle', bending_15, False)
	m, b = np.polyfit(size, area_50, 1)
	y = m*size+b
	#plt.plot(size, area_50, label='Running median')
	#plt.plot(size, y, label='%.0fx%+.0f\nR^2=%.3f' % (m, b, r2_score(area_50,y)))
	#plt.legend()

def rmsd(params, x_param, y_param, plot=True, coll=bending_15):
	x_param_list = x_param.split('.')
	y0_param_list = ['using_peaks', y_param]
	y1_param_list = ['using_contour', y_param]
	x, y0, y1 = [], [], []
	for i in coll.find(params).sort(x_param, 1):
		x.append(i[x_param_list[0]][x_param_list[1]])
		y0.append(i[y0_param_list[0]][y0_param_list[1]])
		y1.append(i[y1_param_list[0]][y1_param_list[1]])
	
	x = np.array(x)
	y0 = np.array(y0)
	y1 = np.array(y1)
	window_size = min(len(x)/10, 100)
	
	if 'WHL.population' in params and params['WHL.population'] == 'central':
		run_x = [0.01, 0.011]
		run_y0 = 2*[np.percentile(y0, 50)]
		run_y1 = 2*[np.percentile(y1, 50)]
		run_rmsd = 2*[np.sqrt(sum((y0-y1)**2)/len(x))]
	else:
		run_x = np.array([np.percentile(x[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
		run_y0 = np.array([np.percentile(y0[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
		run_y1 = np.array([np.percentile(y1[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
		run_rmsd = np.array([np.sqrt(sum((y0[ix:ix+window_size]-y1[ix:ix+window_size])**2)/window_size) for ix in np.arange(len(x)-window_size+1)])
	
	if plot:
		fig, (ax1, ax2) = plt.subplots(2, sharex=True)
		ax1.plot(run_x, run_y0, label=y0_param_list[0])
		ax1.plot(run_x, run_y1, label=y1_param_list[0])
		ax1.legend(loc='best')
		ax1.set_ylabel(y_param)
		ax1.set_title('%s comparison' % y_param)
		ax2.plot(run_x, run_rmsd)
		ax2.set_xlabel(x_param)
		ax2.set_ylabel('rms difference')
		plt.subplots_adjust(hspace=0.05)
		
		fig, ax = plt.subplots(1)
		ax.plot(run_x, run_y0, label='using_peaks.%s' % y_param)
		a, b = np.polyfit(run_x[run_x<1.05], run_rmsd[run_x<1.05], 1)
		ax.plot(run_x, a*run_x+b, label='rms difference linear fit')
		ax.legend(loc='best')
		ax.set_xlabel(x_param)
		ax.set_ylabel('degrees')
		ax.set_title('Bending error comparison')

	return window_size, run_x, run_y0, run_y1, run_rmsd

def z_score(p_vals, tail='two-sided'):
	assert tail in ['one-sided', 'two-sided'], 'tail must be one-sided or two-sided'
	if tail == 'one-sided':
		return stats.norm.ppf(1-p_vals)
	elif tail == 'two-sided':
		return stats.norm.ppf(1-p_vals/2.)

def mass_comp():
	
	params = total_cuts.copy()
	params['WHL.population'] = 'intermediate'
	
	mass = []
	bend = []
	for i in bending_15.find(params):
		mass.append(i['WHL']['M500'])
		bend.append(i['using_peaks']['bending_excess'])
	
	mass = np.array(mass)
	bend = np.array(bend)
	
	high = mass>12
	low = mass<12
	print sum(high), 'high mass sources,', sum(low), 'low mass sources'
	
	ad = stats.anderson_ksamp([bend[high], bend[low]])
	print 'p =', ad.significance_level
	print z_score(ad.significance_level), 'sigma'
	
	plt.hist(bend[high], normed=True, alpha=.8, bins=10*np.arange(15), label='M500 > 15*10^14 M_sun')
	plt.hist(bend[low], normed=True, alpha=.8, bins=10*np.arange(15), label='M500 < 10*10^14 M_sun')
	plt.plot(0, 0, color='w', label='p=%f'%ad.significance_level)
	plt.legend()
	plt.xlabel('using_peaks.bending_excess')
	plt.ylabel('Normalized count')
	plt.title('Excesss bending distribution by mass')

def bcg_comp():
	
	sep = []
	bend = []
	for i in bending_15.find(total_cuts):
		sep.append(i['WHL']['r/r500'])
		bend.append(i['using_peaks']['bending_excess'])
	
	sep = np.array(sep)
	bend = np.array(bend)
	
	bcg = sep<0.01
	outer = np.logical_and(sep>1.5, sep<10)
	print sum(bcg), 'BCGs,', sum(outer), 'outer sources'
	
	ad = stats.anderson_ksamp([bend[bcg], bend[outer]])
	print 'p =', ad.significance_level
	print z_score(ad.significance_level), 'sigma'
	
	plt.hist(bend[bcg], normed=True, alpha=.8, bins=10*np.arange(15), label='r/r500 < 0.01')
	plt.hist(bend[outer], normed=True, alpha=.8, bins=10*np.arange(15), label='1.5 < r/r500 < 10')
	plt.plot(0, 0, color='w', label='p=%.5f'%ad.significance_level)
	plt.legend()
	plt.xlabel('using_peaks.bending_excess')
	plt.ylabel('Normalized count')
	plt.title('Excess bending distribution by separation')

def asym_comp():

	params = total_cuts.copy()
	params['WHL.aligned'] = 'radial'

	sep = []
	asym = []
	p = []
	for i in bending_15.find(params):
		sep.append(i['WHL']['r/r500'])
		asym.append(i['using_contour']['asymmetry'])
		p.append(i['WHL']['P'])

	sep = np.array(sep)
	asym = np.array(asym)
	p = np.array(p)

	sep_i = sep[np.logical_and(sep>0.01, sep<1.5)]
	sep_o = sep[sep>1.5]
	inter = asym[np.logical_and(sep>0.01, sep<1.5)]
	outer = asym[sep>1.5]

	s = stats.spearmanr(sep_i, inter)

	print s, z_score(s.pvalue)
	
	s = stats.spearmanr(p[np.logical_and(sep>0.01, sep<1.5)], inter)
	
	print s, z_score(s.pvalue)

def excess_comp():

	params = total_cuts.copy()
	params['WHL.population'] = 'intermediate'

	sep = []
	bend = []
	for i in bending_15.find(params):
		sep.append(i['WHL']['r/r500'])
		bend.append(i['using_peaks']['bending_excess'])

	sep = np.array(sep)
	bend = np.array(bend)

	s = stats.spearmanr(sep, bend)

	print s, z_score(s.pvalue)

def excess(bending):
	diff = bending**2 - median_bend**2
	return np.sign(diff) * np.sqrt(np.abs(diff))

def get_bending_excess(coll=bending_15):
	for i in coll.find():
		bend = i['using_peaks']['bending_corrected']
		bending_15.update({'_id':i['_id']}, {'$set':{'using_peaks.bending_excess':excess(bend)}})

def get_asymmetry(coll=bending_15):
	for method in ['using_peaks', 'using_contour']:
		print method
		for i in coll.find():
			angle_c = coord.Angle(i['WHL']['position_angle']*u.deg)
			angle_0 = coord.Angle(i[method]['pos_angle_0']*u.deg)
			angle_1 = coord.Angle(i[method]['pos_angle_1']*u.deg)
			diff_0 = min( (angle_c-angle_0).wrap_at('360d'), (angle_0-angle_c).wrap_at('360d') )
			diff_1 = min( (angle_c-angle_1).wrap_at('360d'), (angle_1-angle_c).wrap_at('360d') )
			if diff_0 < diff_1:
				inner = i[method]['tail_deg_0']
				outer = i[method]['tail_deg_1']
			else:
				inner = i[method]['tail_deg_1']
				outer = i[method]['tail_deg_0']
			if (diff_0.degree<45 and diff_1.degree>135) or (diff_1.degree<45 and diff_0.degree>135):
				alignment = 'radial'
			elif (45<diff_0.degree<135) and (45<diff_1.degree<135):
				alignment = 'tangential'
			else:
				alignment = 'other'
			
			if method == 'using_contour':
				bending_15.update({'_id':i['_id']}, {'$set':{method+'.asymmetry':inner/outer, 'WHL.aligned':alignment}})
			else:
				bending_15.update({'_id':i['_id']}, {'$set':{method+'.asymmetry':inner/outer}})

def pressure_calculation():
	'''Calculates the pressure (in keV/cm^3) using Arnaud et al. 2010'''
	
	h70 = float(cosmo.H(0) / (70.*u.km/u.Mpc/u.s) )
	alphaP = 1./0.561 - 5./3.
	P0 = 8.403*pow(h70,-1.5)
	c500, gamma, alpha, beta = 1.177, 0.3081, 1.0510, 5.4905
	
	h = lambda z: float(cosmo.H(z) / cosmo.H(0))
	P500 = lambda z, M500: 1.65e-3 * pow(h(z),8./3.) * pow(M500*h70/3.,2./3.) * pow(h70,2.)
	pp = lambda x: P0 * pow(c500*x,-1*gamma) * pow(1+pow(c500*x,alpha),(gamma-beta)/alpha)
	alphaP_prime = lambda x: 0.10 - (alphaP+0.10) * pow(2*x,3) / (1+pow(2*x,3))
	P = lambda x, z, M500: P500(z,M500) * pp(x) * pow(M500*h70/3.,alphaP+alphaP_prime(x))
	
	for source in bending_15.find():
		x = source['WHL']['r/r500']
		z = source['WHL']['zbest']
		M500 = source['WHL']['M500']
		bending_15.update({'_id':source['_id']}, {'$set':{'WHL.P500':P500(z,M500), 'WHL.P':P(x,z,M500)}})

def get_fits(params, outd):
	import shutil
	
	if os.path.exists(outd):
		count = len(os.listdir(outd))
		cont = raw_input('%s exists and contains %i files; continue? (y/n) ' % (outd,count))
		if cont.lower() == 'n':
			return
	else:
		print 'Initializing', outd
		os.makedirs(outd)
	
	with open('%s/first_fits.txt' % rgz_path) as f:
		lines = f.readlines()
	
	pathdict = {}
	for l in lines:
		spl = l.split(' ')
		pathdict[spl[1].strip()] = '%s/rgz/raw_images/RGZ-full.%i/FIRST-IMGS/%s.fits' % (data_path, int(spl[0]), spl[1].strip())
	
	for source in bending_15.find(params).sort('WHL.r/r500', -1).limit(55):
		print source['WHL']['r/r500']
		zid = source['RGZ']['zooniverse_id']
		fid = catalog.find_one({'zooniverse_id':zid})['first_id']
		fits_loc = pathdict[fid]
		shutil.copy(fits_loc, outd+fid+'.fits')

def flag_dups(coll=bending_15, params={}):
	
	coll.update({}, {'$set': {'RGZ.duplicate_flag':0}}, multi=True)
	
	params['RGZ.duplicate_sources'] = {'$exists':True}
	for source in coll.find(params):
		dups = source['RGZ']['duplicate_sources']
		cid = source['RGZ']['RGZ_id']
		flag = int(cid != min(dups['exact_duplicate']))
		coll.update({'_id':source['_id']}, {'$set':{'RGZ.duplicate_flag':flag}})

if __name__ == '__main__':
	
	logging.basicConfig(filename='%s/bending.log' % rgz_path, level=logging.DEBUG, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	logging.captureWarnings(True)
	logging.info('Bending run from command line')
	
	# Run options
	calculate_bending = True
	match_to_clusters = True
	control = False
	
	# Generate the collection of bent sources
	if calculate_bending:
		done = False
		while not done:
			try:
				if not bent_sources.count():
					bent_sources.create_index('RGZ.RGZ_id', unique=True)
				output('Processing sources from RGZ')
				make_bent_sources()
				output('%i double sources processed' % bent_sources.find({'RGZ.morphology':'double'}).count())
				output('%i triple sources processed' % bent_sources.find({'RGZ.morphology':'triple'}).count())
				done = True
			except pymongo.errors.CursorNotFound as c:
				time.sleep(10)
				output('Cursor timed out; starting again.')
			except BaseException as e:
				logging.exception(e)
				raise
		with open(completed_file, 'w'): pass
	
	# Match the sources in bent_sources to the cluster catalogs
	if match_to_clusters:
		done = False
		while not done:
			try:
				if not bending_15.count():
					bending_15.create_index('RGZ.RGZ_id', unique=True)
				output('Matching sources to WHL')
				make_catalog()
				output('%i double sources matched to WHL' % bending_15.find({'RGZ.morphology':'double'}).count())
				output('%i triple sources matched to WHL' % bending_15.find({'RGZ.morphology':'triple'}).count())
				to_file('%s/csv/bending_catalog_15.csv' % rgz_path, bending_15)
				done = True
			except pymongo.errors.CursorNotFound as c:
				time.sleep(10)
				output('Cursor timed out; starting again.')
			except BaseException as e:
				logging.exception(e)
				raise
	
	# Generate a control sample by shuffling the positions of the sources
	if control:
		try:
			if not bending_control.count():
				bending_control.create_index('RGZ.RGZ_id', unique=True)
			output('Generating control sources')
			random_control()
			output('%i double sources matched to clusters' % bending_control.find({'RGZ.morphology':'double'}).count())
			output('%i triple sources matched to clusters' % bending_control.find({'RGZ.morphology':'double'}).count())
			to_file('%s/csv/bending_control_15.csv' % rgz_path, bending_control)
		except BaseException as e:
			logging.exception(e)
			raise
