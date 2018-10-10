import logging, urllib2, time, json, os
import pymongo, bson
import numpy as np
import StringIO, gzip
import csv
import pandas as pd
from astropy.io import fits
from astropy import wcs, coordinates as coord, units as u
from astropy.cosmology import Planck13 as cosmo
from scipy.optimize import brentq, curve_fit, leastsq
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
import matplotlib
from matplotlib.ticker import FormatStrFormatter
plt.ion()
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 14})

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
sdss_sample = db['sdss_sample']
xmatch = db['sdss_whl_xmatch']
distant_sources = db['distant_sources']

# Final sample cuts
total_cuts = {'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'RGZ.radio_consensus':{'$gte':0.65}, 'using_peaks.bending_corrected':{'$lte':135.}, 'best.redshift':{'$gte':0.02, '$lte':0.8}, 'RGZ.duplicate':0}

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
	
	# Remove the BCG source if it's a triple
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
	Returns the best redshift value and uncertainty for the source (only if SDSS)
	'''
	if 'SDSS' in source and 'spec_redshift' in source['SDSS']:
		return source['SDSS']['spec_redshift'], source['SDSS']['spec_redshift_err']
	elif 'SDSS' in source and 'photo_redshift' in source['SDSS']:
		return source['SDSS']['photo_redshift'], source['SDSS']['photo_redshift_err']
	else:
		return 0, 0
#	elif 'AllWISE' in source and 'photo_redshift' in source['AllWISE']:
#		return source['AllWISE']['photo_redshift'], 0

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
			pop = 'BCG'
		elif c_sep_mpc/cluster_w['r500'] >= 1.5:
			pop = 'outer'
		else:
			pop = 'inner'
		whl_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg, 'r/r500':c_sep_mpc/cluster_w['r500'], 'population':pop, 'zbest':cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot']}
		for key in ['_id', 'N500', 'N500sp', 'RL*500', 'name', 'r500', 'zphot', 'zspec', 'M500']:
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
	rgz_keys = ['RGZ_id', 'RGZ_name', 'zooniverse_id', 'first_id', 'morphology', 'radio_consensus', 'ir_consensus', 'size_arcmin', 'size_kpc', 'solid_angle', 'overedge']
	best_keys = ['ra', 'ra_err', 'dec', 'dec_err', 'redshift', 'redshift_err']
	sdss_keys = ['ra', 'dec', 'objID', 'photo_redshift', 'photo_redshift_err', 'spec_redshift', 'spec_redshift_err', 'u', 'u_err', 'g', 'g_err', 'r', 'r_err', 'i', 'i_err', 'z', 'z_err', 'morphological_class', 'spectral_class', 'number_matches']
	wise_keys = ['ra', 'dec', 'designation', 'photo_redshift', 'w1mpro', 'w1sigmpro', 'w1snr', 'w2mpro', 'w2sigmpro', 'w2snr', 'w3mpro', 'w3sigmpro', 'w3snr', 'w4mpro', 'w4sigmpro', 'w4snr']
	whl_keys = ['name', 'ra', 'dec', 'zphot', 'zspec', 'N500', 'N500sp', 'RL*500', 'M500', 'r500', 'separation_deg', 'separation_Mpc', 'r/r500', 'population', 'position_angle', 'orientation_contour', 'orientation_peaks', 'orientation_folded', 'P', 'P500', 'grad_P', 'alignment']
	rm_keys = ['name', 'ra', 'dec', 'zlambda', 'zspec', 'S', 'lambda', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	amf_keys = ['AMF_id', 'ra', 'dec', 'z', 'r200', 'richness', 'core_radius', 'concentration', 'likelihood', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	bending_keys = ['pos_angle_0', 'pos_angle_1', 'bending_angle', 'bending_corrected', 'bending_err', 'bending_excess', 'bisector', 'asymmetry'] #, 'tail_deg_0', 'tail_deg_1', 'size_deg', 'tail_kpc_0', 'tail_kpc_1', 'size_kpc', 'ratio_1', 'ratio_0']
	
	all_keys = [rgz_keys, best_keys, sdss_keys, whl_keys, bending_keys, bending_keys]
	dict_names = ['RGZ', 'best', 'SDSS', 'WHL', 'using_contour', 'using_peaks']
	
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

def plot_running(x_param, y_param, coll=bending_15, morph=None, pop=None, bin_by=None, bin_count=0, logx=False, logy=False, square=False, bent_cut=0, align=None, combined=False, title=True):
	'''
	Plot the running 1st, 2nd, and 3rd quartiles averaged over window data points
	x_param, y_param, and (optional) bin_by need to be 'category.key' to search for coll[category][key]
	morph can be 'double' or 'triple' to select only that morphology
	pop can be 'BCG', 'non-BCG', 'inner', or 'outer' to select only that population
	  When selecting 'non-BCG', combined=True will combine the 'inner' and 'outer' sources before smoothing
	align can be 'radial' or 'tangential' to select only that alignment
	The axes can be plotted on a log scale by specifying logx and logy to True or False
	  Also select square=True when plotting log-log plots
	Data can be binned with bin_by and bin_count
	bent_cut applies a minimum corrected bending angle
	'''
	
	assert morph in [None, 'double', 'triple'], "morph must be 'double' or 'triple'"
	assert pop in [None, 'BCG', 'inner', 'outer', 'separate', 'non-BCG'], "pop must be 'BCG', 'inner', 'outer', 'separate', or 'non-BCG'"
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
		params['WHL.alignment'] = align
	
	if pop == 'separate':
		pop_list = ['inner', 'outer', 'BCG']
	elif pop == 'non-BCG':
		if combined:
			pop_list = [{'$ne':'BCG'}]
		else:
			pop_list = ['inner', 'outer']
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
				window_size, run_x_50, run_y_50 = get_trends(params, x_param, y_param, coll, True, pop!='BCG')
				windows.add(window_size)
				if needs_labels:
					ax.plot(run_x_50, run_y_50, label=plot_params[bin_by]['fmt'] % (samples[i], samples[i+1]), color='C%i'%i)
					ax.set_position([box.x0, box.y0, box.width * plot_params[bin_by]['width'], box.height])
				else:
					ax.plot(run_x_50, run_y_50, color='C%i'%i)
			needs_labels = False
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		
	else:
		for pop2 in pop_list:
			if pop2 is not None:
				params['WHL.population'] = pop2
			window_size, run_x_50, run_y_25, run_y_50, run_y_75 = get_trends(params, x_param, y_param, coll, False, pop!='BCG')
			windows.add(window_size)
			if logx:
				run_x_50 = np.log10(run_x_50)
			if logy:
				run_y_25, run_y_50, run_y_75 = np.log10(run_y_25), np.log10(run_y_50), np.log10(run_y_75)
			if needs_labels:
				ax.plot(run_x_50, run_y_50, label='50\%', color='C0')
				ax.fill_between(run_x_50, run_y_25, run_y_75, color='C0', alpha=.5)
				needs_labels = False
			else:
				ax.plot(run_x_50, run_y_50, color='C0')
				ax.fill_between(run_x_50, run_y_25, run_y_75, color='C0', alpha=.5)
	
	# Make the plot pretty
	ax.set_xlabel(get_label(x_param, logx))
	ax.set_ylabel(get_label(y_param, logy))
	windows_txt = str(min(windows))
	if len(windows) > 1:
		windows_txt += str(-1*max(windows))
	titletxt = '%s%s%ssources (window size: %s)' % (morph+' ' if type(morph) is str else '', align+' ' if type(align) is str else '', pop+' ' if type(pop) is str and pop!='separate' else '', windows_txt)
	if title:
		ax.set_title('All '+titletxt if titletxt[0:7]=='sources' else titletxt[0].upper()+titletxt[1:])
	#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	if square:
		ax.set_aspect('equal', adjustable='box')
	if y_param == 'using_peaks.bending_corrected':
		params['WHL.population'] = 'outer'
		bend = []
		for i in coll.find(params):
			bend.append(i['using_peaks']['bending_corrected'])
		if logy:
			ax.axhline(np.log10(np.median(bend)), ls=':', c='k')
		else:
			ax.axhline(np.median(bend), ls=':', c='k')
	elif y_param == 'using_peaks.bending_excess':
		ax.axhline(0, ls=':', c='k')
	plt.tight_layout()

def get_trends(params, x_param, y_param, coll, binned, combine_bcg=True):
	x_param_list = x_param.split('.')
	y_param_list = y_param.split('.')
	x, y = [], []
	for i in coll.find(params).sort(x_param, 1):
		x.append(i[x_param_list[0]][x_param_list[1]])
		y.append(i[y_param_list[0]][y_param_list[1]])
	
	x = np.array(x)
	y = np.array(y)
	window_size = min(len(x)/10, 100)
	
	if 'WHL.population' in params and params['WHL.population'] == 'BCG' and combine_bcg:
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

def get_label(param, log=False):
	if param == 'WHL.r/r500':
		label = 'Separation ($r_{500}$)'
	elif param == 'WHL.M500':
		label = 'Cluster mass ($10^{14}~$M$_\odot$)'
	elif param == 'WHL.P':
		label = 'ICM pressure (keV cm$^{-3}$)'
	elif 'bending_angle' in param:
		label = 'Bending angle (deg)'
	elif 'bending_corrected' in param:
		label = 'Corrected bending angle (deg)'
	elif 'bending_excess' in param:
		label = 'Excess bending angle (deg)'
	elif 'asymmetry' in param:
		label = 'Asymmetry'
	elif param == 'RGZ.size_arcmin':
		label = 'Size (armin)'
	elif param == 'RGZ.size_kpc':
		label = 'Size (kpc)'
	elif param == 'WHL.grad_P':
		label = 'Pressure gradient (keV cm$^{-3}$ kpc$^{-1}$)'
	else:
		label = param.replace('_', '\_')
	if log:
		label = '$\log_{10}$ (' + label.replace('(', '[').replace(')', ']') + ')'
	return label

def get_params(param_list, coll=bending_15, morph=None, pop=None):
	'''
	Returns an array of data containing the values of the specified paramters for all sources in the sample
	'''
	
	assert morph in [None, 'double', 'triple'], "morph must be 'double' or 'triple'"
	assert pop in [None, 'BCG', 'inner', 'outer'], "pop must be 'BCG', 'inner', or 'outer'"
	
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
			elif datum=='BCG':
				datum = 0
			elif datum=='inner':
				datum = 1
			elif datum=='outer':
				datum = 2
			data[ix,jx] = datum
	return data

def custom_exp(x, a, b):
	return a*np.exp(b*x)

def bending_correct(coll=bending_15, window_size=100, plot=False, update=False, methods=None, comp_err=False):
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
			params = total_cuts.copy()
			del params['using_peaks.bending_corrected']
			params['RGZ.morphology'] = morph
			params[method+'.bending_angle'] = total_cuts['using_peaks.bending_corrected']
			params['WHL.population'] = 'outer'
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
			perr = np.sqrt(np.diagonal(pcov))
			angles_best = custom_exp(sizes-med_size, *popt)
			angles_best_running = np.array([np.median(angles_best[ix:ix+window_size]) for ix in np.arange(len(angles_best)-window_size+1)])
		
			# Save fits for comparison
			if method == 'using_peaks':
				if morph == 'double':
					d_x = sizes_running
					d_y = angles_running
					d_popt = popt
					d_err = perr
					d_x0 = med_size
				else:
					t_x = sizes_running
					t_y = angles_running
					t_popt = popt
					t_err = perr
					t_x0 = med_size
			
			# Plot fits
			if plot:
				fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
				ax1.plot(sizes_running, angles_running, label='Running median')
				ax1.plot(sizes_running, angles_best_running, ls='--', label='$%.1fe^{%+.1f(x-%.2f)}$\n$R^2 = %.3f$' % (popt[0], popt[1], med_size, r2_score(angles_running, angles_best_running)))
				ax1.legend(loc='upper right')
				ax1.set_ylabel('Bending angle (deg)')#('$\\mathrm{%s.bending_angle}$' % method).replace('_', '\_'))
				#ax1.set_title('Best fit %s sources' % morph)
				ax2.plot(sizes_running, angles_running-angles_best_running, label='residual')
				ax2.axhline(0, color='k', lw=1)
				ax2.set_xlabel('Size (arcmin)')
				ax2.set_ylabel('Residual')
				ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
				ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
				fig.tight_layout()
				plt.subplots_adjust(hspace=0.05)
	
			# Insert to Mongo
			if update:
				for gal in coll.find():
					best_angle = custom_exp(gal['RGZ']['size_arcmin']-med_size, *popt)
					corr_angle = popt[0] * gal[method]['bending_angle'] / best_angle
					coll.update({'_id':gal['_id']}, {'$set': {method+'.bending_corrected':corr_angle}})
				get_bending_excess()
	
	# Plot comparison of fits
	if plot:
		fig, ax = plt.subplots(1)
		ax.plot(d_x, d_y, label='Double sources')
		d_rms = np.sqrt(sum((d_y-custom_exp(d_x-d_x0, *d_popt))**2)/len(d_x))
		ax.fill_between(d_x, custom_exp(d_x-d_x0, *d_popt)+d_rms, custom_exp(d_x-d_x0, *d_popt)-d_rms, color='C0', alpha=.6)#, label='Best fit doubles')
		ax.plot(t_x, t_y, ls='--', label='Triple sources')
		t_rms = np.sqrt(sum((t_y-custom_exp(t_x-t_x0, *t_popt))**2)/len(t_x))
		ax.fill_between(t_x, custom_exp(t_x-t_x0, *t_popt)+t_rms, custom_exp(t_x-t_x0, *t_popt)-t_rms, color='C1', alpha=.6)#, label='Best fit triples')
		ax.legend(loc='upper right')
		ax.set_xlabel('Size (arcmin)')
		ax.set_ylabel('Bending angle (deg)')
		fig.tight_layout()

def pca_analysis(param_space, coll=bending_15, morph=None, pop=None):
	
	assert param_space in ['bending', 'WHL']
	assert morph in [None, 'double', 'triple'], "morph must be 'double' or 'triple'"
	assert pop in [None, 'BCG', 'inner', 'outer', 'non-BCG'], "pop must be 'BCG', 'inner', 'outer', or 'non-BCG'"
	
	if param_space == 'bending':
		param_list = np.array([['using_peaks', 'bending_angle'], ['RGZ', 'size_arcmin'], ['RGZ', 'size_kpc'], ['WHL', 'r/r500'], ['WHL', 'M500'], ['WHL', 'orientation_folded'], ['best', 'redshift']])
		names = np.array(['log(bending)', 'size_arcmin', 'size_kpc', 'r/r500', 'log(M500)', 'orientation', 'redshift'])
	else:
		param_list = np.array([['WHL', 'RL*500'], ['WHL', 'r500'], ['WHL', 'M500'], ['WHL', 'zbest']])
		names = np.array(['log(RL*500)', 'r500', 'log(M500)', 'redshift'])
	
	# Get data
	if pop == 'non-BCG':
		data1 = get_params(param_list, coll=coll, morph=morph, pop='inner')
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
	
	param_list = [{'RGZ.morphology':'double', 'RGZ.duplicate':0, 'best.redshift':{'$exists':True}}, \
		{'RGZ.morphology':'double', 'RGZ.duplicate':0, 'best.redshift':{'$exists':True}}, \
		{'RGZ.morphology':'double', 'RGZ.duplicate':0, 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_angle':{'$lte':135.}, 'best.redshift':{'$gte':0.02, '$lte':0.8}}, \
		{'RGZ.morphology':'double', 'RGZ.duplicate':0, 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_corrected':{'$lte':135.}, 'best.redshift':{'$gte':0.02, '$lte':0.8}}, \
		{'RGZ.morphology':'double', 'RGZ.duplicate':0, 'RGZ.radio_consensus':{'$gte':.65}, 'RGZ.size_arcmin':{'$lte':1.5}, 'RGZ.overedge':0, 'using_peaks.bending_corrected':{'$lte':135.}, 'best.redshift':{'$gte':0.02, '$lte':0.8}}]
	coll_list = [bent_sources, bending_15, bending_15, bent_sources, bending_15]
	label_list = ['No cuts', 'Matched', 'Initial cuts', 'Final (all)', 'Final (matched)']
	
	for params, coll, label in zip(param_list, coll_list, label_list):
		double = coll.find(params).count()
		params['RGZ.morphology'] = 'triple'
		triple = coll.find(params).count()
		print label, '(double, triple, total):', double, triple, double+triple

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
	n, bins, patches = ax.hist(dz_norm[np.abs(dz_norm)<=.1], 30)
	ax.axvline(0.04, c='r', lw=2, label=r'$\Delta z$ threshold')
	ax.axvline(-0.04, c='r', lw=2)
	ax.axhline(np.median(n[(np.abs(bins)>0.04)[:-1]]), c='k', lw=2, ls=':', label='Background')
	#ax.plot(0, 0, c='w', label='Completeness: %.3f\nContamination: %.3f' % (comp, cont))
	ax.legend()
	ax.set_xlim(-0.1, 0.1)
	ax.set_xlabel('$\Delta z / (1+z)$')
	ax.set_ylabel('Count')
	#plt.title('Full RGZ-WH15 cross-match')
	plt.tight_layout()

	sep, n, bins = r500_hist(bending_15, total_cuts, 20)
	mean_cont = cont * sum(n) / (np.pi*np.square(bins[-1])) # Number of contaminating sources per square r500
	bin_area = np.pi * np.array( [np.square(bins[0])] + [np.square(bins[i+1])-np.square(bins[i]) for i in np.arange(len(bins)-1)] )
	cont_per_bin = mean_cont * bin_area
	fractional_cont = cont_per_bin / n
	fc = interp1d(bins, fractional_cont, kind='slinear')
	print '%.3f%% at x=%.2f' % (fc(1.5), 1.5)
	#for i in np.logspace(-1, 2, 20):
	#	xx = brentq(lambda x: fc(x)-i/100., bins[0], bins[-1])
	#	print '%.1f%% at x=%.2f' % (i,xx)

	'''f, ax = plt.subplots(1)
	density = n/bin_area
	err = np.sqrt(n)/bin_area
	ax.errorbar(bins, density, yerr=err, fmt='o', ms=4, label='Source density')
	logx = np.log10(bins[1:-3])
	logy = np.log10(density[1:-3])
	logyerr = err[1:-3] / density[1:-3]
	fitfunc = lambda p, x: p[0]*x + p[1]
	errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
	out = leastsq(errfunc, [-1,1], args=(logx, logy, logyerr), full_output=1)
	index, amp = out[0]
	index_err = np.sqrt(out[1][1][1])
	amp_err = np.sqrt(out[1][0][0])*amp
	ax.plot(bins, pow(10, index*np.log10(bins)+amp), c='k', label='$\sim(r/r_{500})^{%.2f\pm%.2f}$'%(index,index_err))
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.legend(bbox_to_anchor=(1, 1))
	ax.set_xlabel('Separation ($r_{500}$)')
	ax.set_ylabel('Count / area of annulus')
	#ax.set_title('Count density vs. separation')
	ax.set_aspect('equal', adjustable='box')
	plt.tight_layout()'''

	f, (ax1, ax2) = plt.subplots(2, sharex=True)
	ax1.plot(np.log10(bins), np.log10(n), label='Observed')
	ax1.plot(np.log10(bins), np.log10(cont_per_bin), ls='--', label='Contamination')
	ax1.legend(loc='lower right')
	ax1.set_ylim(-2.24, 3.24)
	ax1.set_yticks(np.arange(-2,4))
	ax1.set_ylabel('$\log_{10}$ (Count)')
	#ax1.set_title('Contamination')
	ax2.plot(np.log10(bins), np.log10(fractional_cont), c='C2', label='Contamination\nfraction')
	ax2.legend(loc='lower right')
	ax2.axhline(0, ls=':', c='k')
	ax2.set_xlabel(get_label('WHL.r/r500', True))
	ax2.set_ylabel('$\log_{10}$ (Fraction)')
	plt.tight_layout()

def r500_hist(coll, params, bin_count=20):
	fig, ax = plt.subplots(1)
	sep = []
	for i in coll.find(params):
		sep.append(i['WHL']['r/r500'])
	sep = np.array(sep)
	min_sep = min(np.log10(sep))
	sep = np.clip(sep, .01, None) # Combine everything less than 0.01 into one bin
	n, bins, patches = ax.hist(np.log10(sep), bins=bin_count)
	bins0 = bins[0]
	bins[0] = min_sep
	n, bins, patches = ax.hist(sep, bins=pow(10,bins))
	ax.set_xscale('log')
	bins[0] = pow(10,bins0)
	return sep, n, (bins[:-1]+bins[1:])/2.

def orientation(bent_cutoff=None, folded=True, r_min=0.01, r_max=10):
	if bent_cutoff is None:
		bent_cutoff = get_bent_cut()
	
	print 'Bending excess cutoff: %.2f' % bent_cutoff
	
	sep = []
	bend = []
	ori = []
	for i in bending_15.find(total_cuts):
		sep.append(i['WHL']['r/r500'])
		bend.append(i['using_peaks']['bending_excess'])
		if folded:
			ori.append(i['WHL']['orientation_folded'])
		else:
			ori.append(i['WHL']['orientation_peaks'])

	sep = np.array(sep)
	bend = np.array(bend)
	ori = np.array(ori)

	inner = np.logical_and(0.01<sep, sep<1.5)
	outer = np.logical_and(1.5<sep, sep<10)
	seps = np.vstack([inner, outer])
	print sum(inner), 'inner sources,', sum(outer), 'outer sources'

	straight = bend<bent_cutoff
	bent = bend>=bent_cutoff
	bends = np.vstack([bent, straight])
	print sum(straight), 'straight sources,', sum(bent), 'bent sources'
	
	if folded:
		max_ori = 90.
	else:
		max_ori = 180.
	
	f, ax = plt.subplots(1)
	#ax.set_title('Orientation distribution ($%.2f<r/r_{500}<%.2f$; $\Delta\\theta<%.0f$ deg)' % (r_min, r_max, excess_cutoff))
	data = ori[np.logical_and(straight, np.logical_and(sep>r_min, sep<r_max))]
	ad = ad_test(data, stats.uniform(0,max_ori))
	print ad.pvalue, z_score(ad.pvalue), len(data)
	ax.hist(data, bins=6, fill=False, hatch='//', label='$n=%i;~p=%.3f$'%(len(data), ad.pvalue) )
	ax.set_ylim(0, 358)
	ax.legend(loc='upper center')
	ax.set_ylabel('Count')
	ax.set_xlabel('Orientation angle (deg)')
	plt.tight_layout()

	f, ax = plt.subplots(1)
	#ax.set_title('Orientation distribution ($%.2f<r/r_{500}<%.2f$; $\Delta\\theta>%.0f$ deg)' % (r_min, r_max, excess_cutoff))
	data = ori[np.logical_and(bent, np.logical_and(sep>r_min, sep<r_max))]
	if folded:
		ad = ad_test(data, stats.uniform(0,max_ori))
		print ad.pvalue, z_score(ad.pvalue), len(data)
	else:
		towards = data[data<90]
		away = data[data>90]
		sig = stats.anderson_ksamp([towards, max_ori-away]).significance_level
		print 'Symmetric: p=%.4f' % sig
	
	n, bins, _ = ax.hist(data, bins=6, fill=False, hatch='//', label='$n=%i;~p=%.3f$'%(len(data), ad.pvalue if folded else sig) )
	ax.set_ylim(0, 75)
	ax.legend(loc='upper center')
	ax.set_ylabel('Count')
	ax.set_xlabel('Orientation angle (deg)')
	plt.tight_layout()
	
	unfolded, sep = [], []
	params = total_cuts.copy()
	params['using_peaks.bending_excess'] = {'$gte': bent_cutoff}
	params['WHL.r/r500'] = {'$gte': r_min, '$lte': r_max}
	for source in bending_15.find(params):
		unfolded.append(source['using_peaks']['bisector'] - source['WHL']['position_angle'])
		sep.append(source['WHL']['r/r500'])

	sep = np.array(sep)
	unfolded = np.array(unfolded)
	tangential = np.sin(unfolded*np.pi/180)
	radial = np.cos(unfolded*np.pi/180)
	beta = 1 - np.var(tangential) / np.var(radial)
	print 'beta', beta

	return None
	
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
	ax.set_xlabel('Separation ($r_{500}$)')
	ax.set_ylabel('Count')
	ax.set_title('Excess of radially-moving sources')

def orientation_test():
	excess_cutoff = get_bent_cut()
	r_min, r_max = 0.01, 10

	sep, bend, folded, unfolded = [], [], [], []
	for i in bending_15.find(total_cuts):
		sep.append(i['WHL']['r/r500'])
		bend.append(i['using_peaks']['bending_excess'])
		folded.append(i['WHL']['orientation_folded'])
		unfolded.append(i['WHL']['orientation_peaks'])

	sep = np.array(sep)
	bend = np.array(bend)
	folded = np.array(folded)[np.logical_and(bend>excess_cutoff, np.logical_and(sep>r_min, sep<r_max))]
	unfolded = np.array(unfolded)[np.logical_and(bend>excess_cutoff, np.logical_and(sep>r_min, sep<r_max))]

	n_f, bins_f, _ = plt.hist(folded, bins=6, alpha=.8, normed=True, fill=False, edgecolor='r', label='Folded')
	plt.figure()
	n_u, bins_u, _ = plt.hist(unfolded, bins=12, alpha=.8, normed=True, fill=False, hatch='//', label='Unfolded')

	# model 1: count goes to 0 at 180 deg
	class quad(stats.rv_continuous):
		def _argcheck(self, a, b, c):
			return np.isfinite(a) and np.isfinite(b) and np.isfinite(c)
		def _pdf(self, x, a, b, c):
			x0 = 180.
			norm = a*x0**3/3. + b*x0**2/2. + c*x0
			if type(x) is float:
				return max(a*x**2 + b*x + c, 0) / norm
			elif type(x) is np.ndarray:
				return np.max([a*x**2 + b*x + c, np.zeros(len(x))], axis=0) / norm
			else:
				raise TypeError('Got %s instead' % str(type(x)))

	a, b, c = n_f[0], n_f[-1]/2., 0
	popt = np.polyfit([0,90,180], [a,b,c], 2)
	case1 = quad(a=0, b=180, shapes='a, b, c')(*popt)
	ad1 = ad_test(unfolded, case1)
	print 'Model 1: p=%.2g (%.2f sigma)' % (ad1.pvalue, z_score(ad1.pvalue))

	# model 2: count plateaus at 90 deg
	class piecewise_lin(stats.rv_continuous):
		def _pdf(self, x, a, b, c):
			norm = 45.*a + 135.*c
			m = (c-a) / 90.
			if type(x) is float:
				return (c + m*min([x-90, 0])) / norm
			elif type(x) is np.ndarray:
				return (c + m*np.min([x-90, np.zeros(len(x))], axis=0)) / norm
			else:
				raise TypeError('Got %s instead' % str(type(x)))

	a, b, c = n_f[0]-n_f[-1]/2., n_f[-1]/2., n_f[-1]/2.
	case2 = piecewise_lin(a=0, b=180, shapes='a, b, c')(a, b, c)
	ad2 = ad_test(unfolded, case2)
	print 'Model 1: p=%.2g (%.2f sigma)' % (ad2.pvalue, z_score(ad2.pvalue))
	
	x = np.arange(181)
	plt.plot(x, case1.pdf(x), c='C0', label='Model 1')
	plt.plot(x, case2.pdf(x), c='C1', ls=':', label='Model 2')
	plt.legend()
	plt.ylabel('Normalized count')
	plt.xlabel('Orientation angle (deg)')

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

def rmsd(params=total_cuts.copy(), x_param='RGZ.size_arcmin', y_param='bending_angle', plot=True, coll=bending_15):
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
	
	if 'WHL.population' in params and params['WHL.population'] == 'BCG':
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
		ax1.plot(run_x, run_y0, label='Peak method')
		ax1.plot(run_x, run_y1, label='Contour method')
		ax1.legend(loc='best')
		ax1.set_ylabel(get_label(y_param))
		#ax1.set_title('%s comparison' % y_param)
		ax2.plot(run_x, run_rmsd)
		ax2.set_xlabel(get_label(x_param))
		ax2.set_ylabel('RMS difference')
		fig.tight_layout()
		
		fig, ax = plt.subplots(1)
		ax.plot(run_x, run_y0, label='Bending angle (peak method)')
		#a, b = np.polyfit(run_x[run_x<1.05], run_rmsd[run_x<1.05], 1)
		#ax.plot(run_x, a*run_x+b, label='rms difference linear fit')
		ax.plot(run_x, run_rmsd, label='RMS difference', ls='--')
		ax.legend(loc='best')
		ax.set_xlabel(get_label(x_param))
		ax.set_ylabel('Angle (deg)')
		#ax.set_title('Bending error comparison')
		fig.tight_layout()
	
	return window_size, run_x, run_y0, run_y1, run_rmsd

def rmsd_debug():
	params = total_cuts.copy()
	params['using_peaks.bending_angle'] = params['using_peaks.bending_corrected']
	del params['using_peaks.bending_corrected']
	
	x_param = 'RGZ.size_arcmin'
	y_param = 'bending_angle'
	x_param_list = x_param.split('.')
	y0_param_list = ['using_peaks', y_param]
	y1_param_list = ['using_contour', y_param]
	x, y0, y1, zid = [], [], [], []
	for i in bending_15.find(params).sort(x_param, 1):
		x.append(i[x_param_list[0]][x_param_list[1]])
		y0.append(i[y0_param_list[0]][y0_param_list[1]])
		y1.append(i[y1_param_list[0]][y1_param_list[1]])
		zid.append(i['RGZ']['zooniverse_id'])
	
	x = np.array(x)
	y0 = np.array(y0)
	y1 = np.array(y1)
	window_size = min(len(x)/10, 100)
	print 'Original sample:', len(x)
	
	run_x = np.array([np.percentile(x[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_y0 = np.array([np.percentile(y0[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_y1 = np.array([np.percentile(y1[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_rmsd = np.array([np.sqrt(sum((y0[ix:ix+window_size]-y1[ix:ix+window_size])**2)/window_size) for ix in np.arange(len(x)-window_size+1)])
	
	fig, ax = plt.subplots(1)
	plt.scatter(x, y0, s=1, label='using peaks')
	plt.scatter(x, y1, s=1, label='using contour')
	plt.plot(run_x, run_rmsd, c='k', label='rmsd')
	plt.xlabel(get_label(x_param))
	plt.ylabel(get_label('%s.%s' % tuple(y0_param_list)))
	
	'''outlier = np.logical_and(np.logical_and(x>1.17, x<1.19), np.logical_and(y1>105, y1<109))
	x = x[np.logical_not(outlier)]
	y1 = y1[np.logical_not(outlier)]
	print np.where(outlier, zid, False)[outlier][0]
	
	run_x = np.array([np.percentile(x[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_y0 = np.array([np.percentile(y0[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_y1 = np.array([np.percentile(y1[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_rmsd = np.array([np.sqrt(sum((y0[ix:ix+window_size]-y1[ix:ix+window_size])**2)/window_size) for ix in np.arange(len(x)-window_size+1)])
	plt.plot(run_x, run_rmsd, c='g', label='outlier\nremoved')
	
	ell = matplotlib.patches.Ellipse(xy=(1.18,107), width=0.06, height=13, fill=False, color='r', label='outlier')
	ax.add_artist(ell)'''
	
	logdy = np.log10(np.abs(y0-y1))
	mask = logdy < 3*np.std(logdy)
	outliers = np.logical_not(mask)
	plt.scatter(x[outliers], y0[outliers], c='g', label='outliers', s=1)
	plt.scatter(x[outliers], y1[outliers], c='g', s=1)
	x = x[mask]
	y0 = y0[mask]
	y1 = y1[mask]
	print 'Outliers:', sum(outliers)
	
	run_x = np.array([np.percentile(x[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_y0 = np.array([np.percentile(y0[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_y1 = np.array([np.percentile(y1[ix:ix+window_size], 50) for ix in np.arange(len(x)-window_size+1)])
	run_rmsd = np.array([np.sqrt(sum((y0[ix:ix+window_size]-y1[ix:ix+window_size])**2)/window_size) for ix in np.arange(len(x)-window_size+1)])
	plt.plot(run_x, run_rmsd, c='g', label='$<3\sigma$')
	
	plt.legend()

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
	
    modified_z_score = 0.6745 * diff / med_abs_deviation
	
    return modified_z_score > thresh

def z_score(p_vals, tail='two-sided'):
	assert tail in ['one-sided', 'two-sided'], 'tail must be one-sided or two-sided'
	if tail == 'one-sided':
		return stats.norm.ppf(1-p_vals)
	elif tail == 'two-sided':
		return stats.norm.ppf(1-p_vals/2.)

def mass_comp(region):
	
	params = total_cuts.copy()
	params['WHL.population'] = region

	mass = []
	bend = []
	for i in bending_15.find(params):
		mass.append(i['WHL']['M500'])
		bend.append(i['using_peaks']['bending_excess'])

	mass = np.array(mass)
	bend = np.array(bend)
	
	rho = stats.spearmanr(mass, bend)
	print 'p = ', rho.pvalue
	print z_score(rho.pvalue), 'sigma'
	
	w, run_mass, run_bend = get_trends(params, 'WHL.M500', 'using_peaks.bending_excess', bending_15, True, False)
	popt = np.polyfit(run_mass, run_bend, 1)
	lower, upper = 5, 20
	print 'Rise of %.2f deg between %i and %i x 10^14 M_sun' % (popt[0]*(upper-lower), lower, upper)
	
	return None
	
	high = mass>12
	low = mass<12
	print sum(high), 'high mass sources,', sum(low), 'low mass sources'
	
	print 'bending difference:', np.median(bend[high])-np.median(bend[low])
	
	ad = stats.anderson_ksamp([bend[high], bend[low]])
	print 'p =', ad.significance_level
	print z_score(ad.significance_level), 'sigma'
	
	plt.hist(bend[high], normed=True, alpha=.8, bins=11*(np.arange(13)-1), label='$M_{500} > 1.2\\times 10^{15}~M_\\odot}$')
	plt.hist(bend[low], normed=True, alpha=.8, bins=11*(np.arange(13)-1), label='$M_{500} < 1.2\\times 10^{15}~M_\\odot}$')
	plt.plot(0, 0, color='w', label='$p=%f$'%ad.significance_level)
	plt.legend()
	plt.xlabel('Excess bending angle (deg)')
	plt.ylabel('Normalized count')

def mass_pop_comp():

	bent_cut = get_bent_cut()
	
	bent = []
	straight = []
	for i in bending_15.find(total_cuts):
		if i['using_peaks']['bending_excess']>bent_cut:
			bent.append(i['WHL']['M500'])
		else:
			straight.append(i['WHL']['M500'])
	
	ad = stats.anderson_ksamp([bent, straight])
	print 'p =', ad.significance_level
	print z_score(ad.significance_level), 'sigma'

def bcg_comp():
	
	sep = []
	bend = []
	for i in bending_15.find(total_cuts.copy()):
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
	
	plt.hist(bend[bcg], normed=True, alpha=.8, bins=10*np.arange(15), label='$r/r_{500} < 0.01$')
	plt.hist(bend[outer], normed=True, alpha=.8, bins=10*np.arange(15), label='$1.5 < r/r_{500} < 10$')
	plt.plot(0, 0, color='w', label='p=%.5f'%ad.significance_level)
	plt.legend()
	plt.xlabel('Excess bending angle (deg)')
	plt.ylabel('Normalized count')

def asym_comp(align='radial', update_cut=None):

	if update_cut is not None:
		get_asymmetry(cut)

	params = total_cuts.copy()
	params['WHL.alignment'] = align
	params['WHL.population'] = 'inner'
	print 'Alignment:', align

	sep = []
	asym = []
	p = []
	grad_p = []
	for i in bending_15.find(params):
		sep.append(i['WHL']['r/r500'])
		asym.append(i['using_contour']['asymmetry'])
		p.append(i['WHL']['P'])
		grad_p.append(i['WHL']['grad_P'])

	sep = np.array(sep)
	asym = np.array(asym)
	p = np.array(p)
	grad_p = np.array(grad_p)

	print 'n:', len(sep)
	
	popt, pcov = curve_fit(lambda x, a, b: a*x+b, sep, asym)
	perr = np.sqrt(np.diag(pcov))
	print 'slope:', popt[0], '+-', perr[0]

	s = stats.spearmanr(sep, asym)
	print 'r/r500:', s.pvalue, z_score(s.pvalue)

	s = stats.spearmanr(p, asym)
	#s = stats.spearmanr(p[p>5e-4)], asym[p>5e-4)])
	print 'P:', s.pvalue, z_score(s.pvalue)
	
	s = stats.spearmanr(grad_p, asym)
	print 'grad_P:', s.pvalue, z_score(s.pvalue)

def excess_comp():

	params = total_cuts.copy()
	params['WHL.population'] = 'inner'

	sep = []
	bend = []
	for i in bending_15.find(params):
		sep.append(i['WHL']['r/r500'])
		bend.append(i['using_peaks']['bending_excess'])

	sep = np.array(sep)
	bend = np.array(bend)

	s = stats.spearmanr(sep, bend)

	print s, z_score(s.pvalue)

def bend_corr_comp():

	params = total_cuts.copy()
	params['WHL.population'] = 'outer'
	
	morph = []
	bend = []
	size = []
	for source in bending_15.find(params):
		morph.append(2 if source['RGZ']['morphology']=='double' else 3)
		bend.append(source['using_peaks']['bending_angle'])
		size.append(source['RGZ']['size_arcmin'])
	
	morph = np.array(morph)
	bend = np.array(bend)
	size = np.array(size)
	
	ad = stats.anderson_ksamp([bend[morph==2], bend[morph==3]])
	print 'p =', ad.significance_level
	print z_score(ad.significance_level), 'sigma'
	
	plt.hist(bend[morph==2], normed=True, alpha=0.8, label='double', bins=20)
	plt.hist(bend[morph==3], normed=True, alpha=0.8, label='triple', bins=20)
	plt.xlabel(get_label('using_peaks.bending_angle'))
	plt.ylabel('Normalized count')
	plt.legend(loc='upper right')
	
	min_size = max(min(size[morph==2]), min(size[morph==3]))
	max_size = min(max(size[morph==2]), max(size[morph==3]))
	x = np.random.uniform(min_size, max_size, 1000)
	f1 = lambda x: 9.7*np.exp(-1.5*(x-0.54))
	f2 = lambda x: 7.6*np.exp(-1.3*(x-0.77))

def excess(bending, baseline):
	diff = bending**2 - baseline**2
	return np.sign(diff) * np.sqrt(np.abs(diff))

def get_median_bend(params=total_cuts.copy()):
	params['WHL.population'] = 'outer'
	bend = []
	for i in bending_15.find(params):
		bend.append(i['using_peaks']['bending_corrected'])
	return np.median(bend)

def get_bending_excess(params=total_cuts.copy()):
	median_bend = get_median_bend(params)
	for i in bending_15.find(params):
		bend = i['using_peaks']['bending_corrected']
		bending_15.update({'_id':i['_id']}, {'$set':{'using_peaks.bending_excess':excess(bend, median_bend)}})
	
	for i in bent_sources.find(params):
		bend = i['using_peaks']['bending_corrected']
		bent_sources.update({'_id':i['_id']}, {'$set':{'using_peaks.bending_excess':excess(bend, median_bend)}})

def get_asymmetry(cut, coll=bending_15, params=total_cuts.copy()):
	for method in ['using_peaks', 'using_contour']:
		print method
		if 'WHL' in coll.find_one(params):
			for i in coll.find(params):
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
				if (diff_0.degree<cut and diff_1.degree>180-cut) or (diff_1.degree<cut and diff_0.degree>180-cut):
					alignment = 'radial'
				elif cut<diff_0.degree and diff_0.degree<180-cut and cut<diff_1.degree and diff_1.degree<180-cut:
					alignment = 'tangential'
				else:
					alignment = 'other'
			
				if method == 'using_contour':
					coll.update({'_id':i['_id']}, {'$set':{method+'.asymmetry':inner/outer, 'WHL.alignment':alignment}})
				else:
					coll.update({'_id':i['_id']}, {'$set':{method+'.asymmetry':inner/outer}})
		else:
			for i in coll.find(params):
				if np.random.randint(0,2):
					inner = i[method]['tail_deg_0']
					outer = i[method]['tail_deg_1']
				else:
					inner = i[method]['tail_deg_1']
					outer = i[method]['tail_deg_0']
				coll.update({'_id':i['_id']}, {'$set':{method+'.asymmetry':inner/outer}})

def pressure_calc():
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

		theta0 = np.abs(source['using_peaks']['pos_angle_0'] - source['WHL']['position_angle'])
		theta1 = np.abs(source['using_peaks']['pos_angle_1'] - source['WHL']['position_angle'])
		if theta0>90:
			theta0 = 180 - theta0
		if theta1>90:
			theta1 = 180 - theta1
		dx0 = source['using_contour']['tail_kpc_0'] * np.cos(theta0) / 1000. / source['WHL']['r500']
		dx1 = source['using_contour']['tail_kpc_1'] * np.cos(theta1) / 1000. / source['WHL']['r500']
		dP = np.abs( (P(np.abs(x+dx0),z,M500) - P(np.abs(x+dx1),z,M500)) / (dx0 - dx1) )
		
		bending_15.update({'_id':source['_id']}, {'$set':{'WHL.P500':P500(z,M500), 'WHL.P':P(x,z,M500), 'WHL.grad_P':dP}})

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

def make_sdss_sample():
	
	if not sdss_sample.count():
		sdss_sample.create_index('bestObjID', unique=True)
	else:
		print '%i entries already in catalog' % sdss_sample.count()
	
	params = total_cuts.copy()
	params['best.ra'] = {'$gt':100, '$lt':275}

	zs, ras, decs = [], [], []
	for source in bending_15.find(params):
		zs.append(source['best']['redshift'])
		ras.append(source['best']['ra'])
		decs.append(source['best']['dec'])

	zs = np.array(zs)
	ras = np.array(ras)
	decs = np.array(decs)

	for x in [zs, ras, decs]:
		x.sort()

	zs_cdf = np.arange(len(zs))/float(len(zs)-1)
	zs_cdf_inv = interp1d(zs_cdf, zs)

	ras_cdf = np.arange(len(ras))/float(len(ras)-1)
	ras_cdf_inv = interp1d(ras_cdf, ras)

	decs_cdf = np.arange(len(decs))/float(len(decs)-1)
	decs_cdf_inv = interp1d(decs_cdf, decs)

	for i in range(20000):

		probs = np.random.uniform(0,1,3)
		z = zs_cdf_inv(probs[0])
		ra = ras_cdf_inv(probs[1])
		dec = decs_cdf_inv(probs[2])

		query = '''select top 1 s.bestObjID, s.ra, s.dec, s.z, s.zErr,
			   		g.cModelMag_u, g.cModelMag_g, g.cModelMag_r, g.cModelMag_i, g.cModelMag_z
			   from SpecObj as s
			       join GalaxyTag as g on s.bestObjID=g.objID
			   where s.class="GALAXY" and s.zWarning=0 and (s.z between %f and %f) and (s.ra between %f and %f) and 
				(s.dec between %f and %f)''' % (0.99*z, 1.01*z, 0.99*ra, 1.01*ra, 0.99*dec, 1.01*dec)
		df = SDSS_select(query)
		if len(df):
			entry = {}
			for key in df:
				entry[key] = df[key][0]
			try:
				sdss_sample.insert(entry)
			except pymongo.errors.DuplicateKeyError as e:
				print e
				pass

		print '%i/%i' % (sdss_sample.count(), i+1)

def sdss_patch():
	for s in sdss_sample.find({'cModelMag_u':{'$exists':False}}).batch_size(50):
		query = '''select cModelMag_u, cModelMag_g, cModelMag_r, cModelMag_i, cModelMag_z
				   from GalaxyTag as g
				   where objID=%s''' % s['bestObjID']
		df = SDSS_select(query)
		update = {}
		try:
			for key in df:
				update[key] = df[key][0]
			sdss_sample.update({'_id':s['_id']}, {'$set':update})
		except IndexError:
			pass
	
	for s in sdss_sample.find():
		xmatch.update({'SDSS._id':s['_id']}, {'$set':{'SDSS':s}})

def plot_sdss_sample():

	params = total_cuts.copy()
	params['best.ra'] = {'$gt':100, '$lt':275}
	params['SDSS.i'] = {'$exists':True}

	zs, ras, decs, mags = [], [], [], []
	for source in bending_15.find(params):
		zs.append(source['best']['redshift'])
		ras.append(source['best']['ra'])
		decs.append(source['best']['dec'])
		mags.append(source['SDSS']['i'])

	z_r, ra_r, dec_r, mag_r = [], [], [], []
	for source in sdss_sample.find({'cModelMag_i':{'$exists':True}}):
		z_r.append(source['z'])
		ra_r.append(source['ra'])
		dec_r.append(source['dec'])
		mag_r.append(source['cModelMag_i'])

	fig, ax = plt.subplots(1)
	ax.hist(zs, bins=15, alpha=.8, normed=True, label='Bending sample')
	ax.hist(z_r, bins=15, alpha=.8, normed=True, label='SDSS sample')
	ax.legend()
	ax.set_xlabel('z')
	ax.set_ylabel('Normalized count')
	ax.set_title('Redshift distribution (northern region)')
	
	fig, ax = plt.subplots(1)
	ax.hist(mags, bins=15, alpha=.8, normed=True, label='Bending sample')
	ax.hist(mag_r, bins=15, alpha=.8, normed=True, label='SDSS sample')
	ax.legend()
	ax.set_xlabel('i-band (mag)')
	ax.set_ylabel('Normalized count')
	ax.set_title('i-band magnitude distribution (northern region)')

	fig, ax = plt.subplots(1)
	ax.scatter(ras, decs, s=1, alpha=.8, label='Bending sample')
	ax.scatter(ra_r, dec_r, s=1, alpha=.8, label='SDSS sample')
	ax.legend()
	ax.set_xlabel('RA (deg)')
	ax.set_ylabel('Dec (deg)')
	ax.set_title('Skymap of sources (northern region)')

def sdss_xmatch():

	count = 0
	for source in sdss_sample.find().batch_size(50):
		count += 1
		ir = coord.SkyCoord(source['ra'], source['dec'], unit=(u.deg,u.deg), frame='icrs')
		cluster_w = get_whl(ir, source['z'], source['zErr'], 15, 0.04*(1+source['z']))
		if cluster_w is not None:
			whl_prop = {}
			c_pos = coord.SkyCoord(cluster_w['RAdeg'], cluster_w['DEdeg'], unit=(u.deg,u.deg), frame='icrs')
			c_sep_arc = c_pos.separation(ir)
			zbest = cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot']
			c_sep_mpc = float(cosmo.angular_diameter_distance(zbest)/u.Mpc * c_sep_arc.to(u.rad)/u.rad)
			c_pos_angle = c_pos.position_angle(ir)
			r = c_sep_mpc/cluster_w['r500']
			if r < 0.01:
				pop = 'BCG'
			elif r >= 1.5:
				pop = 'outer'
			else:
				pop = 'inner'
			whl_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'r/r500':r, 'population':pop, 'zbest':zbest}
			for key in ['_id', 'N500', 'N500sp', 'RL*500', 'name', 'r500', 'zphot', 'zspec', 'M500']:
				if key in cluster_w:
					whl_prop[key] = cluster_w[key]
			entry = {'SDSS':source, 'WHL':whl_prop}
			print '%i/%i' % (xmatch.count(), count)
			xmatch.insert(entry)

def sdss_density():
	params = total_cuts.copy()
	#params['SDSS.spec_redshift'] = {'$exists':True}
	
	f1, (ax1, ax2) = plt.subplots(1, 2, sharey=True, subplot_kw=dict(adjustable='datalim', aspect='equal'))
	ax = f1.add_subplot(111, frameon=False)
	ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
	f2, ax0 = plt.subplots(1)

	sep, n, bins = r500_hist(bending_15, params, 20)
	bin_area = np.pi * np.array( [np.square(bins[0])] + [np.square(bins[i+1])-np.square(bins[i]) for i in np.arange(len(bins)-1)] )
	density = n/bin_area
	err = np.sqrt(n)/bin_area
	ax1.errorbar(np.log10(bins), np.log10(density), yerr=np.vstack((np.log10(density+err) - np.log10(density), np.log10(density) - np.log10(density-err))), fmt='o', ms=4, label='RGZ')
	logx = np.log10(bins[1:-3])
	logy = np.log10(density[1:-3])
	logyerr = err[1:-3] / density[1:-3]
	if sum(np.isnan(logyerr)):
		logx = logx[np.logical_not(np.isnan(logyerr))]
		logy = logy[np.logical_not(np.isnan(logyerr))]
		logyerr = logyerr[np.logical_not(np.isnan(logyerr))]
	
	fitfunc = lambda p, x: p[0]*x + p[1]
	errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
	out = leastsq(errfunc, [-1,1], args=(logx, logy, logyerr), full_output=1)
	index, amp = out[0]
	index_err = np.sqrt(out[1][1][1])
	amp_err = np.sqrt(out[1][0][0])*amp
	ax1.plot(np.log10(bins), index*np.log10(bins)+amp, c='k', label='$\\alpha = $\n$  %.2f\pm%.2f$'%(index,index_err))
	err_max = np.maximum( (index+index_err)*np.log10(bins)+(amp+amp_err), (index-index_err)*np.log10(bins)+(amp+amp_err) )
	err_min = np.minimum( (index+index_err)*np.log10(bins)+(amp-amp_err), (index-index_err)*np.log10(bins)+(amp-amp_err) )
	#ax1.fill_between(np.log10(bins), err_min, err_max, color='k', lw=0, alpha=.5, label='$\\alpha = $\n$ %.2f\pm%.2f$'%(index,index_err))
	ax1.legend(loc='upper right')
	
	ax0.errorbar(np.log10(bins), np.log10(density), yerr=np.vstack((np.log10(density+err) - np.log10(density), np.log10(density) - np.log10(density-err))), fmt='o', ms=4, label='RGZ')
	ax0.plot(np.log10(bins), index*np.log10(bins)+amp, c='k', label='$\\alpha = $\n$  %.2f\pm%.2f$'%(index,index_err))
	
	r_bcgs = n[0]
	r_cluster = sum(n[bins<10])
	r_bcg_ratio = 1.*r_bcgs/r_cluster
	r_bcg_err = r_bcg_ratio * (1./np.sqrt(r_bcgs) + 1./np.sqrt(r_cluster))
	print 'RGZ BCGs: %i/%i (%.1f pm %.1f)%%' % (r_bcgs, r_cluster, 100*r_bcg_ratio, 100*r_bcg_err)
	r_within_r500 = sum(n[bins<1][1:])
	r_cluster_no_bcgs = sum(n[bins<10][1:])
	r_within_ratio = 1.*r_within_r500/r_cluster_no_bcgs
	r_within_err = r_within_ratio * (1./np.sqrt(r_within_r500) + 1./np.sqrt(r_cluster_no_bcgs))
	print 'RGZ r500: %i/%i (%.1f pm %.1f)%%' % (r_within_r500, r_cluster_no_bcgs, 100*r_within_ratio, 100*r_within_err)

	sep, n, bins = r500_hist(xmatch, {}, 20)
	bin_area = np.pi * np.array( [np.square(bins[0])] + [np.square(bins[i+1])-np.square(bins[i]) for i in np.arange(len(bins)-1)] )
	density = n/bin_area
	err = np.sqrt(n)/bin_area
	ax2.errorbar(np.log10(bins), np.log10(density), yerr=np.vstack((np.log10(density+err) - np.log10(density), np.log10(density) - np.log10(density-err))), fmt='o', ms=4, label='SDSS')
	logx = np.log10(bins[1:-3])
	logy = np.log10(density[1:-3])
	logyerr = err[1:-3] / density[1:-3]
	out = leastsq(errfunc, [-1,1], args=(logx, logy, logyerr), full_output=1)
	index, amp = out[0]
	index_err = np.sqrt(out[1][1][1])
	amp_err = np.sqrt(out[1][0][0])*amp
	ax2.plot(np.log10(bins), index*np.log10(bins)+amp, c='k', label='$\\alpha = $\n$  %.2f\pm%.2f$'%(index,index_err))
	err_max = np.maximum( (index+index_err)*np.log10(bins)+(amp+amp_err), (index-index_err)*np.log10(bins)+(amp+amp_err) )
	err_min = np.minimum( (index+index_err)*np.log10(bins)+(amp-amp_err), (index-index_err)*np.log10(bins)+(amp-amp_err) )
	#ax2.fill_between(np.log10(bins), err_min, err_max, color='k', lw=0, alpha=.5, label='$\\alpha = $\n$ %.2f\pm%.2f$'%(index,index_err))
	ax2.legend(loc='upper right')
	
	ax0.errorbar(np.log10(bins), np.log10(density), yerr=np.vstack((np.log10(density+err) - np.log10(density), np.log10(density) - np.log10(density-err))), fmt='o', ms=4, label='SDSS')
	ax0.plot(np.log10(bins), index*np.log10(bins)+amp, c='k', ls='--', label='$\\alpha = $\n$  %.2f\pm%.2f$'%(index,index_err))
	ax0.legend(loc='upper right')
	
	s_bcgs = n[0]
	s_cluster = sum(n[bins<10])
	s_bcg_ratio = 1.*s_bcgs/s_cluster
	s_bcg_err = s_bcg_ratio * (1./np.sqrt(s_bcgs) + 1./np.sqrt(s_cluster))
	print 'SDSS BCGs: %i/%i (%.1f pm %.1f)%%' % (s_bcgs, s_cluster, 100*s_bcg_ratio, 100*s_bcg_err)
	s_within_r500 = sum(n[bins<1][1:])
	s_cluster_no_bcgs = sum(n[bins<10][1:])
	s_within_ratio = 1.*s_within_r500/s_cluster_no_bcgs
	s_within_err = s_within_ratio * (1./np.sqrt(s_within_r500) + 1./np.sqrt(s_cluster_no_bcgs))
	print 'SDSS r500: %i/%i (%.1f pm %.1f)%%' % (s_within_r500, s_cluster_no_bcgs, 100*s_within_ratio, 100*s_within_err)
	
	ax.set_xlabel(get_label('WHL.r/r500', True))
	ax1.set_ylabel('$\\log_{10}$ (Surface density of galaxies [$r_{500}^{-2}$])')
	#ax.set_title('Surface density vs. separation')
	f1.tight_layout()
	
	ax0.set_xlabel(get_label('WHL.r/r500', True))
	ax0.set_ylabel('$\\log_{10}$ (Surface density of galaxies [$r_{500}^{-2}$])')
	#ax0.set_title('Surface density vs. separation')
	f2.tight_layout()
	
	r_s_bcgs = r_bcg_ratio / s_bcg_ratio
	r_s_bcgs_err = r_s_bcgs * (r_bcg_err/r_bcg_ratio + s_bcg_err/s_bcg_ratio)
	r_s_within = r_within_ratio / s_within_ratio
	r_s_within_err = r_s_within * (r_within_err/r_within_ratio + s_within_err/s_within_ratio)
	print 'BCG excess: %.3f pm %.3f\nr500 excess: %.3f pm %.3f' % (r_s_bcgs, r_s_bcgs_err, r_s_within, r_s_within_err)

def fractional_bent():
	
	params = total_cuts.copy()
	params['best.redshift']['$lte'] = 0.6

	sep, bend = [], []	
	for source in bending_15.find(total_cuts.copy()):
		sep.append(source['WHL']['r/r500'])
		bend.append(source['using_peaks']['bending_excess'])

	sep_field, bend_field = [], []
	for source in bent_sources.find(total_cuts.copy()):
		if bending_15.find_one({'RGZ.RGZ_id':source['RGZ']['RGZ_id']}) is None:
			sep_field.append(99)
			bend_field.append(source['using_peaks']['bending_excess'])

	sep = np.array(sep+sep_field)
	bend = np.array(bend+bend_field)

	straight = bend<0
	marginal = np.logical_and(0<bend, bend<5)
	bent = bend>5
	inner = sep<1.5#np.logical_and(sep>0.01, sep<1.5)
	cluster = sep<10
	
	num, denom = sum(sep<1.5), len(sep)
	print 'All: %.2f +- %.2f' % (1.*num/denom, 1.*num/denom * (1./np.sqrt(num) + 1./np.sqrt(denom)))
	num, denom = sum(np.logical_and(sep<1.5, bend>get_bent_cut())), sum(bend>get_bent_cut())
	print 'Highly bent: %.2f +- %.2f' % (1.*num/denom, 1.*num/denom * (1./np.sqrt(num) + 1./np.sqrt(denom)))
	num, denom = sum(np.logical_and(sep<1.5, bend>50)), sum(bend>50)
	print 'Max: %.2f +- %.2f' % (1.*num/denom, 1.*num/denom * (1./np.sqrt(num) + 1./np.sqrt(denom)))

	fig, ax = plt.subplots(1)
	n_all, bins, _ = ax.hist(np.log10(bend[bent]), bins=15, label='All sources')
	n_all = np.insert(n_all, 0, [sum(straight), sum(marginal)])
	n_inner, _, _ = ax.hist(np.log10(bend[np.logical_and(bent,inner)]), bins=bins, label='Sources within 1.5 r500')
	n_inner = np.insert(n_inner, 0, [sum(np.logical_and(straight,inner)), sum(np.logical_and(marginal,inner))])
	ax.set_xlabel('log(Bending excess [deg])')
	ax.set_ylabel('Count')
	#ax.set_title('Distributions of bending excess')
	ax.legend()
	plt.tight_layout()

	fig, ax = plt.subplots(1)
	lin_bins = np.insert(pow(10,(bins[:-1]+bins[1:])/2.), 0, [0, np.median(bend[marginal])])
	n_inner = np.array([ sum(bend[inner]>i) for i in lin_bins ], dtype=float)
	n_all = np.array([ sum(bend>i) for i in lin_bins ], dtype=float)
	ratio = n_inner/n_all
	err = ratio * np.sqrt(1/n_inner + 1/n_all)
	plt.errorbar(lin_bins, ratio, yerr=err)
	ax.set_xlabel('Excess bending angle (deg)')
	ax.set_ylabel('Fraction within $1.5~r_{500}$')
	#ax.set_title('Fraction of sources in BCG/inner regions of clusters')
	plt.tight_layout()

def histedges_equalN(x, nbin):
	npt = len(x)
	return np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x))

def import_env_csv1():
	params = total_cuts.copy()
	count = 0
	env = {}
	dirname = '/home/garon/Documents/RGZdata/bending/'
	for filename in ['distantstraight_env.csv', 'distantbent_env.csv']:
		with open(dirname+filename, 'r') as f:
			r = csv.reader(f)
			r.next()
			for objID, ra, dec, z, n_objID, n_ra, n_dec, n_z, z_type, cModelMag_r, fracDeV_r, dist in r:
				ra, dec, z = float(ra), float(dec), float(z)
				params['best.ra'] = {'$gte':ra-1e-3, '$lte':ra+1e-3}
				params['best.dec'] = {'$gte':dec-1e-3, '$lte':dec+1e-3}
				params['best.redshift'] = {'$gte':z-1e-3, '$lte':z+1e-3}
				source = bending_15.find_one(params)
				if source is not None:
					objID = source['SDSS']['objID']
					if objID not in env:
						env[objID] = {'objID':objID, 'ra':ra, 'dec':dec, 'z':z, 'bending_excess':source['using_peaks']['bending_excess'], 'neighbors':{}, 'r/r500':source['WHL']['r/r500']}
					if np.abs(z-float(n_z))/(1+z) <= 0.04:
						ADD = float(cosmo.angular_diameter_distance(z)/u.Mpc)
						env[objID]['neighbors'][n_objID] = {'ra':float(n_ra), 'dec':float(n_dec), 'z':float(n_z), 'z_type':z_type, 'cModelMag_r':float(cModelMag_r), 'fracDeV_r':float(fracDeV_r), 'dist_deg':float(dist), 'dist_Mpc':float(dist)*np.pi/180*ADD}
					if len(env)>count:
						count += 1
						print count
	
	db.drop_collection('distant_sources')
	distant_sources.create_index('objID', unique=True)
	for key in env:
		distant_sources.insert(env[key])

def import_env_csv2():
	count = 0
	env = {}
	dirname = '/home/garon/Documents/RGZdata/bending/'
	for filename in ['more_source_envPhotoz.csv']:
		with open(dirname+filename, 'r') as f:
			r = csv.reader(f)
			r.next()
			for objID, ra, dec, z, bending, r500, n_objID, n_ra, n_dec, n_z, z_type, cModelMag_r, fracDeV_r, dist in r:
				ra, dec, z, bending, r500 = float(ra), float(dec), float(z), float(bending), float(r500)
				if objID not in env:
					env[objID] = {'objID':objID, 'ra':ra, 'dec':dec, 'z':z, 'bending_excess':bending, 'r/r500':r500, 'neighbors':{}}
				if objID != n_objID and np.abs(z-float(n_z))/(1+z) <= 0.04:
					ADD = float(cosmo.angular_diameter_distance(z)/u.Mpc)
					env[objID]['neighbors'][n_objID] = {'ra':float(n_ra), 'dec':float(n_dec), 'z':float(n_z), 'z_type':z_type, 'cModelMag_r':float(cModelMag_r), 'fracDeV_r':float(fracDeV_r), 'dist_deg':float(dist), 'dist_Mpc':float(dist)*np.pi/180*ADD}
				if len(env)>count:
					count += 1
					print count
	
	#db.drop_collection('distant_sources')
	#distant_sources.create_index('objID', unique=True)
	for key in env:
		distant_sources.insert(env[key])

def import_env():
	import_env_csv1()
	import_env_csv2()
	for source in distant_sources.find({'r/r500':{'$exists':False}}):
		s2 = bending_15.find_one({'SDSS.objID':source['objID']})
		distant_sources.update({'_id':source['_id']}, {'$set':{'r/r500':s2['WHL']['r/r500']}})

def find_more_sources():
	bent_cut = get_bent_cut()
	with open('more_sources.csv', 'w') as f:
		print >> f, 'objID,ra,dec,z,bending,r/r500'
		count = distant_sources.find({'bending_excess':{'$gt':bent_cut}}).count()
		for source in bending_15.find(total_cuts):
			if source['using_peaks']['bending_excess']>30 and source['WHL']['r/r500']>6 and not distant_sources.find({'objID':source['SDSS']['objID']}).count():
				print >> f, '%s,%s,%s,%s,%s,%s' % (source['SDSS']['objID'], source['best']['ra'], source['best']['dec'], source['best']['redshift'], source['using_peaks']['bending_excess'], source['WHL']['r/r500'])
				count += 1
				if count>=150:
					break
	print count, 'bent'
	
	with open('more_sources.csv', 'a') as f:
		count = distant_sources.find({'bending_excess':{'$lt':bent_cut}}).count()
		for source in bending_15.find(total_cuts):
			if source['using_peaks']['bending_excess']<10 and 6<source['WHL']['r/r500']<15 and not distant_sources.find({'objID':source['SDSS']['objID']}).count():
				print >> f, '%s,%s,%s,%s,%s,%s' % (source['SDSS']['objID'], source['best']['ra'], source['best']['dec'], source['best']['redshift'], source['using_peaks']['bending_excess'], source['WHL']['r/r500'])
				count += 1
				if count>=150:
					break
	print count, 'straight'

def distant_env(bin_count=8, inner_cut=0.05, outer_cut=0.35):
	bent_cut = get_bent_cut()
	dist_bent, dist_straight, sep = [], [], []
	z_bent, z_straight = [], []
	for source in distant_sources.find():
		if source['bending_excess']>bent_cut:
			for neighbor in source['neighbors']:
				dist_bent.append(source['neighbors'][neighbor]['dist_Mpc'])
			z_bent.append(source['z'])
		else:
			for neighbor in source['neighbors']:
				dist_straight.append(source['neighbors'][neighbor]['dist_Mpc'])
			z_straight.append(source['z'])
		sep.append(source['r/r500'])
	
	dist_bent = np.array(dist_bent)
	dist_bent = dist_bent[dist_bent>inner_cut]
	dist_straight = np.array(dist_straight)
	dist_straight = dist_straight[dist_straight>inner_cut]
	for dist in [dist_bent, dist_straight]:
		dist.sort()

	n_bent = distant_sources.find({'bending_excess':{'$gt':bent_cut}}).count()
	n_straight = distant_sources.find({'bending_excess':{'$lt':bent_cut}}).count()
	density_bent = stats.rankdata(dist_bent) / (np.pi*dist_bent**2) / n_bent
	density_straight = stats.rankdata(dist_straight) / (np.pi*dist_straight**2) / n_straight

	plt.figure()
	plt.plot(dist_bent, density_bent, label='Bent sources')
	plt.plot(dist_straight, density_straight, ls='--', label='Straight sources')
	plt.legend()
	plt.xlabel('Projected distance from radio galaxy (Mpc)')
	plt.ylabel('Companion density (Mpc$^{-2}$)')
	plt.tight_layout()
	
	inner, outer = 2., 3.
	mask = np.logical_and(inner<dist_bent, dist_bent<outer)
	bg_count_bent = sum(mask)
	bg_bent = bg_count_bent / np.pi / (outer**2-inner**2)# / n_bent
	mask = np.logical_and(2<dist_straight, dist_straight<3)
	bg_count_straight = sum(mask)
	bg_straight = bg_count_straight / np.pi / (outer**2-inner**2)# / n_straight

	plt.figure()
	per_bin_bent, bins, _ = plt.hist(dist_bent, bins=np.linspace(inner_cut, 2, bin_count+1))
	bin_bg_bent = bg_bent*np.pi*(bins[1:]**2-bins[:-1]**2)
	frac_bent = per_bin_bent / bin_bg_bent
	frac_err_bent = frac_bent * (1./np.sqrt(per_bin_bent) + 1./np.sqrt(bg_bent))
	per_bin_straight, _, _ = plt.hist(dist_straight, bins=bins)
	bin_bg_straight = bg_straight*np.pi*(bins[1:]**2-bins[:-1]**2)
	frac_straight = per_bin_straight / bin_bg_straight
	frac_err_straight = frac_straight * (1./np.sqrt(per_bin_straight) + 1./np.sqrt(bg_straight))
	
	plt.figure()
	plt.errorbar(bins[1:], frac_bent, frac_err_bent, label='High bending sources')
	plt.errorbar(bins[1:]-0.01, frac_straight, frac_err_straight, ls='--', label='Low bending sources')
	plt.legend(loc='upper right')
	plt.xlabel('Projected distance from radio galaxy (Mpc)')
	plt.ylabel('Normalized surface density')
	plt.tight_layout()
	
	# 2D statistics
	area = np.pi * (outer_cut**2-inner_cut**2)
	
	cum_bent = sum(dist_bent<outer_cut)
	cum_excess_bent = cum_bent - bg_bent*area
	err_bent = np.sqrt(cum_bent+bg_bent*area**2)
	frac_bent = (cum_excess_bent/area) / bg_bent
	frac_excess_bent = frac_bent - 1.
	frac_err_bent = frac_bent * (err_bent/cum_excess_bent + 1./np.sqrt(bg_count_bent))
	
	cum_straight = sum(dist_straight<outer_cut)
	cum_excess_straight = cum_straight - bg_straight*area
	err_straight = np.sqrt(cum_straight+bg_straight*area**2)
	frac_straight = (cum_excess_straight/area) / bg_straight
	frac_excess_straight = frac_straight - 1.
	frac_err_straight = frac_straight * (err_straight/cum_excess_straight + 1./np.sqrt(bg_count_straight))
	
	print 'Surface density between %i to %i kpc' % (int(inner_cut*1000), int(outer_cut*1000))
	print 'Bent: %.3f +- %.3f (%.2f sigma)' % (frac_bent, frac_err_bent, frac_bent / frac_err_bent)
	print 'Straight: %.3f +- %.3f (%.2f sigma)' % (frac_straight, frac_err_straight, frac_straight / frac_err_straight)
	print 'Difference: %.3f +- %.3f (%.2f sigma)' % (frac_bent-frac_straight, np.sqrt(frac_err_bent**2 + frac_err_straight**2), (frac_bent-frac_straight) / np.sqrt(frac_err_bent**2 + frac_err_straight**2))
	
	# 3D statistics
	vol_sph = 4./3. * np.pi * (outer_cut**3-inner_cut**3)
	vol_cyl = np.pi * outer * (outer_cut**2-inner_cut**2)
	vol_bg = np.pi * outer * (outer**2-inner**2)

	bg_bent_vol = bg_count_bent / vol_bg
	cum_bent = sum(dist_bent<outer_cut)
	cum_excess_bent = cum_bent - bg_bent_vol*vol_cyl
	err_bent = np.sqrt(cum_bent + bg_bent_vol*vol_cyl**2)
	frac_bent = (cum_excess_bent/vol_sph) / bg_bent_vol
	frac_excess_bent = frac_bent - 1.
	frac_err_bent = frac_bent * (err_bent/cum_excess_bent + 1./np.sqrt(bg_count_bent))

	bg_straight_vol = bg_count_straight / vol_bg
	cum_straight = sum(dist_straight<outer_cut)
	cum_excess_straight = cum_straight - bg_straight_vol*vol_cyl
	err_straight = np.sqrt(cum_straight + bg_straight_vol*vol_cyl**2)
	frac_straight = (cum_excess_straight/vol_sph) / bg_straight_vol
	frac_excess_straight = frac_straight - 1.
	frac_err_straight = frac_straight * (err_straight/cum_excess_straight + 1./np.sqrt(bg_count_straight))

	print '\nVolume density between %i to %i kpc' % (int(inner_cut*1000), int(outer_cut*1000))
	print 'Bent: %.3f +- %.3f (%.2f sigma)' % (frac_bent, frac_err_bent, frac_bent / frac_err_bent)
	print 'Straight: %.3f +- %.3f (%.2f sigma)' % (frac_straight, frac_err_straight, frac_straight / frac_err_straight)
	print 'Difference: %.3f +- %.3f (%.2f sigma)' % (frac_bent-frac_straight, np.sqrt(frac_err_bent**2 + frac_err_straight**2), (frac_bent-frac_straight) / np.sqrt(frac_err_bent**2 + frac_err_straight**2))
	
	med_sep = np.median(sep)
	rho_bg = 500. / med_sep**2
	print '%.1f, %.1f rho_crit around high and low bending, respectively' % (frac_bent*rho_bg, frac_straight*rho_bg)
	
	return None
	
	step = inner_cut
	x = np.arange(step, 1, step)+step
	cum_excess_bent, cum_excess_straight = [], []
	frac_excess_bent, frac_excess_straight = [], []
	err_bent, err_straight = [], []
	frac_err_bent, frac_err_straight = [], []
	for i in x:
		cum_bent = sum(dist_bent<i)
		cum_excess_bent.append(cum_bent - bg_bent*np.pi*i**2)
		err_bent.append(np.sqrt(cum_bent+bg_bent*(np.pi*i**2)**2))
		frac_bent = cum_bent/(bg_bent*np.pi*i**2)
		frac_excess_bent.append(frac_bent - 1.)
		frac_err_bent.append(frac_bent * (1./np.sqrt(cum_bent) + 1./np.sqrt(bg_bent)))
		cum_straight = sum(dist_straight<i)
		cum_excess_straight.append(cum_straight - bg_straight*np.pi*i**2)
		err_straight.append(np.sqrt(cum_straight+bg_straight*(np.pi*i**2)**2))
		frac_straight = cum_straight/(bg_straight*np.pi*i**2)
		frac_excess_straight.append(frac_straight - 1.)
		frac_err_straight.append(frac_straight * (1./np.sqrt(cum_straight) + 1./np.sqrt(bg_straight)))

	# normalize
	cum_excess_bent_norm = np.array(cum_excess_bent)/n_bent
	err_bent_norm = cum_excess_bent_norm * (np.array(err_bent)/np.array(cum_excess_bent) + 1./np.sqrt(n_bent))
	cum_excess_straight_norm = np.array(cum_excess_straight)/n_straight
	err_straight_norm = cum_excess_straight_norm * (np.array(err_straight)/np.array(cum_excess_straight) + 1./np.sqrt(n_straight))

	# cumulative excess plot
	plt.figure()
	plt.errorbar(x, cum_excess_bent_norm, yerr=err_bent_norm, label='Bent sources')
	plt.errorbar(x+0.15*step, cum_excess_straight_norm, yerr=err_straight_norm, ls='--', label='Straight sources')
	#plt.axhline(0, ls=':', c='k')
	plt.legend(loc='upper left')
	plt.xlabel('Radius around radio galaxy (Mpc)')
	plt.ylabel('Cumulative excess density (Mpc$^{-2}$)')
	plt.tight_layout()
	
	# fractional excess plot
	plt.figure()
	plt.errorbar(x, frac_excess_bent, yerr=frac_err_bent, label='Bent sources')
	plt.errorbar(x+0.15*step, frac_excess_straight, yerr=frac_err_straight, ls='--', label='Straight sources')
	#plt.axhline(0, ls=':', c='k')
	plt.legend(loc='upper right')
	plt.xlabel('Radius around radio galaxy (Mpc)')
	plt.ylabel('Fractional excess density')
	plt.tight_layout()
	
	return None
	
	def sig(i, j):
		frac_b = cum_bent / bg_bent
		frac_s = cum_straight / bg_straight
		dfrac_b = frac_b * (1/np.sqrt(cum_bent) + 1/np.sqrt(bg_bent)) / np.sqrt(1.*i/n_bent)
		dfrac_s = frac_s * (1/np.sqrt(cum_straight) + 1/np.sqrt(bg_straight)) / np.sqrt(1.*j/n_straight)
		return (frac_b - frac_s) / max(dfrac_b, dfrac_s)

	# statistical excess plot
	plt.figure()
	plt.plot(x, cum_excess_bent_norm/err_bent_norm, label='Bent sources')
	plt.plot(x+0.15*step, cum_excess_straight_norm/err_straight_norm, ls='--', label='Straight sources')
	plt.legend(loc='upper right')
	plt.xlabel('Radius around RG (Mpc)')
	plt.ylabel('Standard deviations above background')
	plt.tight_layout()

	# statistical difference plot
	plt.figure()
	plt.plot(x, (cum_excess_bent_norm-cum_excess_straight_norm)/np.max([err_bent_norm,err_straight_norm], 0))
	plt.xlabel('Radius around RG (Mpc)')
	plt.ylabel('Standard deviations difference')
	plt.tight_layout()

	# distribution of bending vs density
	bent_companions, straight_companions = [], []
	for source in distant_sources.find():
		count = 0
		for neighbor in source['neighbors']:
			if 0.05 < source['neighbors'][neighbor]['dist_Mpc'] < 0.25:
				count += 1
		if source['bending_excess']>bent_cut:
			bent_companions.append(count)
		else:
			straight_companions.append(count)

	plt.figure()
	plt.hist(bent_companions, bins=range(13), alpha=.8, normed=True, label='Bent sources')
	plt.hist(straight_companions, bins=range(13), alpha=.8, normed=True, label='Straight sources')
	plt.scatter(6, .1, c='w', alpha=0, label='$p=%.3f$'%stats.anderson_ksamp([bent_companions, straight_companions]).significance_level)
	plt.xlabel('Number of companions between 50 and 250 kpc')
	plt.ylabel('Normalized count')
	plt.legend()
	plt.tight_layout()
	
	# bending vs density
	bending, neighbors = [], []
	for source in distant_sources.find():
		bending.append(source['bending_excess'])
		count = 0
		for neighbor in source['neighbors']:
			if 0.05 < source['neighbors'][neighbor]['dist_Mpc'] < 0.25:
				count += 1
		neighbors.append(count)

	bending = np.array(bending)
	neighbors = np.array(neighbors)
	order = bending.argsort()
	bg = (bg_bent+bg_straight)/(n_bent+n_straight)*np.pi*(0.25**2-0.05**2)

	plt.figure()
	n, bins, patches = plt.hist(bending)
	plt.figure()
	cut = 1.5
	plt.hist(bending[neighbors/bg>cut], bins=bins, alpha=.8, normed=True, label='Density $>%.1f\\times$background\n$n=%i$'%(cut,sum(neighbors/bg>cut)))
	plt.hist(bending[neighbors/bg<cut], bins=bins, alpha=.8, normed=True, label='Density $<%.1f\\times$background\n$n=%i$'%(cut,sum(neighbors/bg<cut)))
	plt.xlabel('Excess bending angle (deg)')
	plt.ylabel('Normalized count')
	plt.legend()
	plt.tight_layout()

def cluster_spotting(compress=False):
	n_bent = np.array([22, 53, 81, 101, 112, 121, 113, 162, 146, 207, 208, 210, 244, 228, 270, 275, 317, 271, 330, 328, 301, 326, 354, 378, 394, 414, 419, 435, 444, 529])
	n_straight = np.array([35, 51, 61, 94, 100, 114, 117, 154, 136, 169, 167, 216, 188, 234, 241, 255, 267, 282, 304, 309, 357, 335, 335, 366, 368, 411, 433, 484, 463, 516]) * 90./96. # normalized
	base_step_size = 0.1 # Mpc
	base_step_count = 30
	if compress:
		step_size = 2.*base_step_size
		step_count = base_step_count/2.
		n_bent = n_bent[::2] + n_bent[1::2]
		n_straight = n_straight[::2] + n_straight[1::2]
	else:
		step_size = base_step_size
		step_count = base_step_count
	
	bg_bent = 254
	bg_straight = 241
	sep = step_size*(np.arange(step_count)+1)
	excess_bent = n_bent - bg_bent * np.pi * (sep**2-(sep-step_size)**2)
	excess_straight = n_straight - bg_straight * np.pi * (sep**2-(sep-step_size)**2)
	cum_excess_bent = np.cumsum(excess_bent/90.)
	cum_excess_straight = np.cumsum(excess_straight/90.)
	err_bent = np.sqrt(n_bent + bg_bent)
	err_straight = np.sqrt(n_straight + bg_straight)
	cum_err_bent = np.sqrt(np.cumsum(err_bent/90.))
	cum_err_straight = np.sqrt(np.cumsum(err_straight/90.))
	print 'Bent excess within 2 Mpc: %.2f\pm%.2f\nStraight excess within 2 Mpc: %.2f\pm%.2f' % (cum_excess_bent[sep==2], cum_err_bent[sep==2], cum_excess_straight[sep==2], cum_err_straight[sep==2])
	
	diff = cum_excess_bent[sep==2] - cum_excess_straight[sep==2]
	diff_err = np.sqrt(cum_err_bent[sep==2] + cum_err_straight[sep==2])
	sig_level = cum_excess_bent[sep==2] - cum_excess_straight[sep==2]
	print 'Difference within 2 Mpc: %.2f\pm%.2f' % (diff, diff_err)
	
	fig, ax = plt.subplots(1)
	ax.errorbar(sep, cum_excess_bent, yerr=cum_err_bent, label='Bent sources')
	ax.errorbar(sep+0.2*base_step_size, cum_excess_straight, yerr=cum_err_straight, ls='--', label='Straight sources')
	ax.legend(loc='upper left')
	ax.set_xlabel('Distance from radio source (Mpc)')
	ax.set_ylabel('Cumulative excess per radio source')
	ax.set_title('Companions per Mpc$^2$ around non-cluster radio sources')

def bending_limit():
	h = 1. # jet width in kpc
	p_min = np.mean([0.9, 0.6, 1.4, 1.7, 0.4, 0.6, 1.4])*0.0062415 # minimum synchrotron pressure in keV cm^-3
	for source in bending_15.find():
		size = source['RGZ']['size_kpc']
		p_ram = source['WHL']['P']
		theta = np.arcsin(0.5 * size/h * p_ram/p_min)
		bending_15.update({'_id':source['_id']}, {'$set':{'WHL.bending_max':theta*180./np.pi if not np.isnan(theta) else 90.}})

def get_bent_cut():
	bend = []
	for i in bending_15.find(total_cuts):
		bend.append(i['using_peaks']['bending_excess'])

	bend = np.array(bend)
	sigma = 0.682689492137
	return np.percentile(bend, 100*(1+sigma)/2)

def remove_AllWISE(drop=False):
	count = bending_15.count({'SDSS':{'$exists':True}})
	if drop:
		print bending_15.update({}, {'$unset':{'best':None}}, multi=True)
	for source in bending_15.find({'best':{'$exists':False}}).batch_size(500):
		z, z_err = get_z(source)
		if z:
			sql = 'select raerr, decerr from photoprimary where objid=%i' % source['SDSS']['objID']
			df = SDSS_select(sql)
			ra_err = df['raerr'][0]
			dec_err = df['decerr'][0]
			best = {'redshift':z, 'redshift_err':z_err, 'ra':source['SDSS']['ra'], 'ra_err':ra_err, 'dec':source['SDSS']['dec'], 'dec_err':dec_err}
			bending_15.update({'_id':source['_id']}, {'$set':{'best':best}})
			print '%i/%i' % (bending_15.find({'best':{'$exists':True}}).count(), count)

def sep_hist():
	params = total_cuts.copy()
	params['using_peaks.bending_angle'] = params['using_peaks.bending_corrected']
	del params['using_peaks.bending_corrected']

	sep = []
	for source in bending_15.find(total_cuts):
		sep.append(source['WHL']['r/r500'])

	sep = np.array(sep)
	plt.figure()
	n, bins, _ = plt.hist(np.log10(sep), bins=20, fill=False, hatch='//')

	area = np.pi * (pow(10, bins)[1:]**2 - pow(10, bins)[:-1]**2)
	outer = np.percentile(sep, 95)
	density = sum(sep<outer)/(np.pi*outer**2)
	x = (pow(10, bins)[1:]+pow(10, bins)[:-1])/2.
	y = density*area
	popt = np.polyfit(x, y, 2)
	newx = np.logspace(np.log10(min(sep)), np.log10(outer), 100)
	plt.plot(np.log10(newx), popt[0]*newx**2 + popt[1]*newx + popt[2], label='Assuming uniform \ndistribution on the sky', c='b', ls=':')
	plt.ylim(ymax=1300)
	
	#plt.figure()
	#plt.hist(sep, bins=pow(10,bins), fill=False, hatch='//')
	#plt.xscale('log')
	plt.xlabel(get_label('WHL.r/r500', True))
	plt.ylabel('Count')
	plt.legend(loc='upper left')
	plt.tight_layout()

def density_ratio(x1, x2):
	c500, gamma, alpha, beta = 1.177, 0.3081, 1.0510, 5.4905
	rho1 = (c500*x1)**gamma * (1 + (c500*x1)**alpha)**(1.*(beta-gamma)/alpha)
	rho2 = (c500*x2)**gamma * (1 + (c500*x2)**alpha)**(1.*(beta-gamma)/alpha)
	return rho2 / rho1

def v500():
	# orbital velocity at r500
	c = np.sqrt(1.327e25) # km^1.5 s^-1
	m500, r500 = [], []
	for source in bending_15.find(total_cuts):
		m500.append(source['WHL']['M500'])
		r500.append(source['WHL']['r500'])
	
	m500 = np.array(m500)
	r500 = np.array(r500)
	v = c * np.sqrt(m500/r500/3.086e19) # km s^-1
	return np.median(v)

def scatter():
	bending, mass, sep, pop = [], [], [], []
	for source in bending_15.find(total_cuts):
		bending.append(source['using_peaks']['bending_corrected'])
		mass.append(source['WHL']['M500'])
		sep.append(source['WHL']['r/r500'])
		pop.append(source['WHL']['population'])
	
	bending = np.array(bending)
	mass = np.array(mass)
	sep = np.array(sep)
	pop = np.array(pop)
	
	plt.scatter(mass[pop=='outer'], bending[pop=='outer'], s=1, label='outer region')
	plt.scatter(mass[pop=='BCG'], bending[pop=='BCG'], s=1, label='BCGs')
	plt.scatter(mass[pop=='inner'], bending[pop=='inner'], s=1, label='inner region')
	plt.axhline(np.median(bending), ls=':', c='k')
	plt.legend()
	plt.xlabel(get_label('WHL.M500'))
	plt.ylabel(get_label('using_peaks.bending_corrected'))
	plt.tight_layout()

def get_phot_matches():
	z_r, z_w, cid = [], [], []
	for source in bending_15.find(total_cuts):
		z_r.append(source['best']['redshift'])
		z_w.append(source['WHL']['zbest'])
		cid.append(source['RGZ']['RGZ_id'])
	z_r = np.array(z_r)
	z_w = np.array(z_w)
	cid = np.array(cid)
	bad_cid = list(cid[np.abs(z_w - z_r) > .04*(1+z_r)])
	return bad_cid

def flag_dups(coll=bending_15):
	coll.update({}, {'$set':{'RGZ.duplicate':0}}, multi=True)
	cid, name = [], []
	for source in coll.find(total_cuts):
		name.append(source['RGZ']['RGZ_name'])
	for n in set(name):
		name.remove(n)
	for source in coll.find({'RGZ.RGZ_name':{'$in':name}}):
		cid.append(source['RGZ']['RGZ_id'])
	for source in coll.find({'RGZ.RGZ_id':{'$in':cid}}):
		if source['RGZ']['RGZ_name'] in name:
			coll.update({'RGZ.RGZ_id':source['RGZ']['RGZ_id']}, {'$set':{'RGZ.duplicate':1}})
			name.remove(source['RGZ']['RGZ_name'])

def print_supplement(filename='/home/garon/Documents/RGZdata/bending/data_supplement.csv'):
	bad_cid = get_phot_matches()
	cids = []
	with open(filename[:-4] + '_table1' + filename[-4:], 'w') as f1:
		with open(filename[:-4] + '_sample1' + filename[-4:], 'w') as g1:
			for source in bending_15.find(total_cuts).sort('RGZ.RGZ_name', 1):
				cids.append(source['RGZ']['RGZ_id'])
				morph = 2 if source['RGZ']['morphology'] == 'double' else 3
				ztype = 's' if 'spec_redshift' in source['SDSS'] else 'p'
				align = 'r' if source['WHL']['alignment'] == 'radial' else ('t' if source['WHL']['alignment'] == 'tangential' else 'o')
				if source['RGZ']['RGZ_id'] in bad_cid:
					w_z = source['WHL']['zphot']
					w_ztype = 'p'
				else:
					w_z = source['WHL']['zbest']
					w_ztype = 's' if 'zspec' in source['WHL'] else 'p'
				datastr = '%s,%s,%i,%.3f,%.5f,%.5f,%.4f,%.4f,%s,%.1f,%.1f,%.1f,%.2f,%s,%.5f,%.5f,%.4f,%s,%.2f,%.2f,%.2f,%.2f,%.1f,%s' % (source['RGZ']['RGZ_name'][3:], source['RGZ']['zooniverse_id'], morph, source['RGZ']['size_arcmin'], source['best']['ra'], source['best']['dec'], source['best']['redshift'], source['best']['redshift_err'], ztype, source['using_peaks']['bending_angle'], source['using_peaks']['bending_corrected'], source['using_peaks']['bending_excess'], source['using_contour']['asymmetry'], source['WHL']['name'], source['WHL']['ra'], source['WHL']['dec'], w_z, w_ztype, source['WHL']['r/r500'], source['WHL']['r500'], source['WHL']['M500'], np.log10(source['WHL']['P']), source['WHL']['orientation_peaks'], align)
				print >> f1, datastr
				print >> g1, to_tex(datastr)
	
	with open(filename[:-4] + '_table2' + filename[-4:], 'w') as f2:
		with open(filename[:-4] + '_sample2' + filename[-4:], 'w') as g2:
			params = total_cuts.copy()
			params['RGZ.RGZ_id'] = {'$nin':cids}
			for source in bent_sources.find(params).sort('RGZ.RGZ_name', 1):
				morph = 2 if source['RGZ']['morphology'] == 'double' else 3
				ztype = 's' if 'spec_redshift' in source['SDSS'] else 'p'
				datastr = '%s,%s,%i,%.3f,%.5f,%.5f,%.4f,%.4f,%s,%.1f,%.1f,%.1f,%.2f' % (source['RGZ']['RGZ_name'][3:], source['RGZ']['zooniverse_id'], morph, source['RGZ']['size_arcmin'], source['best']['ra'], source['best']['dec'], source['best']['redshift'], source['best']['redshift_err'], ztype, source['using_peaks']['bending_angle'], source['using_peaks']['bending_corrected'], source['using_peaks']['bending_excess'], source['using_contour']['asymmetry'])
				print >> f2, datastr
				print >> g2, to_tex(datastr)

def to_tex(s):
	keylist = [(',',' & '), ('& s','& $s$'), ('& p','& $p$'), ('& r','& $r$'), ('& t','& $t$'), ('& o','& $o$')]
	for key in keylist:
		s = s.replace(*key)
	return s + ' \\\\'

def make_all_figs():
	contamination()
	sep_hist()
	sdss_density()
	rmsd()
	bending_correct(plot=True, methods='using_peaks')
	plot_running('WHL.r/r500', 'using_peaks.bending_excess', logx=True, pop='separate', title=False)
	plot_running('WHL.M500', 'using_peaks.bending_excess', pop='inner', title=False)
	plot_running('WHL.M500', 'using_peaks.bending_excess', pop='BCG', title=False)
	plot_running('WHL.P', 'using_peaks.bending_excess', logx=True, pop='non-BCG', combined=True, title=False)
	orientation()
	orientation_test()
	fractional_bent()
	distant_env()

def update_SDSS():
	for entry in catalog.find({'SDSS':{'$exists':False}, 'consensus.ir_ra':{'$exists':True}}).batch_size(100):
		print entry['catalog_id']
		sdss_match = getSDSS(entry)
		if sdss_match is not None:
			z = sdss_match['spec_redshift'] if 'spec_redshift' in sdss_match else (sdss_match['photo_redshift'] if 'photo_redshift' in sdss_match else 0)
			if z > 0:
				radio = entry['radio']
				radio_nested = {'radio':radio}
				physical = getPhysical(z, radio_nested)
				radio.update(physical)
				catalog.update({'_id':entry['_id']}, {'$set':{'radio':radio, 'SDSS':sdss_match}})
			else:
				catalog.update({'_id':entry['_id']}, {'$set':{'SDSS':sdss_match}})

def vienna(data_in=False, data_out=False):
	infile = '/home/garon/Downloads/vienna.csv'
	outfile = '/home/garon/Downloads/vienna_matches.csv'
	
	if data_in:
		id, ra, dec, z = [], [], [], []
		with open(infile, 'r') as f:
			r = csv.reader(f)
			r.next()
			for row in r:
				id.append(row[0])
				ra.append(row[1])
				dec.append(row[2])
				z.append(row[3])
		id = np.array(id, dtype=int)
		ra = np.array(ra, dtype=float)
		dec = np.array(dec, dtype=float)
		z = np.array(z, dtype=float)
	
	if data_out:
		with open(outfile, 'w') as f:
			print >> f, 'source#,radeg,decdeg,z,sep_mpc,sep_r500'
			for i in range(len(id)):
				loc = coord.SkyCoord(ra[i], dec[i], unit=(u.deg,u.deg))
				w = get_whl(loc, z[i], 0, 15, 0.04*(1+z[i]))
				if w is not None:
					loc2 = coord.SkyCoord(w['RAdeg'], w['DEdeg'], unit=(u.deg,u.deg))
					sep_deg = loc2.separation(loc)
					sep_mpc = float(cosmo.angular_diameter_distance(w['zspec'] if 'zspec' in w else w['zphot'])/u.Mpc * sep_deg.to(u.rad)/u.rad)
					sep_r500 = sep_mpc / w['r500']
				else:
					sep_mpc, sep_r500 = 99., 99.
				print >> f, '%i,%f,%f,%f,%f,%f' % (id[i], ra[i], dec[i], z[i], sep_mpc, sep_r500)

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
