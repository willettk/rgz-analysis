import logging, urllib2, time, json, os
import pymongo
import numpy as np
import StringIO, gzip
from astropy.io import fits
from astropy import wcs, coordinates as coord, units as u
from astropy.cosmology import Planck13 as cosmo
from scipy.optimize import brentq
from scipy.interpolate import interp1d

from consensus import rgz_path, data_path, db, version
import contour_path_object as cpo
completed_file = '%s/bending_completed.txt' % rgz_path

# Connect to Mongo database
subjects = db['radio_subjects']
consensus = db['consensus{}'.format(version)]
catalog = db['catalog{}'.format(version)]
whl = db['WHL15']
rm_m = db['redmapper_members']
rm_c = db['redmapper_clusters']
amf = db['AMFDR9']
bent_sources = db['bent_sources']
bending_15 = db['bending_15']
bending_control = db['bending_control']

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
	if opening_angle > np.pi*u.rad:
		opening_angle = coord.Angle(2*np.pi*u.rad - opening_angle)
	bisector = (pos_angle_1+pos_angle_0)/2.
	if np.abs(bisector-pos_angle_0) > np.pi/2*u.rad:
		bisector += np.pi*u.rad
	bending_angles = {'pos_angle_0':pos_angle_0, 'pos_angle_1':pos_angle_1, 'opening_angle':opening_angle, 'bisector':bisector.wrap_at(2*np.pi*u.rad)}
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
	
	# If they cross exactly once, find the global solution
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
	ir = coord.SkyCoord(source['SDSS']['ra'], source['SDSS']['dec'], unit=(u.deg,u.deg), frame='icrs') if 'SDSS' in source else \
		 coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
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
	bending_angles = get_angles(w, ir, 'contour', contour_list)
	tail_lengths_apparent = get_tail_lengths(w, ir, 'contour', contour_list)
	tail_lengths_physical = []
	for tail in tail_lengths_apparent:
			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
	ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
	asymmetry = ratios[1]/ratios[0]
	using_contour = {'tail_deg_0':tail_lengths_apparent[0], 'tail_deg_1':tail_lengths_apparent[1], 'size_deg':sum(tail_lengths_apparent), 'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
	using_contour.update(bending_angles)
	for key in using_contour.keys():
		if type(using_contour[key]) is coord.angles.Angle:
			using_contour[key] = using_contour[key].deg
	
	# Using the 'peak' method
	bending_angles = get_angles(w, ir, 'peak', peaks)
	tail_lengths_apparent = get_tail_lengths(w, ir, 'peak', contour_list, peaks)
	tail_lengths_physical = []
	for tail in tail_lengths_apparent:
			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
	ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
	asymmetry = ratios[1]/ratios[0]
	using_peaks = {'tail_deg_0':tail_lengths_apparent[0], 'tail_deg_1':tail_lengths_apparent[1], 'size_deg':sum(tail_lengths_apparent), 'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
	using_peaks.update(bending_angles)
	for key in using_peaks.keys():
		if type(using_peaks[key]) is coord.angles.Angle:
			using_peaks[key] = using_peaks[key].deg
	
	morphology = 'double' if peak_count == 2 else 'triple'
	rgz = {'RGZ_id':source['catalog_id'], 'zooniverse_id':source['zooniverse_id'], 'ir_consensus':source['consensus']['ir_level'], 'radio_consensus':source['consensus']['radio_level'], 'peaks':peaks, 'components':source['radio']['components'], 'morphology':morphology}
	entry = {'RGZ':rgz, 'using_contour':using_contour, 'using_peaks':using_peaks}
	if 'SDSS' in source:
		entry['SDSS'] = source['SDSS']
	if 'AllWISE' in source:
		entry['AllWISE'] = source['AllWISE']
	
	best_ra = entry['SDSS']['ra'] if 'SDSS' in entry else entry['AllWISE']['ra']
	best_dec = entry['SDSS']['dec'] if 'SDSS' in entry else entry['AllWISE']['dec']
	best = {'ra':best_ra, 'dec':best_dec, 'redshift':z}
	entry.update({'best':best})
	
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
	double_args = {'$and': [{'ignore_bending':False, 'overedge':0, 'catalog_id':{'$nin':completed}}, \
							{'$or':  [{'radio.number_peaks':2, 'radio.number_components':1}, \
						 			  {'radio.number_components':2}]}, \
							{'$or':  [{'SDSS.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
									  {'SDSS.spec_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
									  {'AllWISE.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}] }]}
	triple_args = {'$and': [{'ignore_bending':False, 'overedge':0, 'catalog_id':{'$nin':completed}}, \
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
				print >> f, source['catalog_id']
				if entry is not None:
					count += 1
					output('%i %s' % (count, source['zooniverse_id']))
					bent_sources.insert(entry)

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
		whl_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg, 'r/r500':c_sep_mpc/whl_prop['r500']}
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
		tail_lengths_apparent = [source['using_contour']['tail_deg_0']*u.deg, source['using_contour']['tail_deg_1']*u.deg]
		tail_lengths_physical = []
		for tail in tail_lengths_apparent:
			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
		using_contour = {'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc)}
		source['using_contour'].update(using_contour)
		
		# Using the 'peak' method
		tail_lengths_apparent = [source['using_peaks']['tail_deg_0']*u.deg, source['using_peaks']['tail_deg_1']*u.deg]
		tail_lengths_physical = []
		for tail in tail_lengths_apparent:
			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.to(u.rad)/u.rad)
		using_peaks = {'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'size_kpc':float(sum(tail_lengths_physical)/u.kpc)}
		source['using_peaks'].update(using_peaks)
		
		# Calculate the orientation angle
#		for cluster,prop in zip([cluster_w,cluster_r,cluster_a], [whl_prop,rm_prop,amf_prop]):
#			if cluster is not None:
		cluster, prop = cluster_w, whl_prop
		for method,using in zip(['contour','peaks'], [source['using_contour'],source['using_peaks']]):
			orientation = coord.angles.Angle(prop['position_angle'] - using['bisector'], unit=u.deg).wrap_at(360.*u.deg).deg
			if orientation > 180.:
				orientation = 360. - orientation
			prop['orientation_%s' % method] = orientation
		
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
				print >> f, source['RGZ']['RGZ_id']
				if entry is not None:
					count += 1
					output('%i %s' % (count, source['RGZ']['zooniverse_id']))
					bending_15.insert(entry)

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

def to_file(filename, collection):
	'''
	Print the bending collection to a csv file for analysis
	'''
	rgz_keys = ['RGZ_id', 'zooniverse_id', 'morphology', 'radio_consensus', 'ir_consensus']
	best_keys = ['ra', 'dec', 'redshift']
	sdss_keys = ['ra', 'dec', 'objID', 'photo_redshift', 'photo_redshift_err', 'spec_redshift', 'spec_redshift_err', 'u', 'u_err', 'g', 'g_err', 'r', 'r_err', 'i', 'i_err', 'z', 'z_err', 'number_matches']
	wise_keys = ['ra', 'dec', 'designation', 'photo_redshift', 'w1mpro', 'w1sigmpro', 'w1snr', 'w2mpro', 'w2sigmpro', 'w2snr', 'w3mpro', 'w3sigmpro', 'w3snr', 'w4mpro', 'w4sigmpro', 'w4snr']
	whl_keys = ['name', 'ra', 'dec', 'zphot', 'zspec', 'N500', 'N500sp', 'RL*500', 'r500', 'separation_deg', 'separation_Mpc', 'r/r500', 'position_angle', 'orientation_contour', 'orientation_peaks']
	rm_keys = ['name', 'ra', 'dec', 'zlambda', 'zspec', 'S', 'lambda', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	amf_keys = ['AMF_id', 'ra', 'dec', 'z', 'r200', 'richness', 'core_radius', 'concentration', 'likelihood', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	bending_keys = ['pos_angle_0', 'pos_angle_1', 'opening_angle', 'bisector', 'tail_deg_0', 'tail_deg_1', 'size_deg', 'tail_kpc_0', 'tail_kpc_1', 'size_kpc', 'ratio_1', 'ratio_0', 'asymmetry']
	
	all_keys = [rgz_keys, best_keys, sdss_keys, wise_keys, whl_keys, bending_keys, bending_keys]
	dict_names = ['RGZ', 'best', 'SDSS', 'AllWISE', 'WHL', 'using_contour', 'using_peaks']
	
	success = 0
	with open(filename, 'w') as f:
		
		header = ''
		for superkey, key_list in zip(dict_names, all_keys):
			for key in key_list:
				header += '%s.%s,' % (str(superkey), str(key))
		print >> f, header[:-1]
		
		for entry in collection.find():
			try:
				row = ''
				for superkey, key_list in zip(dict_names, all_keys):
					for key in key_list:
						if superkey in entry and key in entry[superkey]:
							row += '%s,' % str(entry[superkey][key])
						else:
							row += '-99,'
				print >> f, row[:-1]
				success += 1
			except BaseException as e:
				output(e, logging.exception)
		
		output('%i/%i successfully printed to %s' % (success, collection.count(), filename))

if __name__ == '__main__':
	
	logging.basicConfig(filename='%s/bending.log' % rgz_path, level=logging.DEBUG, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	logging.captureWarnings(True)
	logging.info('Bending run from command line')
	
	# Run options
	calculate_bending = False
	match_to_clusters = False
	control = True
	
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
