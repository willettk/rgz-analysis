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
completed_file = '%s/bending_control_15.txt' % rgz_path

# Connect to Mongo database
subjects = db['radio_subjects']
consensus = db['consensus{}'.format(version)]
catalog = db['catalog{}'.format(version)]
whl = db['WHL15']
rm_m = db['redmapper_members']
rm_c = db['redmapper_clusters']
amf = db['AMFDR9']
xmatch = db['cluster_xmatch']
bending_coll = db['bending_control']

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

def get_contours(w, ir, peaks, data, peak_count):
	'''
	Returns a list of Path objects corresponding to each outer contour in the data, in RA and dec coordinates
	Removes outer layers until there are two components and a disjoint IR
	Removes any contours that aren't in this source
	'''
	assert (peak_count in [2,3]), 'Not a valid morphology'
	ir_pos = w.wcs_world2pix(np.array([[ir.ra.deg,ir.dec.deg]]), 1)
	peak_pos = w.wcs_world2pix(np.array([ [peak['ra'],peak['dec']] for peak in peaks ]), 1)
	
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

def get_bending_angles(w, ir, method, method_data):
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

def get_global_peaks(w, peaks, contour_tree):
	'''
	Determines the position of the global maximum for each component in the contour
	'''
	peak_pos = w.wcs_world2pix(np.array([ [peak['ra'],peak['dec']] for peak in peaks ]), 1)
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
		intersections = np.zeros(len(limits)-2)
		for ix in range(len(intersections)):
			intersections[ix] = brentq(diff, limits[ix], limits[ix+1])
		return intersections

def get_colinear_separation(w, ir, peak, contour):
	'''
	Finds the distance from the host to the edge of the contour, passing through the peak
	'''
	peak_pix = np.array(w.wcs_world2pix(peak['ra'], peak['dec'], 1))
	ir_pix = np.array(w.wcs_world2pix(ir.ra.deg, ir.dec.deg, 1))
	
	# Extrapolate the line connecting the peak to the IR position
	slope = (peak_pix[1]-ir_pix[1])/(peak_pix[0]-ir_pix[0])
	extrap_pos = ir_pix + w._naxis1*np.array([1.,slope])
	extrap_neg = ir_pix - w._naxis1*np.array([1.,slope])
	
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
	
	if peak_pix[0] > ir_pix[0]:
		tail = np.vstack((extrap_neg, ir_pix, peak_pix, extrap_pos))
	else:
		tail = np.vstack((extrap_pos, ir_pix, peak_pix, extrap_neg))
	
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

def z_hist(r_max, pathname='.', z_range=[-np.inf,np.inf]):
	'''
	Plot the histogram of (cluster_z-galaxy_z) for galaxies matched within a radius r_max of their corresponding cluster
	'''
	z_whl, z_rm = [], []
	dz_whl, dz_rm = [], []
	for source in catalog.find({'ignore_bending':False, \
								'$or': [{'SDSS.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
										{'SDSS.spec_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
										{'AllWISE.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}] }):
		ir = coord.SkyCoord(source['SDSS']['ra'], source['SDSS']['dec'], unit=(u.deg,u.deg), frame='icrs') if 'SDSS' in source else \
			 coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
		z, z_err = get_z(source)
		
		# Using WHL15
		whl_match = get_whl(ir, z, z_err, r_max, 1)
		if whl_match is not None:
			cluster_z = whl_match['zspec'] if 'zspec' in whl_match else whl_match['zphot']
			z_whl.append(z)
			dz_whl.append(cluster_z-z)
		
		# Using redMaPPer
		rm_match = get_redmapper(source['SDSS']['objID'] if 'SDSS' in source else None, ir, z, z_err, r_max, 1)
		if rm_match is not None:
			cluster_z = rm_match['zspec'] if 'zspec' in rm_match else rm_match['zlambda']
			z_rm.append(z)
			dz_rm.append(cluster_z-z)
	
	with open('%s/dz_whl_%.1f-%.1f.csv' % (pathname,z_range[0],z_range[1]), 'w') as f:
		print >> f, 'z,dz'
		for z,dz in zip(z_whl,dz_whl):
			print >> f, z, dz
	
	with open('%s/dz_rm_%.1f-%.1f.csv' % (pathname,z_range[0],z_range[1]), 'w') as f:
		print >> f, 'z,dz'
		for z,dz in zip(z_rm,dz_rm):
			print >> f, z, dz
	
	# Read it back from the file
	if False:
		z_whl, z_rm = [], []
		dz_whl, dz_rm = [], []
		with open('%s/dz_whl_%.1f-%.1f.csv' % (pathname,z_range[0],z_range[1]), 'r') as f:
			l = f.readline()
			l = f.readline()
			while(l != ''):
				vals = l.split(' ')
				z_whl.append(float(vals[0].strip()))
				dz_whl.append(float(vals[1].strip()))
				l = f.readline()
		with open('%s/dz_rm_%.1f-%.1f.csv' % (pathname,z_range[0],z_range[1]), 'r') as f:
			l = f.readline()
			l = f.readline()
			while(l != ''):
				vals = l.split(' ')
				z_rm.append(float(vals[0].strip()))
				dz_rm.append(float(vals[1].strip()))
				l = f.readline()

def cluster_xmatch():
	'''
	Cross-match two cluster catalogs, or a new cluster catalog to the existing collection
	Definition: centers within a 0.5 Mpc square and z within 0.04(1+z) (WHL def.)
	'''
	
	# Cross-match the catalogs (currently redMaPPer to WHL15)
	for cluster1 in rm_c.find():
		z1 = cluster1['zspec'] if 'zspec' in cluster1 else cluster1['zlambda']
		ra1 = cluster1['RAdeg']
		dec1 = cluster1['DEdeg']
		
		# Check if WHL15 has overlap
		max_sep = float(0.5 * u.Mpc / cosmo.angular_diameter_distance(z1) * u.rad / u.deg)
		delta_z = 0.04*(1+z1)
		cluster2 = whl.find_one({'RAdeg':{'$gt':ra1-max_sep, '$lt':ra1+max_sep}, 'DEdeg':{'$gt':dec1-max_sep, '$lt':dec1+max_sep}, \
								 '$or':[{'zspec':{'$gt':z1-delta_z, '$lt':z1+delta_z}}, {'zphot':{'$gt':z1-delta_z, '$lt':z1+delta_z}}]})
		if cluster2 is not None:
			
			# Insert these clusters to the xmatch collection
			z2 = cluster2['zspec'] if 'zspec' in cluster2 else cluster2['zphot']
			ra2 = cluster2['RAdeg']
			dec2 = cluster2['DEdeg']
			entry = {'redMaPPer': {'ra':ra1, 'dec':dec1, 'z':z1, '_id':cluster1['_id']}, \
					 'WHL':		  {'ra':ra2, 'dec':dec2, 'z':z2, '_id':cluster2['_id']} }
			xmatch.insert(entry)

def richness():
	'''
	Get an analytical formula to convert between WHL and redMaPPer richness values
	'''
	rich_whl, rich_rm = [], []
	#for pair in xmatch.find():
	#	rich_whl.append(whl.find_one({'_id':pair['WHL']['_id']})['RL*500'])
	#	rich_rm.append(rm_c.find_one({'_id':pair['redMaPPer']['_id']})['lambda'])
	with open('/home/garon/Documents/RGZdata/bending/richness.csv', 'r') as f:
		l = f.readline()
		l = f.readline()
		while(l != ''):
			vals = l.split(',')
			rich_whl.append(float(vals[0].strip()))
			rich_rm.append(float(vals[1].strip()))
			l = f.readline()
	
	rich_whl = np.array(rich_whl)
	rich_rm = np.array(rich_rm)
	
	# Truncated regression model
	from sklearn.linear_model import SGDRegressor
	model = SGDRegressor(eta0=1e-6, n_iter=55)
	model.fit(rich_whl.reshape(-1,1), rich_rm)
	xx = np.unique(rich_whl)
	yy = model.predict(xx.reshape(-1,1))
	a,b = np.polyfit(xx,yy,1)

def get_entry(source, peak_count, method=None):
	'''
	Given a source from RGZ and a method of measuring the tail (either 'contour' or 'peak'),
	calculate the bending parameters and match it to a cluster
	'''
	assert (peak_count in [2, 3]), 'Not a valid morphology'
	assert (method in ['contour', 'peak', None]), 'Not a valid method'
	subject = subjects.find_one({'zooniverse_id':source['zooniverse_id']})
	
	# Get pixel-to-WCS conversion
	fid = subject['metadata']['source']
	fits_loc = pathdict[fid]
	w = wcs.WCS(fits.getheader(fits_loc, 0))
	
	# Get the location of the source
	ir = coord.SkyCoord(source['SDSS']['ra'], source['SDSS']['dec'], unit=(u.deg,u.deg), frame='icrs') if 'SDSS' in source else \
		 coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
	z, z_err = get_z(source)
	
	# Get image parameters for this source
	data = get_data(subject)
	contour_tree = get_contours(w, ir, source['radio']['peaks'], data, peak_count)
	peaks = get_global_peaks(w, source['radio']['peaks'], contour_tree)
	if len(peaks) != 2:
		output("%s didn't have 2 tails" % source['zooniverse_id'])
		return
	peak_pos = w.wcs_world2pix(np.array([ [peak['ra'],peak['dec']] for peak in peaks ]), 1)
	contour_list = [child.path for child in contour_tree.children if any(child.contains(peak_pos))]
	
	# Match to cluster catalogs
	cluster_w = get_whl(ir, z, z_err, 15, 10.04*(1+z))
	whl_prop = {}
	if cluster_w is not None:
		c_pos = coord.SkyCoord(cluster_w['RAdeg'], cluster_w['DEdeg'], unit=(u.deg,u.deg), frame='icrs')
		c_sep_arc = c_pos.separation(ir)
		c_sep_mpc = float(cosmo.angular_diameter_distance(cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot'])/u.Mpc * c_sep_arc.to(u.rad)/u.rad)
		c_pos_angle = c_pos.position_angle(ir)
		whl_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg}
		for key in ['_id', 'N500', 'N500sp', 'RL*500', 'name', 'r500', 'zphot', 'zspec']:
			if key in cluster_w:
				whl_prop[key] = cluster_w[key]
	
	objID = source['SDSS']['objID'] if 'SDSS' in source else None
	cluster_r = get_redmapper(objID, ir, z, z_err, 15, 10.04*(1+z))
	rm_prop = {}
	if cluster_r is not None:
		c_pos = coord.SkyCoord(cluster_r['RAdeg'], cluster_r['DEdeg'], unit=(u.deg,u.deg), frame='icrs')
		c_sep_arc = c_pos.separation(ir)
		c_sep_mpc = float(cosmo.angular_diameter_distance(cluster_r['zspec'] if 'zspec' in cluster_r else cluster_r['zlambda'])/u.Mpc * c_sep_arc.to(u.rad)/u.rad)
		c_pos_angle = c_pos.position_angle(ir)
		rm_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg, 'objID':cluster_r['ObjID'], 'name':cluster_r['Name']}
		for key in ['_id', 'zlambda', 'zspec', 'lambda', 'S']:
			if key in cluster_r:
				rm_prop[key] = cluster_r[key]
	
	cluster_a = get_amf(ir, z, z_err, 15, 10.04*(1+z))
	amf_prop = {}
	if cluster_a is not None:
		c_pos = coord.SkyCoord(cluster_a['ra'], cluster_a['dec'], unit=(u.deg,u.deg), frame='icrs')
		c_sep_arc = c_pos.separation(ir)
		c_sep_mpc = float(cosmo.angular_diameter_distance(cluster_a['z'])/u.Mpc * c_sep_arc.to(u.rad)/u.rad)
		c_pos_angle = c_pos.position_angle(ir)
		amf_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg}
		for key in ['_id', 'r200', 'richness', 'core_radius', 'AMF_id', 'concentration', 'likelihood']:
			if key in cluster_a:
				amf_prop[key] = cluster_a[key]
	
	# Only continue if a cluster was matched
	if cluster_w is None and cluster_r is None and cluster_a is None:
		
		output("%s didn't match to a cluster" % source['zooniverse_id'])
		return
	
	else:
	
		# Using the 'contour' method
		#if method == 'contour':
		bending_angles = get_bending_angles(w, ir, 'contour', contour_list)
		tail_lengths_apparent = get_tail_lengths(w, ir, 'contour', contour_list)
		tail_lengths_physical = []
		for tail in tail_lengths_apparent:
			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.rad)
		ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
		asymmetry = ratios[1]/ratios[0]
		using_contour = {'tail_deg_0':tail_lengths_apparent[0], 'tail_deg_1':tail_lengths_apparent[1], 'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
		using_contour.update(bending_angles)
		for key in using_contour.keys():
			if type(using_contour[key]) is coord.angles.Angle:
				using_contour[key] = using_contour[key].deg
		
		# Using the 'peak' method
		#elif method == 'peak':
		bending_angles = get_bending_angles(w, ir, 'peak', peaks)
		tail_lengths_apparent = get_tail_lengths(w, ir, 'peak', contour_list, peaks)
		tail_lengths_physical = []
		for tail in tail_lengths_apparent:
			tail_lengths_physical.append(cosmo.angular_diameter_distance(z) * tail.rad)
		ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
		asymmetry = ratios[1]/ratios[0]
		using_peaks = {'tail_deg_0':tail_lengths_apparent[0], 'tail_deg_1':tail_lengths_apparent[1], 'tail_kpc_0':float(tail_lengths_physical[0]/u.kpc), 'tail_kpc_1':float(tail_lengths_physical[1]/u.kpc), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
		using_peaks.update(bending_angles)
		using_peaks.update(bending_angles)
		for key in using_peaks.keys():
			if type(using_peaks[key]) is coord.angles.Angle:
				using_peaks[key] = using_peaks[key].deg
		
		# Tail length properties
		#ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
		#asymmetry = ratios[1]/ratios[0]
		
		# Calculate the orientation angle
		for cluster,prop in zip([cluster_w,cluster_r,cluster_a], [whl_prop,rm_prop,amf_prop]):
			if cluster is not None:
				for method,using in zip(['contour','peaks'], [using_contour,using_peaks]):
					orientation = coord.angles.Angle(prop['position_angle'] - using['bisector'], unit=u.deg).wrap_at(360.*u.deg).deg
					if orientation > 180.:
						orientation = 360. - orientation
					prop['orientation_%s' % method] = orientation
		
		# Compile all results
		morphology = 'double' if peak_count == 2 else 'triple'
		rgz = {'RGZ_id':source['catalog_id'], 'zooniverse_id':source['zooniverse_id'], 'ir_consensus':source['consensus']['ir_level'], 'radio_consensus':source['consensus']['radio_level'], 'peaks':peaks, 'components':source['radio']['components'], 'morphology':morphology}
		entry = {'RGZ':rgz, 'using_contour':using_contour, 'using_peaks':using_peaks}
		if 'SDSS' in source:
			entry['SDSS'] = source['SDSS']
		if 'AllWISE' in source:
			entry['AllWISE'] = source['AllWISE']
		if cluster_w is not None:
			entry['WHL'] = whl_prop
		if cluster_r is not None:
			entry['redMaPPer'] = rm_prop
		if cluster_a is not None:
			entry['AMFDR9'] = amf_prop
		
		return entry

def bending(args, peak_count):
	'''
	Find the sources matching the given search parameters and morphology and run the processing pipeline
	'''
	
	morphology = 'double' if peak_count == 2 else 'triple'
	
	# While loop lasts as long as cursor keeps timing out
	done = False
	while not done:
		
		# Determine which sources have already been processed
		completed = []
		if os.path.exists(completed_file):
			with open(completed_file, 'r') as f:
				lines = f.readlines()
				for line in lines:
					completed.append(int(line))
		
		previous = len(lines)
		args['$and'][0]['catalog_id'] = {'$nin':completed}
		
		# Find the bending and cluster results for each source that matches
		try:
			count = bending_coll.find({'RGZ.morphology':morphology}).count()
			with open(completed_file, 'a') as f:
				for source in catalog.find(args).batch_size(100):
					entry = get_entry(source, peak_count)
					print >> f, source['catalog_id']
					if entry is not None:
						count += 1
						output('%i %s' % (count, source['zooniverse_id']))
						bending_coll.insert(entry)
		
		except pymongo.errors.CursorNotFound as e:
			output(e, logging.exception)
			with open(completed_file, 'r') as f:
				lines = f.readlines()
			new = len(lines) - previous
			if not new:
				raise
		
		# If the try block finishes, everything is processed
		done = True
	
	return count

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
	sdss_keys = ['ra', 'dec', 'objID', 'photo_redshift', 'photo_redshift_err', 'spec_redshift', 'spec_redshift_err', 'u', 'u_err', 'g', 'g_err', 'r', 'r_err', 'i', 'i_err', 'z', 'z_err']
	wise_keys = ['ra', 'dec', 'designation', 'photo_redshift', 'w1mpro', 'w1sigmpro', 'w1snr', 'w2mpro', 'w2sigmpro', 'w2snr', 'w3mpro', 'w3sigmpro', 'w3snr', 'w4mpro', 'w4sigmpro', 'w4snr']
	whl_keys = ['name', 'ra', 'dec', 'zphot', 'zspec', 'N500', 'N500sp', 'RL*500', 'r500', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	#rm_keys = ['name', 'ra', 'dec', 'zlambda', 'zspec', 'S', 'lambda', 'separation_deg', 'separation_Mpc', 'position_angle', 'orientation_contour', 'orientation_peaks']
	bending_keys = ['pos_angle_0', 'pos_angle_1', 'opening_angle', 'bisector', 'tail_deg_0', 'tail_deg_1', 'tail_kpc_0', 'tail_kpc_1', 'ratio_1', 'ratio_0', 'asymmetry']
	
	all_keys = [rgz_keys, sdss_keys, wise_keys, whl_keys, bending_keys, bending_keys]
	dict_names = ['RGZ', 'SDSS', 'AllWISE', 'WHL', 'using_contour', 'using_peaks']
	
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
	
	# When changing parameters, check lines 13, 23, 517, 530, 542, and 724
	
	logging.basicConfig(filename='%s/bending.log' % rgz_path, level=logging.DEBUG, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	logging.captureWarnings(True)
	logging.info('Bending run from command line')
	
	z_range = [0.01, 0.8]
	double_args = {'$and': [{'ignore_bending':False}, \
				  {'$or': [{'radio.number_peaks':2, 'radio.number_components':1}, \
				 		   {'radio.number_components':2}]}, \
				  {'$or': [{'SDSS.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
						   {'SDSS.spec_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
						   {'AllWISE.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}] }]}
	triple_args = {'$and': [{'ignore_bending':False}, \
				  {'$or': [{'radio.number_peaks':3, 'radio.number_components':{'$in':[1,2]}}, \
				 		   {'radio.number_components':3}]}, \
				  {'$or': [{'SDSS.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
						   {'SDSS.spec_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}, \
						   {'AllWISE.photo_redshift':{'$gte':z_range[0], '$lt':z_range[1]}}] }]}
	
	try:
		output('Processing double sources from RGZ')
		output('%i double sources matched to clusters' % bending(double_args, 2))
		output('Processing triple sources from RGZ')
		output('%i triple sources matched to clusters' % bending(triple_args, 3))
		to_file('%s/csv/bending_control_15.csv' % rgz_path, bending_coll)
	except BaseException as e:
		output(e, logging.exception)
		raise
