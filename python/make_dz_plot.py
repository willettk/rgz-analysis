from bending_analysis import *
import csv
from scipy.stats import linregress

def get_cluster_match2(source):
	'''
	Given a source from RGZ, match it to a cluster and calculate the redshift-dependent bending parameters
	'''
	
	ir = coord.SkyCoord(source['SDSS']['ra'], source['SDSS']['dec'], unit=(u.deg,u.deg), frame='icrs') if 'SDSS' in source else \
		 coord.SkyCoord(source['AllWISE']['ra'], source['AllWISE']['dec'], unit=(u.deg,u.deg), frame='icrs')
	z, z_err = get_z(source)
	
	# Match to cluster catalogs
	cluster_w = get_whl(ir, z, z_err, 15, 1)
	whl_prop = {}
	if cluster_w is not None:
		c_pos = coord.SkyCoord(cluster_w['RAdeg'], cluster_w['DEdeg'], unit=(u.deg,u.deg), frame='icrs')
		c_sep_arc = c_pos.separation(ir)
		c_sep_mpc = float(cosmo.angular_diameter_distance(cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot'])/u.Mpc * c_sep_arc.to(u.rad)/u.rad)
		c_pos_angle = c_pos.position_angle(ir)
		whl_prop = {'ra':c_pos.ra.deg, 'dec':c_pos.dec.deg, 'separation_deg':c_sep_arc.deg, 'separation_Mpc':c_sep_mpc, 'position_angle':c_pos_angle.wrap_at(2*np.pi*u.rad).deg, 'r/r500':c_sep_mpc/cluster_w['r500'], 'M500':np.exp(1.08*np.log(cluster_w['RL*500'])-1.37), 'zbest':cluster_w['zspec'] if 'zspec' in cluster_w else cluster_w['zphot']}
		for key in ['_id', 'N500', 'N500sp', 'RL*500', 'name', 'r500', 'zphot', 'zspec']:
			if key in cluster_w:
				whl_prop[key] = cluster_w[key]
		
		source['WHL'] = whl_prop
		return source

def print_data(filename='dz_data.csv'):
	if not os.path.isfile(filename) or os.path.getsize(filename)==0: # Add header is file is empty
		with open(filename, 'w') as f:
			print >> f, 'rgz_id,z_gal,z_whl'
	
	rgz_ids, z_gal, z_whl = load_data()
	
	with open(filename, 'a') as f:
		params = {'RGZ.RGZ_id':{'$nin':list(rgz_ids)}, 'RGZ.overedge': 0, 'RGZ.radio_consensus': {'$gte': 0.65}, 'using_peaks.size_deg': {'$lte': 0.025}}
		for source in bent_sources.find(params).batch_size(50):
			match = get_cluster_match2(source)
			if match is not None:
				print >> f, '%i,%f,%f' % (match['RGZ']['RGZ_id'], match['best']['redshift'], match['WHL']['zbest'])

def load_data(filename='dz_data.csv'):
	rgz_ids, z_gal, z_whl = [], [], []
	with open(filename, 'r') as f:
		reader = csv.reader(f, delimiter=',')
		next(reader, None)
		for row in reader:
			#rgz_ids.append(int(row[0]))
			z_gal.append(float(row[0]))
			z_whl.append(float(row[1]))
	
	return np.array(rgz_ids), np.array(z_gal), np.array(z_whl)

def make_plot():
	#rgz_ids, z_gal, z_whl = load_data()
	#dz = (z_whl-z_gal) / (1+z_gal)
	
	null, z, dz = load_data('/home/garon/Documents/RGZdata/bending/dz_whl_with_z.csv')
	dz = -dz / (1+z)
	
	cut = 0.04
	n, bins, patches = plt.hist(dz[np.abs(dz)<.25], bins=30)
	#plt.xlim([-2.5*cut,2.5*cut])
	
	plt.rc('text', usetex=True)
	plt.xlabel(r'$\Delta$z / (1+z)')
	plt.ylabel('Count')
	plt.axvline(-cut, color='r', label=r'$\Delta$z cut')
	plt.axvline(cut, color='r')
	
	bin_centers = (bins[:-1]+bins[1:])/2.
	wings = np.abs(bin_centers)>cut
	wing_pos = bin_centers[wings]
	wing_n = n[wings]
	popt, pcov = curve_fit(lambda x,c: c, wing_pos, wing_n)
	plt.axhline(popt[0], color='k', label='Background')
	
	false_rejected = sum(wing_n[wing_n>popt[0]]-popt[0])
	false_approved = popt[0] * len(n[np.logical_not(wings)])
	correct_rejected = sum(wing_n) - false_rejected
	correct_approved = sum(n[np.logical_not(wings)]) - false_approved
	completion = 1 - (false_rejected / correct_approved)
	contamination = false_approved / correct_approved
	
	plt.plot(0, 0, c='w', label='Completion: %.3f\nContamination: %.3f' % (completion, contamination))
	#plt.title('Bending sample cross-match')
	plt.title('Full RGZ-WHL cross-match')
	plt.legend()

def using_kde():
	null, z, dz = load_data('/home/garon/Documents/RGZdata/bending/dz_whl_with_z.csv')
	dz = -dz / (1+z)
	
	from scipy.stats import gaussian_kde as g_kde
	kernel = g_kde(dz)
	z_range = np.arange(np.min(dz), np.max(dz), (np.max(dz)-np.min(dz))/10000.)
	dist = kernel(z_range)
	plt.plot(z_range, dist, label='Full')

	crop = np.abs(dz)<0.1
	kernel_crop = g_kde(dz[crop])
	z_range_crop = np.arange(np.min(dz[crop]), np.max(dz[crop]), (np.max(dz[crop])-np.min(dz[crop]))/1000.)
	dist_crop = kernel_crop(z_range_crop)
	plt.plot(z_range_crop, dist_crop, label='Crop')

	plt.legend()

