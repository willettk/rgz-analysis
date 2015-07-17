from pymongo import MongoClient
from collections import Counter
import datetime
import numpy as np
import operator
from scipy.ndimage.filters import maximum_filter
from scipy import stats
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from astropy.io import fits
from astropy import wcs
import time
from scipy.linalg.basic import LinAlgError

start = time.clock()

# Trim this down to just get the consensus answers

client = MongoClient('localhost', 27017)
db = client['radio'] 

main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)
IMG_HEIGHT_OLD = 424.0			# number of pixels in the original JPG image along the y axis
IMG_WIDTH_OLD = 424.0			# number of pixels in the original JPG image along the x axis
IMG_HEIGHT_NEW = 500.0			# number of pixels in the downloaded JPG image along the y axis
IMG_WIDTH_NEW = 500.0			# number of pixels in the downloaded JPG image along the x axis
FITS_HEIGHT = 301.0			# number of pixels in the FITS image (?) along the y axis
FITS_WIDTH = 301.0			# number of pixels in the FITS image (?) along the x axis
FIRST_FITS_HEIGHT = 132.0       # number of pixels in the FITS image along the y axis
FIRST_FITS_WIDTH = 132.0       # number of pixels in the FITS image along the y axis
PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image
xmin = 1.
xmax = IMG_HEIGHT_NEW
ymin = 1.
ymax = IMG_WIDTH_NEW

subjects = db['radio_subjects'] 		# subjects = images
classifications = db['radio_classifications']	# classifications = classifications of each subject per user

# Running as part of loop

N = 0
class_lim = {'$or':[{'state':'complete','classification_count':{'$gt':19}},{'state':'complete','classification_count':{'$gt':4},'metadata.contour_count':1}]}

# Get dictionary for finding the path to FITS files and WCS headers

with open('/Volumes/3TB/rgz/first_fits.txt') as f:
    lines = f.readlines()

pathdict = {}
for l in lines:
    spl = l.split(' ')
    pathdict[spl[1].strip()] = '/Volumes/3TB/rgz/raw_images/RGZ-full.%i/FIRST-IMGS/%s.fits' % (int(spl[0]),spl[1].strip())

radec_writefile = open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/radec_completed_nusers.csv','w')
nosources_writefile = open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/no_sources.csv','w')

print >> radec_writefile,'first_id              ra            dec          n_class     n_users'

start_ind = 0   # Restart list from last completed classification if desired
subjects_list = list(subjects.find(class_lim))[start_ind:]
for sub in subjects_list:

    N += 1
    
    imgid = sub['_id']
    zid = sub['zooniverse_id']
    srcid = sub['metadata']['source']
    if srcid == 'tutorial':
        print 'Tutorial image - skipping due to lack of WCS.'
        continue

    clist = list(classifications.find({"subject_ids": imgid, "updated_at": {"$gt": main_release_date}}))
    
    bad_keys = ('finished_at','started_at','user_agent','lang','pending')
    
    cdict = {}
    checksum_list = []
    
    for c in clist:
        # Want most popular combo for each NUMBER of galaxies
        
        sumlist = []    # List of the checksums over all possible combinations
        goodann = [x for x in c['annotations'] if x.keys()[0] not in bad_keys]
        n_galaxies = len(goodann)
    
        for idx,ann in enumerate(goodann):

            xmaxlist = []
            try:
                radio_comps = ann['radio']
    
                # loop over all the radio components within an galaxy
                if radio_comps != 'No Contours':
                    for rc in radio_comps:
                        xmaxlist.append(float(radio_comps[rc]['xmax']))
                # or make the value -99 if there are no contours
                else:
                    xmaxlist.append(-99)
    
                product = reduce(operator.mul, xmaxlist, 1)
                sumlist.append(round(product,3))
            except KeyError:
                print zid,goodann
    
        checksum = round(sum(sumlist),3)
        checksum_list.append(checksum)
        c['checksum'] = checksum
    
        if cdict.has_key(n_galaxies):
            cdict[n_galaxies].append(checksum)
        else:
            cdict[n_galaxies] = [checksum]
    
    # If all entries contain only annotations for "No Contours" (eg FIRSTJ005717.8+101702)
    # or at least one galaxy with "No Contours" (FIRSTJ131852.2+642707), skip remainder of classification.
    # This should never happen: all radio images should have at least 1 component, but there
    # is at least one example where a single user's (mpetruk) errant classifications were all for "No IR", "No Contours"
    if max(checksum_list) < 0:
        print >> nosources_writefile, imgid,zid,srcid
        continue

    maxval=0
    mc_checksum = 0.
    ngals = 0
    for k,v in cdict.iteritems():
        mc = Counter(v).most_common()
        if mc[0][0] == -99.0:
            if len(mc) > 1:
                mc_best = mc[1]
            else:
                continue
        else:
            mc_best = mc[0]
        if mc_best[1] > maxval:
            maxval = mc_best[1]
            mc_checksum = mc_best[0]
            ngals = k
    
    # Find a galaxy that matches the checksum (easier to keep track as a list)
    
    cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
   
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in bad_keys]
    
    # I know how many galaxies there should be. Need to match the correct IR hosts with radio components.
    # Find the sum of the xmax coordinates for each galaxy. This gives the index to search on.
    
    gal_xmax = {}
    ir_x,ir_y = {},{}
    for k,gal in enumerate(goodann):
        xmax_temp = []
        for v in gal['radio'].itervalues():
            xmax_temp.append(float(v['xmax']))
        checksum2 = round(sum(xmax_temp),3)
        gal_xmax[checksum2] = {}
        gal_xmax[checksum2]['ind'] = k
        gal_xmax[checksum2]['xmax'] = xmax_temp
    
        # Make empty copy of next dict in same loop
        ir_x[k] = []
        ir_y[k] = []
    
    # Now loop over the galaxies themselves
    for c in clist:
        if c['checksum'] == mc_checksum:
    
            annlist = [ann for ann in c['annotations'] if ann.keys()[0] not in bad_keys]
            for ann in annlist:
                if 'ir' in ann.keys():
                    # Find the index k that this corresponds to
                    xmax_checksum = round(sum([float(ann['radio'][a]['xmax']) for a in ann['radio']]),3)
                    k = gal_xmax[xmax_checksum]['ind']
    
                    if ann['ir'] == 'No Sources':
                        ir_x[k].append(-99)
                        ir_y[k].append(-99)
                    else:
                        # Only takes the first IR source right now; NEEDS TO BE MODIFIED.
    
                        ir_x[k].append(float(ann['ir']['0']['x']))
                        ir_y[k].append(float(ann['ir']['0']['y']))
    
    # Perform a kernel density estimate on the data for each galaxy
    
    scale_ir = IMG_HEIGHT_NEW/float(IMG_HEIGHT_OLD)
    
    peak_data = []
    for x,y in zip(ir_x.itervalues(),ir_y.itervalues()):
        foo = [(xt,yt) for xt,yt in zip(x,y) if xt != -99.0 and yt != -99.0]
        x_exists = [f[0]*scale_ir for f in foo]
        y_exists = [f[1]*scale_ir for f in foo]
    
        pd = {}
    
        # KDE needs at least 3 unique points to run density estimator
        if len(set(x_exists)) > 2:
            X, Y = np.mgrid[xmin:xmax, ymin:ymax]
            positions = np.vstack([X.ravel(), Y.ravel()])
            values = np.vstack([x_exists, y_exists])
            try:
                kernel = stats.gaussian_kde(values)
                Z = np.reshape(kernel(positions).T, X.shape)
                
                # Find the number of peaks
                # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array
                
                #neighborhood = generate_binary_structure(2,2)
                neighborhood = np.ones((10,10))
                local_max = maximum_filter(Z, footprint=neighborhood)==Z
                background = (Z==0)
                eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
                detected_peaks = local_max - eroded_background
            
                npeaks = detected_peaks.sum()
    
            except LinAlgError:
                print 'LinAlgError - only %i non-unique IR peaks labeled for %s' % (len(set(x_exists)),zid)
                continue

            #return X,Y,Z,npeaks
    
            pd['X'] = X
            pd['Y'] = Y
            pd['Z'] = Z
            pd['npeaks'] = npeaks
            peak_data.append(pd)
        
            if npeaks > 0:
                # Find the peak
                xpeak = float(X[Z==Z.max()][0])
                ypeak = float(Y[Z==Z.max()][0])

                # Convert the pixel coordinates into RA,dec using the WCS object from the header

                hdulist = fits.open(pathdict[srcid])
                w = wcs.WCS(hdulist[0].header)

                first_ir_scale_x = FIRST_FITS_WIDTH / IMG_WIDTH_NEW
                first_ir_scale_y = FIRST_FITS_HEIGHT / IMG_HEIGHT_NEW
                pixcrd = np.array([[xpeak * first_ir_scale_x,(IMG_HEIGHT_NEW-ypeak) * first_ir_scale_y]],np.float_)
                world = w.wcs_pix2world(pixcrd,0)

                # Write to file

                print >> radec_writefile, srcid, world[0][0],world[0][1], len(x_exists), len(clist)

        else:
    
            pd['npeaks'] = 0
            peak_data.append(pd)

    # Check progress by printing to screen every 100 classifications
    if not N % 100:
        print N+start_ind, datetime.datetime.now().strftime('%H:%M:%S.%f')


radec_writefile.close()
nosources_writefile.close()

import time

end = time.clock()
print 'Total time: %.3f seconds' % (end - start)
print 'Total time per galaxy (%i subjects): %.3f seconds' % (len(subjects_list),(end - start)/len(subjects_list))
