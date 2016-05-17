from __future__ import division

# Local RGZ modules

import collinearity
from load_contours import get_contours,make_pathdict

# Default packages

import datetime
import operator
from collections import Counter
import cStringIO
import urllib
import json
import os.path
import time
import shutil

# Other packages

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.path import Path
import matplotlib.patches as patches

from scipy.ndimage.filters import maximum_filter
from scipy import stats
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.linalg.basic import LinAlgError

from astropy.io import fits
from astropy import wcs

from pymongo import MongoClient

from PIL import Image

# MongoDB parameters

client = MongoClient('localhost', 27017)
db = client['radio'] 

subjects = db['radio_subjects'] 		# subjects = images
classifications = db['radio_classifications']	# classifications = classifications of each subject per user

# General variables for the RGZ sample

main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)
IMG_HEIGHT_OLD = 424.0			# number of pixels in the original JPG image along the y axis
IMG_WIDTH_OLD = 424.0			# number of pixels in the original JPG image along the x axis
IMG_HEIGHT_NEW = 500.0			# number of pixels in the downloaded JPG image along the y axis
IMG_WIDTH_NEW = 500.0			# number of pixels in the downloaded JPG image along the x axis
FITS_HEIGHT = 301.0			# number of pixels in the FITS image (?) along the y axis
FITS_WIDTH = 301.0			# number of pixels in the FITS image (?) along the x axis
FIRST_FITS_HEIGHT = 132.0       # number of pixels in the FITS image along the y axis
FIRST_FITS_WIDTH = 132.0       # number of pixels in the FITS image along the y axis

# Need to add parameters for ATLAS, both IR and radio.

PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image
xmin = 1.
xmax = IMG_HEIGHT_NEW
ymin = 1.
ymax = IMG_WIDTH_NEW

bad_keys = ('finished_at','started_at','user_agent','lang','pending')

expert_names = [u'42jkb', u'ivywong', u'stasmanian', u'klmasters', u'Kevin', u'akapinska', u'enno.middelberg', u'xDocR', u'vrooje', u'KWillett', u'DocR']

# Paths

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
pathdict = make_pathdict()

# Find the consensus classification for a single subject

@profile
def checksum(zid='ARG000255x',experts_only=False,excluded=[],no_anonymous=False,write_peak_data=False):

    # Find the consensus for all users who have classified a particular galaxy

    sub = subjects.find_one({'zooniverse_id':zid})
    imgid = sub['_id']

    # Classifications for this subject after launch date
    class_params = {"subject_ids": imgid, "updated_at": {"$gt": main_release_date}}
    # Only get the consensus classification for the science team members
    if experts_only:
        class_params['expert'] = True

    # If comparing a particular volunteer (such as an expert), don't include self-comparison
    if len(excluded) > 0:
        class_params['user_name'] = {"$nin":excluded}

    '''
    # To exclude the experts:
    class_params['expert'] = {"$exists":False}
    '''

    # To exclude anonymous classifications (registered users only):
    if no_anonymous:
        if class_params.has_key('user_name'):
            class_params['user_name']["$exists"] = True
        else:
            class_params['user_name'] = {"$exists":True}

    _c = classifications.find(class_params)
    #clist_all = list(classifications.find(class_params))

    # Empty dicts and lists 
    cdict = {}
    checksum_list = []

    unique_users = set()
    
    #clen_start = len(clist_all)
    clen_start = 0
    listcount = []

    # Compute the most popular combination for each NUMBER of galaxies identified in image

    clist_all = []
    #for c in clist_all:
    for c in _c:

        clist_all.append(c)

        clen_start += 1

        # Skip classification if they already did one?

        try:
            user_name = c['user_name']
        except KeyError:
            user_name = 'Anonymous'

        if user_name not in unique_users or user_name is 'Anonymous':

            unique_users.add(user_name)
            listcount.append(True)
        
            sumlist = []    # List of the checksums over all possible combinations

            # Only find data that was an actual marking, not metadata
            goodann = [x for x in c['annotations'] if (x.keys()[0] not in bad_keys)]
            n_galaxies = len(goodann)
    
            if n_galaxies > 0:  # There must be at least one galaxy!
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
                    except KeyError:
                        xmaxlist.append(-99)
    
                    # To create a unique ID for the combination of radio components,
                    # take the product of all the xmax coordinates and sum them together.
                    product = reduce(operator.mul, xmaxlist, 1)
                    sumlist.append(round(product,3))

                checksum = sum(sumlist)
            else:
                checksum = -99

            checksum_list.append(checksum)
            c['checksum'] = checksum
    
            # Insert checksum into dictionary with number of galaxies as the index
            if cdict.has_key(n_galaxies):
                cdict[n_galaxies].append(checksum)
            else:
                cdict[n_galaxies] = [checksum]

        else:
            listcount.append(False)
            checksum_list.append(-99)
            #print 'Removing classification for %s' % user_name
    
    # Remove duplicates and classifications for no object
    #clist = [c for lc,c in zip(listcount,checksum_list) if lc and c != -99]
    clist = [c for lc,c in zip(listcount,clist_all) if lc and c['checksum'] != -99]

    clen_diff = clen_start - len(clist)

    '''
    if clen_diff > 0:
        print '\nSkipping %i duplicated classifications for %s. %i good classifications total.' % (clen_diff,zid,len(clist))
    '''

    maxval=0
    mc_checksum = 0.

    # Find the number of galaxies that has the highest number of consensus classifications

    for k,v in cdict.iteritems():
        mc = Counter(v).most_common()
        # Check if the most common selection coordinate was for no radio contours
        if mc[0][0] == -99.0:
            if len(mc) > 1:
                # If so, take the selection with the next-highest number of counts
                mc_best = mc[1]
            else:
                continue
        # Selection with the highest number of counts
        else:
            mc_best = mc[0]
        # If the new selection has more counts than the previous one, choose it as the best match;
        # if tied or less than this, remain with the current consensus number of galaxies
        if mc_best[1] > maxval:
            maxval = mc_best[1]
            mc_checksum = mc_best[0]
    
    # Find a galaxy that matches the checksum (easier to keep track as a list)

    try:
        cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
    except StopIteration:
        # Necessary for objects like ARG0003par; one classifier recorded 22 "No IR","No Contours" in a short space. Still shouldn't happen.
        print 'No non-zero classifications recorded for %s' % zid
        return None
    '''

    try:
        index = clist.index(mc_checksum)
        cmatch = _c[index]
    except ValueError:
        # Necessary for objects like ARG0003par; one classifier recorded 22 "No IR","No Contours" in a short space. Still shouldn't happen.
        print 'No non-zero classifications recorded for %s' % zid
        return None
   
    '''
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in bad_keys]

    # Find the sum of the xmax coordinates for each galaxy. This gives the index to search on.
    
    cons = {}
    cons['zid'] = zid
    cons['source'] = sub['metadata']['source']
    ir_x,ir_y = {},{}
    cons['answer'] = {}
    cons['n_users'] = maxval
    cons['n_total'] = len(clist)

    answer = cons['answer']

    for k,gal in enumerate(goodann):
        xmax_temp = []
        bbox_temp = []
        try:
            for v in gal['radio'].itervalues():
                xmax_temp.append(float(v['xmax']))
                bbox_temp.append((v['xmax'],v['ymax'],v['xmin'],v['ymin']))
            checksum2 = round(sum(xmax_temp),3)
            answer[checksum2] = {}
            answer[checksum2]['ind'] = k
            answer[checksum2]['xmax'] = xmax_temp
            answer[checksum2]['bbox'] = bbox_temp
        except KeyError:
            print gal, zid
        except AttributeError:
            print 'No Sources, No IR recorded for %s' % zid
    
        # Make empty copy of next dict in same loop
        ir_x[k] = []
        ir_y[k] = []
    
    # Now loop over all sets of classifications to get the IR counterparts
    for c in clist:
        if c['checksum'] == mc_checksum:
    
            annlist = [ann for ann in c['annotations'] if ann.keys()[0] not in bad_keys]
            for ann in annlist:
                if 'ir' in ann.keys():
                    # Find the index k that this corresponds to
                    try:
                        xmax_checksum = round(sum([float(ann['radio'][a]['xmax']) for a in ann['radio']]),3)
                    except TypeError:
                        xmax_checksum = -99

                    try:
                        k = answer[xmax_checksum]['ind']

                        if ann['ir'] == 'No Sources':
                            ir_x[k].append(-99)
                            ir_y[k].append(-99)
                        else:
                            # Only takes the first IR source right now; NEEDS TO BE MODIFIED.

                            ir_x[k].append(float(ann['ir']['0']['x']))
                            ir_y[k].append(float(ann['ir']['0']['y']))
                    except KeyError:
                        print '"No radio" still appearing as valid consensus option.'

    # Perform a kernel density estimate on the data for each galaxy
    
    scale_ir = IMG_HEIGHT_NEW/IMG_HEIGHT_OLD

    peak_data = []

    # Remove empty IR peaks if they exist

    for (xk,xv),(yk,yv) in zip(ir_x.iteritems(),ir_y.iteritems()):
        
        if len(xv) == 0:
            ir_x.pop(xk)
        if len(yv) == 0:
            ir_y.pop(yk)

    assert len(ir_x) == len(ir_y),'Lengths of ir_x (%i) and ir_y (%i) are not the same' % (len(ir_x),len(ir_y))

    for (xk,xv),(yk,yv) in zip(ir_x.iteritems(),ir_y.iteritems()):
        
        if len(xv) == 0:
            irx

        pd = {}
    
        x_exists = [xt * scale_ir for xt in xv if xt != -99.0]
        y_exists = [yt * scale_ir for yt in yv if yt != -99.0]

        x_all = [xt * scale_ir for xt in xv]
        y_all = [yt * scale_ir for yt in yv]
        coords_all = [(xx,yy) for xx,yy in zip(x_all,y_all)]
        ir_Counter = Counter(coords_all)
        most_common_ir = ir_Counter.most_common(1)[0][0]

        if len(Counter(x_exists)) > 2 and len(Counter(y_exists)) > 2 and most_common_ir != (-99,-99):

            # X,Y = grid of uniform coordinates over the IR pixel plane
            X, Y = np.mgrid[xmin:xmax, ymin:ymax]
            positions = np.vstack([X.ravel(), Y.ravel()])
            try:
                values = np.vstack([x_exists, y_exists])
            except ValueError:
                # Breaks on the tutorial subject. Find out why len(x) != len(y)
                print zid
                print 'Length of IR x array: %i; Length of IR y array: %i' % (len(x_exists),len(y_exists))
            try:
                kernel = stats.gaussian_kde(values)
            except LinAlgError:
                print 'LinAlgError in KD estimation for %s' % zid,x_exists,y_exists
                continue

            # Even if there are more than 2 sets of points, if they are mutually co-linear, 
            # matrix can't invert and kernel returns NaNs. 

            kp = kernel(positions)

            if np.isnan(kp).sum() > 0:
                acp = collinearity.collinear(x_exists,y_exists)
                if len(acp) > 0:
                    print 'There are %i unique points for %s (source no. %i in the field), but all are co-linear; KDE estimate does not work.' % (len(Counter(x_exists)),zid,xk)
                else:
                    print 'There are NaNs in the KDE for %s (source no. %i in the field), but points are not co-linear.' % (zid,xk)

                for k,v in answer.iteritems():
                    if v['ind'] == xk:
                        answer[k]['ir'] = (np.mean(x_exists),np.mean(y_exists))
        
            else:

                Z = np.reshape(kp.T, X.shape)
                
                # Find the number of peaks
                # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array
                
                neighborhood = np.ones((10,10))
                local_max = maximum_filter(Z, footprint=neighborhood)==Z
                background = (Z==0)
                eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
                detected_peaks = local_max ^ eroded_background
                
                npeaks = detected_peaks.sum()
    
                #return X,Y,Z,npeaks
    
                pd['X'] = X
                pd['Y'] = Y
                pd['Z'] = Z
                pd['npeaks'] = npeaks

                try:
                    xpeak = float(pd['X'][pd['Z']==pd['Z'].max()][0])
                    ypeak = float(pd['Y'][pd['Z']==pd['Z'].max()][0])
                except IndexError:
                    print pd
                    print zid, clist

                for k,v in answer.iteritems():
                    if v['ind'] == xk:
                        answer[k]['ir_peak'] = (xpeak,ypeak)
                        # Don't write to consensus for serializable JSON object 
                        if write_peak_data:
                            answer[k]['peak_data'] = pd
                            answer[k]['ir_x'] = x_exists
                            answer[k]['ir_y'] = y_exists
        else:

            # Note: need to actually put a limit in if less than half of users selected IR counterpart.
            # Right now it still IDs a sources even if only 1/10 users said it was there.

            for k,v in answer.iteritems():
                if v['ind'] == xk:
                    # Case 1: multiple users selected IR source, but not enough unique points to pinpoint peak
                    if most_common_ir != (-99,-99) and len(x_exists) > 0 and len(y_exists) > 0:
                        answer[k]['ir'] = (x_exists[0],y_exists[0])
                    # Case 2: most users have selected No Sources
                    else:
                        answer[k]['ir'] = (-99,-99)

    return cons
    
if __name__ == "__main__":
    checksum()
