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

import pymongo 
from pymongo import MongoClient

from PIL import Image

# MongoDB parameters

client = MongoClient('localhost', 27017)
db = client['radio'] 

subjects = db['radio_subjects'] # subjects = images
classifications = db['radio_classifications']# classifications = classifications of each subject per user

# Create index on subject IDs so that queries run faster

subindex = classifications.create_index([('subject_ids',pymongo.ASCENDING)],name='subject_ids_1')

# General variables for the RGZ sample
# Need to add separate parameters for ATLAS vs FIRST, both IR and radio.

main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)

img_params = {
    'first':{
        'IMG_HEIGHT_OLD':424.0    ,         # number of pixels in the IR coordinate system (as recorded in Mongo) along the y axis
        'IMG_WIDTH_OLD':424.0     ,         # number of pixels in the IR coordinate system (as recorded in Mongo) along the x axis
        'IMG_HEIGHT_NEW':500.0    ,         # number of pixels in the downloaded JPG image along the y axis
        'IMG_WIDTH_NEW':500.0     ,         # number of pixels in the downloaded JPG image along the x axis
        'FITS_HEIGHT':132.0       ,         # number of pixels in the FITS image along the y axis (radio only)
        'FITS_WIDTH':132.0        ,         # number of pixels in the FITS image along the y axis (radio only)
        'PIXEL_SIZE':1.3748                 # the number of arcseconds per pixel in the radio FITS image (from CDELT1)
    },
    'atlas':{
        'IMG_HEIGHT_OLD':424.0    ,         # number of pixels in the IR coordinate system (as recorded in Mongo) along the y axis
        'IMG_WIDTH_OLD':424.0     ,         # number of pixels in the IR coordinate system (as recorded in Mongo) along the x axis
        'IMG_HEIGHT_NEW':500.0    ,         # number of pixels in the downloaded PNG image along the y axis
        'IMG_WIDTH_NEW':500.0     ,         # number of pixels in the downloaded PNG image along the x axis
        'FITS_HEIGHT':201.0       ,         # number of pixels in the FITS image along the y axis (both IR and radio)
        'FITS_WIDTH':201.0        ,         # number of pixels in the FITS image along the x axis (both IR and radio)
        'PIXEL_SIZE':0.6000                 # the number of arcseconds per pixel in the radio FITS image (from CDELT1)
    }
}

bad_keys = ('finished_at','started_at','user_agent','lang','pending')

expert_names = [u'42jkb', u'ivywong', u'stasmanian', u'klmasters', u'Kevin', u'akapinska', u'enno.middelberg', u'xDocR', u'vrooje', u'KWillett', u'DocR']

# Paths

paths = ('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis','/data/tabernacle/larry/RGZdata/rgz-analysis')
for path in paths:
    if os.path.exists(path):
        rgz_dir = path
if rgz_dir == None:
    print "Unable to find the hardcoded local path to store outputs."

pathdict = make_pathdict()

# Paths v2 - from RGZcatalogCode

def determine_paths(paths):

    found_path = False
    for path in paths:
        if os.path.exists(path):
            found_path = True
            return path

    if found_path == False:
        print "Unable to find the hardcoded local path:"
        print paths
        return None

rgz_path = determine_paths(('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis',
                           '/data/tabernacle/larry/RGZdata/rgz-analysis'))
data_path = determine_paths(('/Volumes/REISEPASS/','/data/extragal/willett'))
plot_path = "{0}/rgz/plots".format(data_path)

# Find the consensus classification for a single subject

def checksum(zid,experts_only=False,excluded=[],no_anonymous=False,include_peak_data=True):

    # Find the consensus for all users who have classified a particular galaxy

    sub = subjects.find_one({'zooniverse_id':zid})
    imgid = sub['_id']
    survey = sub['metadata']['survey']

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

    # Empty dicts and lists 
    cdict = {}

    unique_users = set()
    
    clen_start = 0
    clist_all = []
    listcount = []

    # Compute the most popular combination for each NUMBER of galaxies identified in image
    
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

            c['checksum'] = checksum
    
            # Insert checksum into dictionary with number of galaxies as the index
            if cdict.has_key(n_galaxies):
                cdict[n_galaxies].append(checksum)
            else:
                cdict[n_galaxies] = [checksum]

        else:
            listcount.append(False)
    
    # Remove duplicates and classifications for no object
    clist = [c for lc,c in zip(listcount,clist_all) if lc and c['checksum'] != -99]

    clen_diff = clen_start - len(clist)

    '''
    if clen_diff > 0:
        print '\nSkipping {0:d} duplicated classifications for {1}. {2:d} good classifications total.'.format(clen_diff,zid,len(clist))
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
        print 'No non-zero classifications recorded for {0}'.format(zid)
        return None
   
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in bad_keys]

    # Find the sum of the xmax coordinates for each galaxy. This gives the index to search on.
    
    cons = {}
    cons['zid'] = zid
    cons['source'] = sub['metadata']['source']
    cons['survey'] = survey
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
            print 'No Sources, No IR recorded for {0}'.format(zid)
    
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
    
    scale_ir = img_params[survey]['IMG_HEIGHT_NEW'] * 1./img_params[survey]['IMG_HEIGHT_OLD']

    peak_data = []

    # Remove empty IR peaks if they exist

    for (xk,xv),(yk,yv) in zip(ir_x.iteritems(),ir_y.iteritems()):
        
        if len(xv) == 0:
            ir_x.pop(xk)
        if len(yv) == 0:
            ir_y.pop(yk)

    assert len(ir_x) == len(ir_y),'Lengths of ir_x ({0:d}) and ir_y ({1:d}) are not the same'.format(len(ir_x),len(ir_y))

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

            xmin = 1.
            xmax = img_params[survey]['IMG_HEIGHT_NEW']
            ymin = 1.
            ymax = img_params[survey]['IMG_WIDTH_NEW']

            # X,Y = grid of uniform coordinates over the IR pixel plane
            X, Y = np.mgrid[xmin:xmax, ymin:ymax]
            positions = np.vstack([X.ravel(), Y.ravel()])
            try:
                values = np.vstack([x_exists, y_exists])
            except ValueError:
                # Breaks on the tutorial subject. Find out why len(x) != len(y)
                print zid
                print 'Length of IR x array: {0:d}; Length of IR y array: {1:d}'.format(len(x_exists),len(y_exists))
            try:
                kernel = stats.gaussian_kde(values)
            except LinAlgError:
                print 'LinAlgError in KD estimation for {0}'.format(zid,x_exists,y_exists)
                continue

            # Even if there are more than 2 sets of points, if they are mutually co-linear, 
            # matrix can't invert and kernel returns NaNs. 

            kp = kernel(positions)

            if np.isnan(kp).sum() > 0:
                acp = collinearity.collinear(x_exists,y_exists)
                if len(acp) > 0:
                    print 'There are {0:d} unique points for {1} (source no. {2:d} in the field), but all are co-linear; KDE estimate does not work.'.format(len(Counter(x_exists)),zid,xk)
                else:
                    print 'There are NaNs in the KDE for {0} (source no. {1:d} in the field), but points are not co-linear.'.format(zid,xk)

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
                        if include_peak_data:
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
    
def one_answer(zid,user_name):

    # Find the result for one user who did a particular galaxy

    sub = subjects.find_one({'zooniverse_id':zid})
    imgid = sub['_id']

    # Classifications for this subject after launch date
    class_params = {"subject_ids": imgid, "updated_at": {"$gt": main_release_date},'user_name':user_name}
    clist = list(classifications.find(class_params))
  
    # Empty dicts and lists 
    cdict = {}
    
    for c in clist:
        # Want most popular combination for each NUMBER of galaxies identified in image
        
        sumlist = []    # List of the checksums over all possible combinations
        # Only find data that was an actual marking, not metadata
        goodann = [x for x in c['annotations'] if x.keys()[0] not in bad_keys]
        n_galaxies = len(goodann)
    
        for idx,ann in enumerate(goodann):
    
            xmaxlist = []
            radio_comps = ann['radio']
    
            # loop over all the radio components within an galaxy
            if radio_comps != 'No Contours':
                for rc in radio_comps:
                    xmaxlist.append(float(radio_comps[rc]['xmax']))
            # or make the value -99 if there are no contours
            else:
                xmaxlist.append(-99)
    
            # To create a unique ID for the combination of radio components,
            # take the product of all the xmax coordinates and sum them together.
            product = reduce(operator.mul, xmaxlist, 1)
            sumlist.append(round(product,3))
    
        checksum = round(sum(sumlist),3)
        c['checksum'] = checksum
    
        # Insert checksum into dictionary with number of galaxies as the index
        if cdict.has_key(n_galaxies):
            cdict[n_galaxies].append(checksum)
        else:
            cdict[n_galaxies] = [checksum]
    
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
        # Crude way to check for No Sources and No Contours (mc_checksum = 0.)
        cons = {'zid':zid,'answer':{}}
        return cons
   
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in bad_keys]

    # Find the sum of the xmax coordinates for each galaxy. This gives the index to search on.
    
    cons = {}
    cons['zid'] = zid
    cons['answer'] = {}
    cons['n_users'] = 1
    cons['n_total'] = 1
    answer = cons['answer']

    ir_x,ir_y = {},{}
    for k,gal in enumerate(goodann):
        xmax_temp = []
        try:
            for v in gal['radio'].itervalues():
                xmax_temp.append(float(v['xmax']))
        except AttributeError:
            xmax_temp.append(-99)

        checksum2 = round(sum(xmax_temp),3)
        answer[checksum2] = {}
        answer[checksum2]['ind'] = k
        answer[checksum2]['xmax'] = xmax_temp
    
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
                    try:
                        xmax_checksum = round(sum([float(ann['radio'][a]['xmax']) for a in ann['radio']]),3)
                    except TypeError:
                        xmax_checksum = -99
                    k = answer[xmax_checksum]['ind']
    
                    if ann['ir'] == 'No Sources':
                        ir_x[k].append(-99)
                        ir_y[k].append(-99)
                    else:
                        # Only takes the first IR source right now; NEEDS TO BE MODIFIED.
    
                        ir_x[k].append(float(ann['ir']['0']['x']))
                        ir_y[k].append(float(ann['ir']['0']['y']))


                    for k,v in answer.iteritems():
                        if v['ind'] == k:
                            answer[k]['ir_peak'] = (xpeak,ypeak)
        
    return cons

def grab_image(subject,imgtype='standard'):

    # Import a JPG from the RGZ subjects. Try to find a local version before downloading over the web
    
    url = subject['location'][imgtype]
    filename = "{0}/{1}/{2}".format(data_path,imgtype,url.split('/')[-1])

    if os.path.exists(filename):
        with open(filename) as f:
            im = Image.open(f)
            im.load()
    else:
        im = Image.open(cStringIO.StringIO(urllib.urlopen(url).read()))

    return im

def plot_consensus(consensus,figno=1,savefig=False):

    # Plot 4-panel image of IR, radio, KDE estimate, and consensus
    
    zid = consensus['zid']
    answer = consensus['answer']
    sub = subjects.find_one({'zooniverse_id':zid})
    survey = sub['metadata']['survey']

    # Get contour data
    contours = get_contours(sub,pathdict)
    
    # Key bit that sets the difference between surveys.
    #   contours['width'] = img_params[survey]['FITS_WIDTH']
    #   contours['height'] = img_params[survey]['FITS_HEIGHT']

    sf_x = img_params[survey]['IMG_WIDTH_NEW'] * 1./contours['width']
    sf_y = img_params[survey]['IMG_HEIGHT_NEW'] * 1./contours['height']
    
    verts_all = []
    codes_all = []
    components = contours['contours']

    for comp in components:
    
        # Order of bounding box components is (xmax,ymax,xmin,ymin)
        comp_xmax,comp_ymax,comp_xmin,comp_ymin = comp[0]['bbox']
        
        # Only plot radio components identified by the users as the consensus;
        # check on the xmax value to make sure
        for v in answer.itervalues():
            if comp_xmax in v['xmax']:
    
                for idx,level in enumerate(comp):
                    verts = [((p['x'])*sf_x,(p['y']-1)*sf_y) for p in level['arr']]
                    
                    codes = np.ones(len(verts),int) * Path.LINETO
                    codes[0] = Path.MOVETO
                
                    verts_all.extend(verts)
                    codes_all.extend(codes)
    
    try:
        path = Path(verts_all, codes_all)
        patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)
    except AssertionError:
        print 'Users found no components for consensus match of {0}'.format(zid)
    
    # Plot the infrared results
    
    fig = plt.figure(figno,(15,4))
    fig.clf()
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)
    
    colormaparr = [cm.hot_r,cm.Blues,cm.RdPu,cm.Greens,cm.PuBu,cm.YlGn,cm.Greys][::-1]
    colorarr = ['r','b','m','g','c','y','k'][::-1]
    
    # If, in the rare case, that the consensus has more unique sources than the number of colors:
    if len(answer) > len(colorarr):
        colorarr *= int(len(answer)/len(colorarr))+1
        colormaparr *= int(len(answer)/len(colorarr))+1
    
    if len(answer) > 0: # At least one galaxy was identified
        for idx,ans in enumerate(answer.itervalues()):

            if ans.has_key('peak_data'):

                xmin = 1.
                xmax = img_params[survey]['IMG_HEIGHT_NEW']
                ymin = 1.
                ymax = img_params[survey]['IMG_WIDTH_NEW']

                # Plot the KDE map
                colormap = colormaparr.pop()
                ax3.imshow(np.rot90(ans['peak_data']['Z']), cmap=colormap,extent=[xmin, xmax, ymin, ymax])
        
                # Plot individual sources
                color = colorarr.pop()
                x_plot,y_plot = ans['ir_x'],ans['ir_y']
                ax3.scatter(x_plot, y_plot, c=color, marker='o', s=10, alpha=1./len(x_plot))
                ax4.plot([ans['ir_peak'][0]],[ans['ir_peak'][1]],color=color,marker='*',markersize=12)
    
            elif ans.has_key('ir'):
                color = colorarr.pop()
                x_plot,y_plot = ans['ir']
                ax3.plot([x_plot],[y_plot],color=color,marker='o',markersize=2)
                ax4.plot([x_plot],[y_plot],color=color,marker='*',markersize=12)
            else:
                ax4.text(img_params[survey]['IMG_WIDTH_NEW']+50,idx*25,'#{0:d} - no IR host'.format(idx),fontsize=11)

    ax3.set_xlim([0, img_params[survey]['IMG_WIDTH_NEW']])
    ax3.set_ylim([img_params[survey]['IMG_HEIGHT_NEW'], 0])
    ax3.set_title(zid)
    ax3.set_aspect('equal')
    
    ax4.set_xlim([0, img_params[survey]['IMG_WIDTH_NEW']])
    ax4.set_ylim([img_params[survey]['IMG_HEIGHT_NEW'], 0])
    ax4.set_title('Consensus ({0:d}/{1:d} users)'.format(consensus['n_users'],consensus['n_total']))
    
    ax4.set_aspect('equal')
    
    # Display IR and radio images
    
    im_standard = grab_image(sub,imgtype='standard')

    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax1.set_title('WISE')

    im_radio = grab_image(sub,imgtype='radio')

    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    ax2.set_title(sub['metadata']['source'])
    ax2.get_yaxis().set_ticklabels([])

    ax3.get_yaxis().set_ticklabels([])

    # Plot contours identified as the consensus
    if len(answer) > 0:
        ax4.add_patch(patch_black)
    ax4.yaxis.tick_right()

    nticks = 5
    
    ax1.get_xaxis().set_ticks(np.arange(nticks)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)
    ax2.get_xaxis().set_ticks(np.arange(nticks)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)
    ax3.get_xaxis().set_ticks(np.arange(nticks)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)
    ax4.get_xaxis().set_ticks(np.arange(nticks+1)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)

    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    if savefig == True:
        fig.savefig('{0}/{1}/{2}.pdf'.format(plot_path,consensus['survey'],zid))
    else:
        plt.show()

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    plt.close()

    return None

def check_class(zid):

    # Print list of users who have classified a particular subject

    sid = subjects.find_one({'zooniverse_id':zid})['_id']
    c_all = classifications.find({'subject_ids':sid,'user_name':{'$exists':True,'$nin':expert_names}}).sort([("updated_at", -1)])
    clist = list(c_all)
    for c in clist:
        try:
            name = c['user_name']
        except KeyError:
            name = 'Anonymous'
        print '{0:25} {1:25} {2}'.format(name,c['user_ip'],c['updated_at'])

    return None

def rc(zid):

    # Visually compare the expert and volunteer consensus for a subject
    
    plt.ion()

    check_class(zid)
    cons = checksum(zid,excluded=expert_names,no_anonymous=True)
    plot_consensus(cons,figno=1,savefig=False)
    print '\nVolunteers: {0:d} sources'.format(len(cons['answer']))

    cons_ex = checksum(zid,experts_only=True)
    plot_consensus(cons_ex,figno=2,savefig=False)
    print '   Experts: {0:d} sources'.format(len(cons_ex['answer']))

    return None

def run_sample(survey,update=True,subset=None,do_plot=False):

    # Run the consensus algorithm on the ATLAS subjects

    filestem = "consensus_rgz_{0}".format(survey)
    
    if subset is not None:

        '''
        Only run consensus for classifications of 
            expert100: the sample of 100 galaxies classified by science team
            goldstandard: the gold standard sample of 20 galaxies classified by all users

            This only applies to FIRST subjects; no (explicit) gold standard yet for ATLAS,
            although there are the manual classifications in Norris et al. (2006).
        '''

        assert survey == 'first', \
            "Subsets only exist for the FIRST data set, not {0}.".format(survey)

        assert subset in ('expert100','goldstandard'), \
            "Subset is {0}; must be either 'expert100' or 'goldstandard'".format(subset)

        pathd = {'expert100':'expert/expert_all_zooniverse_ids.txt',
                    'goldstandard':'goldstandard/gs_zids.txt'}
        with open('{0}/{1}'.format(rgz_dir,pathd[subset]),'rb') as f:
            zooniverse_ids = [line.rstrip() for line in f]

        suffix = '_{0}'.format(subset)

    else:
        all_completed_zids = [cz['zooniverse_id'] for cz in subjects.find({'state':'complete','metadata.survey':survey})]

        if update:
            '''
            Check to see which subjects have already been completed --
                only run on subjects without an existing consensus.
            '''

            master_json = '{0}/json/{1}.json'.format(rgz_dir,filestem)

            with open(master_json,'r') as fm:
                jmaster = json.load(fm)

            already_finished_zids = []
            for gal in jmaster:
                already_finished_zids.append(gal['zid'])

            zooniverse_ids = list(set(all_completed_zids) - set(already_finished_zids))

            print "\n{0:d} RGZ subjects already in master catalog".format(len(already_finished_zids))
            print "{0:d} RGZ subjects completed since last consensus catalog generation on {1}".format(len(zooniverse_ids),time.ctime(os.path.getmtime(master_json)))

        else:

            # Rerun consensus for every completed subject in RGZ.
            zooniverse_ids = all_completed_zids

        suffix = ''

    # Remove the tutorial subject
    tutorial_zid = "ARG0003r15"
    try:
        zooniverse_ids.remove(tutorial_zid)
    except ValueError:
        print '\nTutorial subject {0} not in list.'.format(tutorial_zid)
    
    print '\nLoaded data; running consensus algorithm on {0:d} completed RGZ subjects'.format(len(zooniverse_ids))

    # Empty files and objects for CSV, JSON output
    json_output = []

    # CSV header
    if update:
        fc = open('{0}/csv/{1}{2}.csv'.format(rgz_dir,filestem,suffix),'a')
    else:
        fc = open('{0}/csv/{1}{2}.csv'.format(rgz_dir,filestem,suffix),'w')
        fc.write('zooniverse_id,{0}_id,n_users,n_total,consensus_level,n_radio,label,bbox,ir_peak\n'.format(survey))

    for idx,zid in enumerate(zooniverse_ids):
    
        # Check progress to screen
        if not idx % 100:
            print idx, datetime.datetime.now().strftime('%H:%M:%S.%f')

        cons = checksum(zid,include_peak_data=do_plot)
        if do_plot:

            plot_consensus(cons,savefig=True)

        # Save results to files

        if cons is not None:

            cons['consensus_level'] = (cons['n_users']/cons['n_total'])

            # JSON

            # Remove peak data from saved catalog; numpy arrays are not JSON serializable (may want to adjust later).
            # http://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113
            for ans in cons['answer']:
                if cons['answer'][ans].has_key('peak_data'):
                    popvar = cons['answer'][ans].pop('peak_data',None)

            json_output.append(cons)

            # CSV

            for ans in cons['answer'].itervalues():
                try:
                    ir_peak = ans['ir_peak']
                except KeyError:
                    ir_peak = ans['ir'] if ans.has_key('ir') else (-99,-99)

                try:
                    fc.write('{0},{1},{2:4d},{3:4d},{4:.3f},{5:2d},{6},"{7}","{8}"\n'.format( 
                            cons['zid'],cons['source'],
                            cons['n_users'],cons['n_total'],cons['consensus_level'],
                            len(ans['xmax']),alphabet(ans['ind']),bbox_unravel(ans['bbox']),ir_peak
                            )
                    )
                except KeyError:
                    print zid
                    print cons

    # Close the new CSV file
    fc.close()

    # Write and close the new JSON file
    if update:
        jmaster.extend(json_output)
        jfinal = jmaster
    else:
        jfinal = json_output

    with open('{0}/json/{1}{2}.json'.format(rgz_dir,filestem,suffix),'w') as fj:
        json.dump(jfinal,fj)

    # Make 75% version for full catalog

    if subset is None:
        # JSON
        json75 = filter(lambda a: (a['n_users']/a['n_total']) >= 0.75, jfinal)
        with open('{0}/json/{1}_75.json'.format(rgz_dir,filestem),'w') as fj:
            json.dump(json75,fj)
        # CSV
        import pandas as pd
        cmaster = pd.read_csv('{0}/csv/{1}.csv'.format(rgz_dir,filestem))
        cmaster75 = cmaster[cmaster['consensus_level'] >= 0.75]
        cmaster75.to_csv('{0}/csv/{1}_75.csv'.format(rgz_dir,filestem),index=False)
        
    print '\nCompleted consensus.'

    return None

def force_csv_update(survey='first',suffix=''):

    # Force an update of the CSV file from the JSON, in case of errors.
    #
    filestem = 'consensus_rgz_{0}'.format(survey)

    master_json = '{0}/json/{1}.json'.format(rgz_dir,filestem)
    
    with open(master_json,'r') as fm:
        jmaster = json.load(fm)
    
    fc = open('{0}/csv/{1}{2}.csv'.format(rgz_dir,filestem,suffix),'w')
    fc.write('zooniverse_id,{0}_id,n_users,n_total,consensus_level,n_radio,label,bbox,ir_peak\n'.format(survey))

    for gal in jmaster:
        for ans in gal['answer'].itervalues():
            try:
                ir_peak = ans['ir_peak']
            except KeyError:
                ir_peak = ans['ir'] if ans.has_key('ir') else (-99,-99)

    
            fc.write('{0},{1},{2:4d},{3:4d},{4:.3f},{5:2d},{6},"{7}","{8}"\n'.format(
                    gal['zid'],
                    gal['source'],
                    gal['n_users'],
                    gal['n_total'],
                    gal['n_users'] * 1./gal['n_total'],
                    len(ans['xmax']),
                    alphabet(ans['ind']),
                    bbox_unravel(ans['bbox']),
                    ir_peak
                    )
                )

    fc.close()

    return None

def bbox_unravel(bbox):

    # Turn an array of tuple strings into floats

    bboxes = []
    for lobe in bbox:
        t = [float(x) for x in lobe]
        t = tuple(t)
        bboxes.append(t)

    return bboxes

def alphabet(i):

    from string import letters

    # Return a letter (or set of duplicated letters) of the alphabet for a given integer
    # For i > 25, letters begin multiplying:    alphabet(0) = 'a' 
    #                                           alphabet(1) = 'b'
    #                                           ...
    #                                           alphabet(25) = 'z'
    #                                           alphabet(26) = 'aa'
    #                                           alphabet(27) = 'bb'
    #                                           ...
    #
    lowercase = letters[26:]
    
    try:
        letter = letters[26:][i % 26]*int(i/26 + 1)
        return letter
    except TypeError:
        raise AssertionError("Index must be an integer")

def update_experts(classifications): 

    # Add field to classifications made by members of the expert science team. Takes ~1 minute to run.

    import dateutil.parser

    # Load saved data from the test runs
    json_data = open('{0}/expert/expert_params.json'.format(rgz_dir)).read() 
    experts = json.loads(json_data)

    for ex in experts:
        expert_dates = (dateutil.parser.parse(ex['started_at']),dateutil.parser.parse(ex['ended_at']))
        classifications.update({"updated_at": {"$gt": expert_dates[0],"$lt":expert_dates[1]},"user_name":ex['expert_user']},{'$set':{'expert':True}},multi=True)

    return None

def update_gs_subjects(subjects): 

    # Add field to the Mongo database designating the gold standard subjects.

    with open('{0}/goldstandard/gs_zids.txt'.format(rgz_dir),'r') as f:
        for gal in f:
            subjects.update({'zooniverse_id':gal.strip()},{'$set':{'goldstandard':True}})

    return None

if __name__ == "__main__":

    if pathdict != None:

        print 'Starting at',datetime.datetime.now().strftime('%H:%M:%S.%f')

        for survey in ('atlas','first'):
            run_sample(survey,update=True,do_plot=False)

        print 'Finished at',datetime.datetime.now().strftime('%H:%M:%S.%f')
    else:
        print "\nAborting consensus.py - could not locate raw RGZ image data.\n"
