from pymongo import MongoClient
from collections import Counter
import datetime
import numpy as np
import operator

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

from scipy.ndimage.filters import maximum_filter
from scipy import stats
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.linalg.basic import LinAlgError

from astropy.io import fits

import requests
from PIL import Image
import cStringIO
import urllib
from subprocess import call

import collinearity
import consensus
import rgz

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
PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image
xmin = 1.
xmax = IMG_HEIGHT_NEW
ymin = 1.
ymax = IMG_WIDTH_NEW

bad_keys = ('finished_at','started_at','user_agent','lang','pending')

expert_names = [u'42jkb', u'ivywong', u'stasmanian', u'klmasters', u'Kevin', u'akapinska', u'enno.middelberg', u'xDocR', u'vrooje', u'KWillett', u'DocR']

subjects,classifications = rgz.load_rgz_data()

# Paths

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

# Plot 2 volunteer + 1 expert consensus for one galaxy, showing how they differ

def plot_pce(zid = 'ARG000180p',savefig = False):

    cons1 = checksum(zid,no_anonymous=True,rank=1)
    cons2 = checksum(zid,no_anonymous=True,rank=2)
    cons3 = checksum(zid,experts_only=True,rank=1)
    
    #consensus.plot_consensus(cons1,figno=4, savefig=None)
    #consensus.plot_consensus(cons2,figno=5, savefig=None)
    #consensus.plot_consensus(cons3,figno=6, savefig=None)

    # Plot the infrared results
    
    fig = plt.figure(3,(12,12))
    fig.clf()

    gs1 = gridspec.GridSpec(2, 1)
    gs1.update(left=0.05, right=0.48, wspace=0.05,top=0.95,bottom=0.05)
    ax_wise  = fig.add_subplot(gs1[0,:])
    ax_first = fig.add_subplot(gs1[1,:])
    
    gs2 = gridspec.GridSpec(3, 2)
    gs2.update(left=0.45, right=0.98, hspace=0.05,top=0.95,bottom=0.05)
    axarr1 = [fig.add_subplot(gs2[0,0]),fig.add_subplot(gs2[0,1])]
    axarr2 = [fig.add_subplot(gs2[1,0]),fig.add_subplot(gs2[1,1])]
    axarr3 = [fig.add_subplot(gs2[2,0]),fig.add_subplot(gs2[2,1])]

    #fig,axarr = plt.subplots(3,4,num=3,figsize=(15,12),sharex='col', sharey='row')
    #fig.clf()
    title1 = '%3i/%3i volunteers (#1)' % (cons1['n_users'],cons1['n_total'])
    title2 = '%3i/%3i volunteers (#2)' % (cons2['n_users'],cons2['n_total'])
    title3 = '%i/%i experts'    % (cons3['n_users'],cons3['n_total'])

    plot_images_only(cons1,ax_wise,ax_first)
    plot_consensus_only(cons1, axarr1, title1, set_title=True,set_xticks=False, set_yticks=True)
    plot_consensus_only(cons2, axarr2, title2, set_xticks=False, set_yticks=True)
    plot_consensus_only(cons3, axarr3, title3, set_xticks=True, set_yticks=True)
    
    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    # For some reason, LaTeX breaks when saving this as EPS. Save as PNG and then use ImageMagick to convert.
    if savefig:
        writepath = '/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/paper/figures'
        writefile = 'compare_consensus'
        fig.savefig('%s/%s.png' % (writepath,writefile))
        call(["convert",'%s/%s.png' % (writepath,writefile), '%s/%s.eps' % (writepath,writefile))
        plt.close()
    else:
        plt.show()

    return None

def best_checksum(cdict):

    maxval=0
    mc_checksum = 0.
    best_k = 0

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
            best_k = k

    return mc_checksum,best_k,maxval

def checksum(zid,experts_only=False,excluded=[],no_anonymous=False,rank=1):

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

    clist_all = list(classifications.find(class_params))

    # Empty dicts and lists 
    cdict = {}
    checksum_list = []

    unique_users = set()
    
    clen_start = len(clist_all)
    listcount = []

    # Compute the most popular combination for each NUMBER of galaxies identified in image
    for c in clist_all:

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
    
            checksum = round(sum(sumlist),3)
            checksum_list.append(checksum)
            c['checksum'] = checksum
    
            # Insert checksum into dictionary with number of galaxies as the index
            if cdict.has_key(n_galaxies):
                cdict[n_galaxies].append(checksum)
            else:
                cdict[n_galaxies] = [checksum]

        else:
            listcount.append(False)
            #print 'Removing classification for %s' % user_name
    
    # Remove duplicates and classifications for no object
    clist = [c for lc,c in zip(listcount,clist_all) if lc and c['checksum'] != -99]

    clen_diff = clen_start - len(clist)

    '''
    if clen_diff > 0:
        print '\nSkipping %i duplicated classifications for %s. %i good classifications total.' % (clen_diff,zid,len(clist))
    '''

    # SECOND/THIRD BEST - drop the best checksum, then run through again?

    mc_checksum,best_k,maxval = best_checksum(cdict)
    for i in np.arange(rank-1):
        cbest = cdict[best_k]
        cdict[best_k] = filter(lambda x: x != mc_checksum, cbest)
        mc_checksum,best_k,maxval = best_checksum(cdict)

    # Find a galaxy that matches the checksum (easier to keep track as a list)

    try:
        cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
    except StopIteration:
        # Necessary for objects like ARG0003par; one classifier recorded 22 "No IR","No Contours" in a short space. Still shouldn't happen.
        print 'No non-zero classifications recorded for %s' % zid
        return None
   
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
        try:
            for v in gal['radio'].itervalues():
                xmax_temp.append(float(v['xmax']))
            checksum2 = round(sum(xmax_temp),3)
            answer[checksum2] = {}
            answer[checksum2]['ind'] = k
            answer[checksum2]['xmax'] = xmax_temp
        except KeyError:
            print gal, zid
        except AttributeError:
            print 'No Sources, No IR recorded for %s' % zid
    
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


    # Perform a kernel density estimate on the data for each galaxy
    
    scale_ir = 500./424.

    peak_data = []
    for (xk,xv),(yk,yv) in zip(ir_x.iteritems(),ir_y.iteritems()):

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
                print zid
                print x_exists,y_exists
                print len(x_exists),len(y_exists)
            try:
                kernel = stats.gaussian_kde(values)
            except LinAlgError:
                print 'LinAlgError in KD estimation for %s' % zid,x_exists,y_exists
                continue

            # Even if there are more than 2 sets of points, if they are mutually co-linear, 
            # matrix can't invert and kernel returns NaNs. 

            if np.isnan(kernel(positions)).sum() > 0:
                acp = collinearity.collinear(x_exists,y_exists)
                if len(acp) > 0:
                    print 'There are %i unique points for %s (source no. %i in the field), but all are co-linear; KDE estimate does not work.' % (len(Counter(x_exists)),zid,xk)
                else:
                    print 'There are NaNs in the KDE for %s (source no. %i in the field), but points are not co-linear.' % (zid,xk)

                for k,v in answer.iteritems():
                    if v['ind'] == xk:
                        answer[k]['ir'] = (np.mean(x_exists),np.mean(y_exists))
        
            else:

                Z = np.reshape(kernel(positions).T, X.shape)
                
                # Find the number of peaks
                # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array
                
                neighborhood = np.ones((10,10))
                local_max = maximum_filter(Z, footprint=neighborhood)==Z
                background = (Z==0)
                eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
                detected_peaks = local_max - eroded_background
                
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
                        answer[k]['peak_data'] = pd
                        answer[k]['ir_x'] = x_exists
                        answer[k]['ir_y'] = y_exists
        
        else:

            for k,v in answer.iteritems():
                if v['ind'] == xk:
                    # Case 1: multiple users selected IR source, but not enough unique points to pinpoint peak
                    if most_common_ir != (-99,-99) and len(x_exists) > 0 and len(y_exists) > 0:
                        answer[k]['ir'] = (x_exists[0],y_exists[0])
                    # Case 2: most users have selected No Sources
                    else:
                        answer[k]['ir'] = (-99,-99)

    return cons

def plot_images_only(consensus,ax_wise,ax_first):

    # Plot image
    
    zid = consensus['zid']
    answer = consensus['answer']
    sub = subjects.find_one({'zooniverse_id':zid})

    # Display IR and radio images
    
    url_standard = sub['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax_wise.imshow(im_standard,origin='upper')
    ax_wise.set_title('WISE')

    url_radio = sub['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax_first.imshow(im_radio,origin='upper')
    #ax_first.set_title(sub['metadata']['source'])
    ax_first.set_title('FIRST')
    #ax_first.get_yaxis().set_ticklabels([])

    #ax_wise.get_xaxis().set_ticks([0,100,200,300,400])
    #ax_first.get_xaxis().set_ticks([0,100,200,300,400])

    return None

def plot_consensus_only(consensus,axes,title,set_title=False,set_xticks=False,set_yticks=False):

    axis2,axis3 = axes

    # Plot image
    
    zid = consensus['zid']
    answer = consensus['answer']
    sub = subjects.find_one({'zooniverse_id':zid})

    # Download contour data
    
    r = requests.get(sub['location']['contours'])
    contours = r.json()
    
    sf_x = 500./contours['width']
    sf_y = 500./contours['height']
    
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
        print 'Users found no components for consensus match of %s' % zid
    
    colormaparr = [cm.hot_r,cm.Blues,cm.RdPu,cm.Greens,cm.PuBu,cm.YlGn,cm.Greys][::-1]
    colorarr = ['r','b','m','g','c','y','k'][::-1]
    
    if len(answer) > 0: # At least one galaxy was identified
        for idx,ans in enumerate(answer.itervalues()):

            if ans.has_key('peak_data'):

                # Plot the KDE map
                colormap = colormaparr.pop()
                axis2.imshow(np.rot90(ans['peak_data']['Z']), cmap=colormap,extent=[xmin, xmax, ymin, ymax])
        
                # Plot individual sources
                color = colorarr.pop()
                '''
                x_plot = [xt * 500./424 for xt in ans['ir_x'] if xt != -99.0]
                y_plot = [yt * 500./424 for yt in ans['ir_y'] if yt != -99.0]
                '''
                x_plot,y_plot = ans['ir_x'],ans['ir_y']
                axis2.scatter(x_plot, y_plot, c=color, marker='o', s=10, alpha=1./len(x_plot))
                axis3.plot([ans['ir_peak'][0]],[ans['ir_peak'][1]],color=color,marker='*',markersize=12)
    
            elif ans.has_key('ir'):
                color = colorarr.pop()
                x_plot,y_plot = ans['ir']
                axis2.plot([x_plot],[y_plot],color=color,marker='o',markersize=2)
                axis3.plot([x_plot],[y_plot],color=color,marker='*',markersize=12)
            else:
                axis3.text(550,idx*25,'#%i - no IR host' % idx,fontsize=11)

    
    axis2.set_xlim([0, 500])
    axis2.set_ylim([500, 0])
    if set_title:
        axis2.set_title('KDE')
        axis3.set_title('Consensus')
    axis2.set_aspect('equal')
    
    axis3.text(25,50,title,horizontalalignment='left')

    axis3.set_xlim([0, 500])
    axis3.set_ylim([500, 0])
    
    axis3.set_aspect('equal')
    
    if not set_xticks:
        axis2.get_xaxis().set_ticklabels([])
        axis3.get_xaxis().set_ticklabels([])

    axis2.get_yaxis().set_ticklabels([])

    # Plot contours identified as the consensus
    if len(answer) > 0:
        axis3.add_patch(patch_black)
    axis3.yaxis.tick_right()
    
    #axis2.get_xaxis().set_ticks([0,100,200,300,400])
    #axis3.get_xaxis().set_ticks([0,100,200,300,400,500])

    return None
