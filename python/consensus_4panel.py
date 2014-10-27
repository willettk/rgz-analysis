from pymongo import MongoClient
from collections import Counter
import datetime
import numpy as np
import operator
import pprint
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from scipy.ndimage.filters import maximum_filter
from scipy import stats
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.linalg.basic import LinAlgError
import matplotlib.patches as patches
from matplotlib.path import Path
import requests
from StringIO import StringIO
from PIL import Image
import os
import cStringIO
import urllib

# MongoDB parameters

client = MongoClient('localhost', 27017)
db = client['radio'] 

subjects = db['radio_subjects'] 		# subjects = images
classifications = db['radio_classifications']	# classifications = classifications of each subject per user

# General variables for the RGZ sample

main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)
IMG_HEIGHT = 424.0			# number of pixels in the JPG image along the y axis
IMG_WIDTH = 424.0			# number of pixels in the JPG image along the x axis
IMG_HEIGHT = 500.0			# number of pixels in the JPG image along the y axis
IMG_WIDTH = 500.0			# number of pixels in the JPG image along the x axis
FITS_HEIGHT = 301.0			# number of pixels in the FITS image along the y axis
FITS_WIDTH = 301.0			# number of pixels in the FITS image along the x axis
PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image
xmin = 1.
xmax = IMG_HEIGHT
ymin = 1.
ymax = IMG_WIDTH
bad_keys = ('finished_at','started_at','user_agent','lang')

# Paths

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

# Find the consensus classification for a single subject

def checksum(z,save_fig=None):

    sub = subjects.find_one({'zooniverse_id':z})
    imgid = sub['_id']
    zid = sub['zooniverse_id']

    # Classifications for this subject after launch date
    #clist = list(classifications.find({"subject_ids": imgid, "updated_at": {"$gt": main_release_date},"expert":{"$exists":False}}))
    clist = list(classifications.find({"subject_ids": imgid, "updated_at": {"$gt": main_release_date},"expert":True}))
    
    # Empty dicts and lists 
    cdict = {}
    checksum_list = []
    
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
    
        checksum = sum(sumlist)
        checksum_list.append(checksum)
        c['checksum'] = checksum
    
        # Insert checksum into dictionary with number of galaxies as the index
        if cdict.has_key(n_galaxies):
            cdict[n_galaxies].append(checksum)
        else:
            cdict[n_galaxies] = [checksum]
    
    #print cdict,'\n'
    
    maxval=0
    mc_checksum = 0.
    ngals = 0

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
            ngals = k
    
    # Find a galaxy that matches the checksum (easier to keep track as a list)
    
    cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
   
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in bad_keys]

    # Find the sum of the xmax coordinates for each galaxy. This gives the index to search on.
    
    consensus = {}
    ir_x,ir_y = {},{}
    for k,gal in enumerate(goodann):
        xmax_temp = []
        for v in gal['radio'].itervalues():
            xmax_temp.append(float(v['xmax']))
        checksum2 = round(sum(xmax_temp),3)
        consensus[checksum2] = {}
        consensus[checksum2]['ind'] = k
        consensus[checksum2]['xmax'] = xmax_temp
    
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
                    k = consensus[xmax_checksum]['ind']
    
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
        x_exists = [xt * scale_ir for xt in xv if xt != -99.0]
        y_exists = [yt * scale_ir for yt in yv if yt != -99.0]
    
        pd = {'ind':xk}
    
        if len(Counter(x_exists)) > 2:
            X, Y = np.mgrid[xmin:xmax, ymin:ymax]
            positions = np.vstack([X.ravel(), Y.ravel()])
            values = np.vstack([x_exists, y_exists])
            try:
                kernel = stats.gaussian_kde(values)
            except LinAlgError:
                print 'LinAlgError in KD estimation for %s' % zid,x_exists,y_exists
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
    
            #return X,Y,Z,npeaks
    
            pd['X'] = X
            pd['Y'] = Y
            pd['Z'] = Z
            pd['npeaks'] = npeaks
            peak_data.append(pd)
        
        else:
    
            pd['npeaks'] = 0
            peak_data.append(pd)
    
    # Plot image
    
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
        for v in consensus.itervalues():
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
        print 'Users found no components for consensus match of %s' % z
    
    # Plot the infrared results
    
    fig = plt.figure(1,(15,4))
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)
    
    colormaparr = [cm.hot_r,cm.Blues,cm.RdPu,cm.Greens,cm.PuBu,cm.YlGn,cm.Greys][::-1]
    colorarr = ['r','b','m','g','c','y','k'][::-1]
    
    if ngals > 0: # At least one galaxy was identified
        for idx,(xgal,ygal,peak) in enumerate(zip(ir_x.itervalues(),ir_y.itervalues(),peak_data)):

            if peak['npeaks'] > 0:
        
                # Find the peak
                xpeak = float(peak['X'][peak['Z']==peak['Z'].max()][0])
                ypeak = float(peak['Y'][peak['Z']==peak['Z'].max()][0])
        
                # Plot the KDE map
                colormap = colormaparr.pop()
                ax3.imshow(np.rot90(peak['Z']), cmap=colormap,extent=[xmin, xmax, ymin, ymax])
        
                # Plot individual sources
                if peak['npeaks'] > 0:
                    color = colorarr.pop()
                    x_exists = [xt * 500./424 for xt in xgal if xt != -99.0]
                    y_exists = [yt * 500./424 for yt in ygal if yt != -99.0]
                    ax3.scatter(x_exists, y_exists, c=color, marker='o', s=8, alpha=1./ngals)
                    ax4.plot([xpeak],[ypeak],color=color,marker='*',markersize=12)
    
                else:
                    ax4.text(550,idx*25,'No IR host for galaxy #%i' % idx,fontsize=11)

            else:
                xpeak,ypeak=(xgal[0],ygal[0])
    
            # Append the IR location to the final consensus
            for k,v in consensus.iteritems():
                if v['ind'] == peak['ind']:
                    consensus[k]['ir'] = (xpeak,ypeak)
        
    # Scaling factor for FITS to radio files
    
    radio_ir_scaling_factor = 500./132
    
    ax3.set_xlim([0, 500])
    ax3.set_ylim([500, 0])
    ax3.set_title(zid)
    ax3.set_aspect('equal')
    
    ax4.set_xlim([0, 500])
    ax4.set_ylim([500, 0])
    ax4.set_title('Consensus (%i/%i users)' % (maxval,len(clist)))
    
    ax4.set_aspect('equal')
    
    # Display IR and radio images
    
    url_standard = sub['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax1.set_title('WISE')

    url_radio = sub['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    ax2.set_title(sub['metadata']['source'])
    ax2.get_yaxis().set_ticklabels([])

    ax3.get_yaxis().set_ticklabels([])

    # Plot contours identified as the consensus
    if len(consensus) > 0:
        ax4.add_patch(patch_black)
    ax4.yaxis.tick_right()
    
    ax1.get_xaxis().set_ticks([0,100,200,300,400])
    ax2.get_xaxis().set_ticks([0,100,200,300,400])
    ax3.get_xaxis().set_ticks([0,100,200,300,400])
    ax4.get_xaxis().set_ticks([0,100,200,300,400,500])

    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    if save_fig:
        #writefile = '/Volumes/3TB/rgz/plots/expert_%s.pdf' % zid
        writefile = '%s/plots/expert/expert_all/expert_%s.pdf' % (rgz_dir,zid)
        fig.savefig('%s' % (writefile))
    else:
        plt.show()

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    plt.close()

    return consensus
    
def run_sample():

    # Run all galaxies in the expert sample
    
    N = 0
    
    #with open('%s/goldstandard/gs_zids.txt' % rgz_dir,'rb') as f:
    with open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir,'rb') as f:
        zooniverse_ids = [line.rstrip() for line in f]
    
    for z in zooniverse_ids:
    
        N += 1
        consensus = checksum(z,save_fig=True)

        # Check progress by printing to screen every 100 classifications
        if not N % 100:
            print N, datetime.datetime.now().strftime('%H:%M:%S.%f')


