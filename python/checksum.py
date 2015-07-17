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
import matplotlib.patches as patches
from matplotlib.path import Path
import requests
from StringIO import StringIO
from PIL import Image
import os
import cStringIO
import urllib

client = MongoClient('localhost', 27017)
db = client['radio'] 

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

subjects = db['radio_subjects'] 		# subjects = images
classifications = db['radio_classifications']	# classifications = classifications of each subject per user

# Running as part of loop

N = 0
n_subjects = 10
class_lim = {'state':'complete','classification_count':{'$gt':19}}
class_lim = {'zooniverse_id':'ARG0000jqj'}

# Look at just the newly retired ones (single-contour, 5 classifications)
#class_lim = {'state':'complete','metadata.contour_count':1,'classification_count':5}

for sub in list(subjects.find(class_lim).limit(n_subjects)):

    N += 1
    
    #sub = subjects.find_one({'state':'complete','metadata.contour_count':3,'classification_count':20})
    #sub = subjects.find_one({'zooniverse_id':'ARG0000jqj'})
    imgid = sub['_id']
    zid = sub['zooniverse_id']
    clist = list(classifications.find({"subject_ids": imgid, "updated_at": {"$gt": main_release_date}}))
    
    bad_keys = ('finished_at','started_at','user_agent','lang')
    
    cdict = {}
    checksum_list = []
    
    for c in clist:
        # Want most popular combo for each NUMBER of galaxies
        
        sumlist = []    # List of the checksums over all possible combinations
        goodann = [x for x in c['annotations'] if x.keys()[0] not in bad_keys]
        n_galaxies = len(goodann)
    
        for idx,ann in enumerate(goodann):
    
            if ann.has_key('started_at') or ann.has_key('finished_at') or ann.has_key('user_agent') or ann.has_key('lang'):
                continue
    
            xmaxlist = []
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
    
        checksum = sum(sumlist)
        checksum_list.append(checksum)
        c['checksum'] = checksum
    
        if cdict.has_key(n_galaxies):
            cdict[n_galaxies].append(checksum)
        else:
            cdict[n_galaxies] = [checksum]
    
    #print cdict,'\n'
    
    maxval=0
    mc_checksum = 0.
    ngals = 0
    for k,v in cdict.iteritems():
        mc = Counter(v).most_common()
        #print '%i galaxies: %s' % (k,mc)
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
    
    suffix = 'y' if ngals == 1 else 'ies'
    #print 'Most common selection (%i/%i users) is %i galax%s\n' % (maxval,len(clist),ngals,suffix)
    
    # Find a galaxy that matches the checksum (easier to keep track as a list)
    
    cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
    #pprint.pprint(cmatch['annotations'])
    #print ''
    #print 'http://radiotalk.galaxyzoo.org/#/subjects/%s\n' % zid
    
    #print 'imgid: %s\n'%imgid
    
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in bad_keys]
    
    # I know how many galaxies there should be. Need to make sure to match the correct IR hosts with radio components.
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
    
    scale_ir = 500./424.
    
    peak_data = []
    for x,y in zip(ir_x.itervalues(),ir_y.itervalues()):
        x_exists = [xt * scale_ir for xt in x if xt != -99.0]
        y_exists = [yt * scale_ir for yt in y if yt != -99.0]
    
        pd = {}
    
        if len(x_exists) > 2:
            X, Y = np.mgrid[xmin:xmax, ymin:ymax]
            positions = np.vstack([X.ravel(), Y.ravel()])
            values = np.vstack([x_exists, y_exists])
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
    
    writefile = '/Volumes/3TB/rgz/plots/%s.pdf' % zid
        
    # Download contour data
    
    r = requests.get(sub['location']['contours'])
    contours = r.json()
    
    sf_x = 500./contours['width']
    sf_y = 500./contours['height']
    
    verts_all = []
    codes_all = []
    components = contours['contours']
    
    bboxes_xmax = []
    
    for comp in components:
    
        # Order of bounding box components is (xmax,ymax,xmin,ymin)
        bboxes_xmax.append(comp[0]['bbox'][0])
    
        for idx,level in enumerate(comp):
            verts = [((p['x'])*sf_x,(p['y']-1)*sf_y) for p in level['arr']]
            
            codes = np.ones(len(verts),int) * Path.LINETO
            codes[0] = Path.MOVETO
        
            verts_all.extend(verts)
            codes_all.extend(codes)
    
    path = Path(verts_all, codes_all)
    patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)
    
    # Plot the infrared results
    
    fig = plt.figure(1,(15,4))
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)
    
    colormaparr = [cm.hot_r,cm.Blues,cm.RdPu,cm.Greens,cm.PuBu,cm.YlGn,cm.Greys][::-1]
    colorarr = ['r','b','m','g','c','y','k'][::-1]
    
    if ngals > 0: # More than one galaxy was identified
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
                    #ax3.text(50,40+idx*35,r'IR peak: $(%i,%i)$' % (xpeak,ypeak),color='k',fontsize=12)
                    x_exists = [xt * 500./424 for xt in xgal if xt != -99.0]
                    y_exists = [yt * 500./424 for yt in ygal if yt != -99.0]
                    ax3.scatter(x_exists, y_exists, c=color, marker='o', s=8, alpha=1./ngals)
                    ax4.plot([xpeak],[ypeak],color=color,marker='*',markersize=12)
    
                # Plot contours here?
    
                # Would need to individually again find the xmax values that make up the checksum of the radio data
                # That's stored in gal_xmax, along with individual values
    
    
                else:
                    ax4.text(550,idx*25,'No IR host for galaxy #%i' % idx,fontsize=11)
    
        
    
    '''
    # Plot the radio counts
    
    radio_flattened = [item for sublist in all_radio for item in sublist]
    uniques = set(radio_flattened)
    d = dict(zip(uniques,np.arange(len(uniques))))
    c = Counter(all_radio)
    cmc = c.most_common()[::-1]
    
    # Sort by number of components
    
    for idx,(c_xval,n) in enumerate(cmc):
        if len(c_xval) > 1:
            tlist = [str(d[x]) for x in c_xval]
            t = ' and R'.join(sorted(tlist))
        else:
            t = d[c_xval[0]]
        singular = 's' if n != 1 else ''
        ax3.text(550,400-idx*25,'%3i vote%s: R%s' % (n,singular,t),fontsize=11)
    
    
    '''
    # Scaling factor for FITS to radio files
    
    radio_ir_scaling_factor = 500./132
    
    '''
    # Rectangle showing the radio box size
    
    box_counts = Counter(radio_flattened)
    for ru in radio_unique:
    
        x0,x1,y0,y1 = [float(ru_) * radio_ir_scaling_factor for ru_ in ru]
    
        # Assume xmax matching is still good
    
        xmax_index = '%.6f' % float(ru[1])
        component_number = d[xmax_index]
        number_votes = box_counts[xmax_index]
        rectangle = plt.Rectangle((x0,y0), x1-x0, y1-y0, fill=False, linewidth=number_votes/5., edgecolor = 'c')
        ax3.add_patch(rectangle)
        ax3.text(x0-15,y0-15,'R%s' % component_number)
    
    '''
    ax3.set_xlim([0, 500])
    ax3.set_ylim([500, 0])
    #ax3.set_title(zid)
    ax3.set_aspect('equal')
    
    ax4.set_xlim([0, 500])
    ax4.set_ylim([500, 0])
    #ax4.set_title('Consensus (%i/%i users)' % (maxval,len(clist)))
    
    ax4.set_aspect('equal')
    
    # Display IR and radio images
    
    url_standard = sub['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax4.add_patch(patch_black)
    #ax1.set_title('WISE')
    
    url_radio = sub['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    #ax2.set_title(sub['metadata']['source'])

    plt.subplots_adjust(wspace=0.02)

    ax2.get_yaxis().set_ticklabels([])
    ax3.get_yaxis().set_ticklabels([])
    ax4.yaxis.tick_right()

    ax1.get_xaxis().set_ticks([0,100,200,300,400])
    ax2.get_xaxis().set_ticks([0,100,200,300,400])
    ax3.get_xaxis().set_ticks([0,100,200,300,400])
    ax4.get_xaxis().set_ticks([0,100,200,300,400,500])

    
    #fig.show()
    
    # Save hard copy of the figure
    fig.savefig('%s' % (writefile))
    fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/paper/reduced_example.eps')
    
    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    plt.close()
    
    os.system('open ' + writefile)
    
    # Check progress by printing to screen every 100 classifications
    if not N % 100:
        print N, datetime.datetime.now().strftime('%H:%M:%S.%f')

