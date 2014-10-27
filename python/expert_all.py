# import necessary python packages
import numpy as np
import datetime
import os,sys
import urllib
import json

from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import dateutil.parser
from matplotlib import pyplot as plt
from matplotlib import cm
from collections import Counter
from pymongo import MongoClient
from scipy import stats
from PIL import Image
from collections import OrderedDict
from scipy.linalg.basic import LinAlgError

import matplotlib.patches as patches
from matplotlib.path import Path
import requests
from StringIO import StringIO
#------------------------------------------------------------------------------------------------------------

# Setup path locations

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
csv_dir = '%s/csv' % rgz_dir

plot_dir = '%s/plots/expert' % rgz_dir
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

dat_dir = '%s/datfiles/expert/expert_all' % rgz_dir
if not os.path.isdir(dat_dir):
    os.mkdir(dat_dir)

# Set constants
beta_release_date = datetime.datetime(2013, 10, 20, 12, 0, 0, 0)	# date of beta release (YYY,MM,DD,HH,MM,SS,MS)
main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)

IMG_HEIGHT = 424.0			# number of pixels in the JPG image along the y axis
IMG_WIDTH = 424.0			# number of pixels in the JPG image along the x axis
FITS_HEIGHT = 301.0			# number of pixels in the FITS image along the y axis
FITS_WIDTH = 301.0			# number of pixels in the FITS image along the x axis
PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image

xmin = 1.
xmax = IMG_HEIGHT
ymin = 1.
ymax = IMG_WIDTH

xjpg2fits = float(IMG_WIDTH/FITS_WIDTH)		# map the JPG pixels to the FITS pixels in x
yjpg2fits = float(IMG_HEIGHT/FITS_HEIGHT)	# map the JPG pixels to the FITS pixels in y

bad_keys = ('finished_at','started_at','user_agent','lang')

plt.ion()

def list_flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(list_flatten(el))
        else:
            result.append(el)

    return result

def plot_npeaks(volunteers=False):

    if volunteers:
        readfile = 'expert/npeaks_ir_expert_volunteer.csv'
        writefile = 'expert_all/ir_peaks_histogram_expert_all.png'
    else:
        readfile = 'expert/npeaks_ir_expert_all.csv'
        writefile = 'expert_all/ir_peaks_histogram_expert_volunteer.png'
        
    # Read in data

    with open('%s/%s' % (csv_dir,readfile),'rb') as f:
        npeaks = [int(line.rstrip()) for line in f]
    
    # Plot the distribution of the total number of IR sources per image

    fig = plt.figure(figsize=(8,7))
    ax1 = fig.add_subplot(111)

    h = plt.hist(npeaks,bins=np.arange(np.max(npeaks)+1),axes=ax1)

    ax1.set_title('RGZ source distribution')
    ax1.set_xlabel('Number of IR peaks per image')
    ax1.set_ylabel('Count')

    #fig.show()
    fig.tight_layout()

    # Save hard copy of the figure
    fig.savefig('%s/%s' % (plot_dir,writefile))

    plt.close()

    return None

def find_ir_peak(x,y):

    # Perform a kernel density estimate on the data:
    
    X, Y = np.mgrid[xmin:xmax, ymin:ymax]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([x, y])
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
    
    return X,Y,Z,npeaks

def plot_image(ir_x,ir_y,sub,X,Y,Z,npeaks,all_radio,radio_unique,volunteers=False):

    if volunteers:
        writefile = 'expert_volunteers/ir_peaks/%s_ir_peak.png' % sub['zooniverse_id']
    else:
        writefile = 'expert_all/ir_peaks/%s_ir_peak.png' % sub['zooniverse_id']
        
    # Plot the infrared results
    
    fig = plt.figure(1,(15,4))
    ax3 = fig.add_subplot(143)

    assert len(ir_x) == len(ir_y), \
        'Length of IR x- and y-vectors must be the same' 

    if len(ir_x) > 2:

        # Find the peak

        xpeak = X[Z==Z.max()][0]
        ypeak = Y[Z==Z.max()][0]

        # Plot the KDE map
        ax3.imshow(np.rot90(Z), cmap=plt.cm.hot_r,extent=[xmin, xmax, ymin, ymax])

    # Plot the individual sources

    if len(ir_x) > 2:       # At least 3 non-singular IR points means that the KDE map can be generated and peaks estimated
        ax3.text(50,40,r'Main IR peak: $(%i,%i)$' % (xpeak,ypeak),color='k',fontsize=12)
        ax3.text(50,70,r'$N_{IR\/peaks}$ = %i' % npeaks,color='k',fontsize=12)
        ax3.plot(ir_x, ir_y, 'go', markersize=4)
        ax3.plot([xpeak],[ypeak],'c*',markersize=12)
    elif len(ir_x) == 2:    # 2 IR points are simply plotted, with no KDE estimation
        ax3.text(50,40,r'IR peaks: $(%i,%i),(%i,%i)$' % (ir_x[0],ir_y[0],ir_x[1],ir_y[1]),color='k',fontsize=12)
        ax3.text(50,70,r'$N_{IR\/peaks}$ = 2',color='k',fontsize=12)
        ax3.plot(ir_x,ir_y,'c*',markersize=12)
    elif len(ir_x) == 1:    # 1 IR point is simply plotted, with no KDE estimation
        print ir_x,ir_y
        ax3.text(50,40,r'IR peak: $(%i,%i)$' % (ir_x[0],ir_y[0]),color='k',fontsize=12)
        ax3.text(50,70,r'$N_{IR\/peaks}$ = 1',color='k',fontsize=12)
        ax3.plot(ir_x,ir_y,'c*',markersize=12)
    if len(ir_x) == 0:      # 0 IR points identified by any user
        ax3.text(50,70,'No IR sources',color='k',fontsize=12)

    # Plot the radio counts

    radio_flattened = [item for sublist in all_radio for item in sublist]
    uniques = set(radio_flattened)
    d = dict(zip(uniques,np.arange(len(uniques))))
    c = Counter(all_radio)
    cmc = c.most_common()[::-1]

    # Sort by number of components?

    for idx,(c_xval,n) in enumerate(cmc):
        if len(c_xval) > 1:
            tlist = [str(d[x]) for x in c_xval]
            t = ' and R'.join(sorted(tlist))
        else:
            t = d[c_xval[0]]
        singular = 's' if n != 1 else ''
        ax3.text(550,400-idx*25,'%3i vote%s: R%s' % (n,singular,t),fontsize=11)

    # Download contour data

    r = requests.get(sub['location']['contours'])
    contours = r.json()
    
    sf_x = 500./contours['width']
    sf_y = 500./contours['height']

    verts_all = []
    codes_all = []
    components = contours['contours']
    
    for comp in components:
    
        for idx,level in enumerate(comp):
            verts = [((p['x'])*sf_x,(p['y']-1)*sf_y) for p in level['arr']]
            
            codes = np.ones(len(verts),int) * Path.LINETO
            codes[0] = Path.MOVETO
        
            verts_all.extend(verts)
            codes_all.extend(codes)
    
    path = Path(verts_all, codes_all)
    patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)

    # Scaling factor for FITS to radio files

    radio_ir_scaling_factor = 500./132

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

    ax3.set_xlim([0, 500])
    ax3.set_ylim([500, 0])
    ax3.set_title('%s \n %s' % (sub['zooniverse_id'],sub['metadata']['source']))
    ax3.set_aspect('equal')

    # Display IR and radio images

    url_standard = sub['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax1.add_patch(patch_black)
    ax1.set_title('WISE')

    url_radio = sub['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    ax2.set_title('FIRST')

    if volunteers:
        ax2.set_xlabel('volunteers')
    else:
        ax2.set_xlabel('experts')

    # Save hard copy of the figure
    fig.savefig('%s/%s' % (plot_dir,writefile))

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    plt.close()

    return None
    
def find_consensus(sub,classifications,verbose=False,volunteers=False):
    
    Nclass = sub["classification_count"]	# number of classifications made per image
    srcid = sub["metadata"]["source"]	# determine the image source id

    imgid = sub["_id"]			# grab the ObjectId corresponding for this image 

    # locate all the classifications of this image by either the experts or volunteers
    if volunteers:
        user_classifications = classifications.find({"subject_ids": imgid, 'expert':{'$exists':False}})
        Nusers = classifications.find({"subject_ids": imgid, 'expert':{'$exists':False}}).count()
        prefix = 'volunteer'
    else:
        user_classifications = classifications.find({"subject_ids": imgid, 'expert':True})
        Nusers = classifications.find({"subject_ids": imgid, 'expert':True}).count()
        prefix = 'expert'


    # loop over the number of classifications
    if Nusers > 0:			# the number of classifications should equal the number of users who classified

        classfile2 = open('%s/RGZ-%s-%s-classifications.txt' % (dat_dir,prefix,srcid), 'w')
    
        # initialise coordinate variables
        radio_ra = []
        radio_dec = []
        radio_x = []
        radio_y = []
        radio_w = []
        radio_h = []
        ir_ra = []
        ir_dec = []
        ir_radius = []
        ir_x = []
        ir_y = []
        
        radio_comp = []
        ir_comp = []

        all_radio = []
        all_radio_markings = []

        Nuser_id = 0			# User id number
        
        #---------------------------------------------------------------------------------------------------------------------
        #---START: loop through the users who classified the image
        for classification in list(user_classifications):
            compid = 0				                # Component id per image
            rclass = classification["annotations"]  # For now, analyze only the first set of continuous regions selected. 
                                                    # Note that last two fields in annotations are timestamp and user_agent

            Nuser_id += 1			                # Increase the number of users who classified by 1.
            
            #-------------------------------------------------------------------------------------------------------------------
            #---START: loop through the keys in the annotation array, making sure that a classification has been made
            for ann in rclass:
            
                if ann.has_key('started_at') or ann.has_key('finished_at') or ann.has_key('user_agent') or ann.has_key('lang'):
                    continue
                
                Nradio = 0					# counter for the number of radio components per classification
                Nir = 0						# counter for the number of IR components per classification
                
                if (ann.has_key('radio') and ann['radio'] != 'No Contours'):			# get the radio annotations
                    radio = ann["radio"]
                    Nradio = len(radio)				# count the number of radio components per classification
                    '''
                    print 'RADIO:'
                    print radio
                    '''
                
                    compid += 1				# we have a radio source - all components will be id with this number

                    list_radio = []

                    #---------------------------------------------------------------------------------------------------------------
                    #---START: loop through number of radio components in user classification
                    for rr in radio:
                        radio_marking = radio[rr]
                        
                        # Find the location and size of the radio box in pixels
                        
                        list_radio.append('%.6f' % float(radio_marking['xmax']))
                        all_radio_markings.append(radio_marking)

                        print >> classfile2, Nuser_id, compid,'RADIO', radio_marking['xmin'], radio_marking['xmax'], radio_marking['ymin'], radio_marking['ymax']

                    all_radio.append(tuple(sorted(list_radio)))

                    #---END: loop through number of radio components in user classification
                    #---------------------------------------------------------------------------------------------------------------
                
                    # get IR counterpart
                    irkey = ann.has_key('ir')
                    ir_nosources = True if (irkey and ann['ir'] == 'No Sources') else False
                
                    if (irkey and not ir_nosources):		# get the infrared annotation for the radio classification.
                        ir = ann["ir"]
                        Nir = 1 #len(ir)				# number of IR counterparts.  
                        '''
                        print 'IR:'
                        print ir
                        '''
                        #exit()
                        #jj = 0
                        for ii in ir:
                            ir_marking = ir[ii]
                        
                            # write to annotation file
                            print >> classfile2, Nuser_id, compid, 'IR', float(ir_marking['x']), float(ir_marking['y'])
                            ir_x.append(float(ir_marking['x']))
                            ir_y.append(float(ir_marking['y']))
                
                    else:	# user did not classify an infrared source
                        Nir = 0
                        xir = -99.
                        yir = -99.
                        radiusir = -99.
                        print >> classfile2, Nuser_id, compid, 'IR', xir, yir
                    
                else:	# user did not classify a radio source
                    Nradio = 0
                    Nir = 0
                    # there should always be a radio source, bug in program if we reach this part.
                    if not ann.has_key('radio'):
                        print >> classfile2,'%i No radio source - error in processing on image %s' % (Nuser_id, srcid)
                    elif ann['radio'] == 'No Contours':
                        print >> classfile2,'%i No radio source labeled by user for image %s' % (Nuser_id,srcid)
                    else:
                        print >> classfile2,'Unknown error processing radio source'
                
                radio_comp.append( Nradio )			# add the number of radio components per user source to array.
                ir_comp.append( Nir )				# add the number of IR counterparts per user soruce to array.

        #---END: loop through the users who classified the image
        #---------------------------------------------------------------------------------------------------------------------
        
        # Process the radio markings into unique components

        rlist = [(rr['xmin'],rr['xmax'],rr['ymin'],rr['ymax']) for rr in all_radio_markings]
        if len(rlist) > 0:
            radio_unique = [rlist[0]]
            if len(all_radio_markings) > 1:
                for rr in rlist[1:]:
                    if rr not in radio_unique:
                        radio_unique.append(rr)
            nr = False
        else:
            nr = True
            radio_unique = [(0,0,0,0)]

        # Use a 2-D Gaussian kernel to find the center of the IR sources and plot the analysis images

        if len(ir_x) > 2:
            try:
                xpeak,ypeak,Z,npeaks = find_ir_peak(ir_x,ir_y)
                plot_image(ir_x,ir_y,sub,xpeak,ypeak,Z,npeaks,all_radio,radio_unique,volunteers=volunteers)
            except LinAlgError:
                npeaks = len(ir_x)
        else:
            print 'Length of IR vector was less than 2'
            npeaks = len(ir_x)
            xpeak,ypeak,Z = np.zeros((423,423)),np.zeros((423,423)),np.zeros((423,423))
            plot_image(ir_x,ir_y,sub,xpeak,ypeak,Z,npeaks,all_radio,radio_unique,volunteers=volunteers)
    
        '''
        # Plot analysis images

        npeaks = len(ir_x)
        if npeaks == 0:
            ir_x,ir_y = [0],[0]
        plot_image(ir_x[0],ir_y[0],npeaks,sub,all_radio,radio_unique,no_radio = nr)
        '''
        
        # calculate the median number of components for both IR and radio for each object in image.
        radio_med = np.median(radio_comp)				# median number of radio components
        Ncomp_radio = np.size(np.where(radio_comp == radio_med))	# number of classifications = median number
        ir_med = np.median(ir_comp)					# median number of infrared components
        Ncomp_ir = np.size(np.where(ir_comp == ir_med))		# number of classifications = median number
        
        if verbose:
            print ' '
            print 'Source............................................................................................: %s' % srcid
            print 'Number of %s users who classified the object..................................................: %d' % (prefix,Nusers)
            print '................'
            print 'Number of %s users who classified the radio source with the median value of radio components..: %d' % (prefix,Ncomp_radio)
            print 'Median number of radio components per %s user.................................................: %f' % (prefix,radio_med)
            print 'Number of %s users who classified the IR source with the median value of IR components........: %d' % (prefix,Ncomp_ir)
            print 'Median number of IR components per %s user....................................................: %f' % (prefix,ir_med)
            print ' '

        classfile2.close()

    else:
        print '%ss did not classify subject %s.' % (prefix.capitalize(),sub['zooniverse_id'])
        ir_x,ir_y = 0,0

    return None

def load_rgz_data():

    # Connect to Mongo database
    # Make sure to run mongorestore /path/to/database to restore the updated files
    # mongod client must be running locally
    
    client = MongoClient('localhost', 27017)
    db = client['radio'] 
    
    subjects = db['radio_subjects'] 		# subjects = images
    classifications = db['radio_classifications']	# classifications = classifications of each subject per user
    users = db['radio_users']	# volunteers doing each classification (can be anonymous)

    return subjects,classifications,users

def load_expert_parameters():

    expert_path = '%s/expert' % rgz_dir

    # Note all times should be in UTC (Zulu)
    json_data = open('%s/expert_params.json' % expert_path).read() 
    experts = json.loads(json_data)

    return experts

def run_expert_sample(subjects,classifications):

    expert_zid = open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir).read().splitlines()

    N = 0
    with open('%s/expert/npeaks_ir_expert_all.csv' % (csv_dir),'wb') as f:
        for sub in list(subjects.find({'zooniverse_id':{'$in':expert_zid},'classification_count':{'$gt':0}})):
        
            Nclass = sub["classification_count"]	# number of classifications made per image
            if Nclass > 0:			# if no classifications move to next image
                npeak = find_consensus(sub,classifications)
                print >> f, npeak

            N += 1
            
            # Check progress by printing to screen every 10 classifications
            if not N % 10:
                print N, datetime.datetime.now().strftime('%H:%M:%S.%f')
    
    
    return None

def run_volunteer_sample(subjects,classifications):

    expert_zid = open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir).read().splitlines()

    N = 0
    with open('%s/expert/npeaks_ir_expert_volunteer.csv' % (csv_dir),'wb') as f:
        for sub in list(subjects.find({'zooniverse_id':{'$in':expert_zid},'classification_count':{'$gt':0}})):
        
            Nclass = sub["classification_count"]	# number of classifications made per image
            if Nclass > 0:			# if no classifications move to next image
                npeak = find_consensus(sub,classifications,volunteers=True)
                print >> f, npeak

            N += 1
            
            # Check progress by printing to screen every 10 classifications
            if not N % 10:
                print N, datetime.datetime.now().strftime('%H:%M:%S.%f')
    
    
    return None

'''
########################################
TO DO

Find a way to automatically check whether classifications matched for a particular object among the expert sample.

Maybe by reading the .dat files?

########################################
'''
def sdi(data):

    # Shannon diversity index
    
    def p(n, N):
        """ Relative abundance """
        if n is 0:
            return 0
        else:
            return (float(n)/N) * np.log(float(n)/N)
            
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n is not 0)

def compare_expert_consensus():

    with open('%s/expert/expert_all_first_ids.txt' % rgz_dir,'rb') as f:
        first_ids = [line.rstrip() for line in f]

    exlist = load_expert_parameters()

    ir_array = []
    ex_array = []
    for first in first_ids:
        ir_temp = []
        exlist_arr = []
        for ex in exlist:
            with open('%s/datfiles/expert/%s/RGZBETA2-%s-classifications.txt' % (rgz_dir,ex['expert_user'],first)) as f:
                x = [line.rstrip() for line in f]
                try:
                    last_line = x[-1].split()
                    n_ir = int(last_line[1])
                    x_ir = last_line[3]
                    if x_ir == '-99.0':
                        n_ir = 0
                    ir_temp.append(n_ir)
                    exlist_arr.append(ex['expert_user'])
                except:
                    pass
        
        ex_array.append(exlist_arr)
        ir_array.append(ir_temp)

    '''
    # Plot number of users per galaxy

    excount = [len(x) for x in ex_array]
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    ax2.hist(excount,bins=6,range=(5,11))
    ax2.set_xlim(6,11)
    fig2.show()
    fig2.tight_layout()
    '''

    c = [Counter(i) for i in ir_array]

    fig = plt.figure(1,(15,4))
    ax = fig.add_subplot(111)

    larr = []
    varr = []
    sdi_arr = []
    for idx,cc in enumerate(c):
        l,v = zip(*cc.items())
        larr.append(list(l))
        varr.append(list(v))
        if len(l) > 1:
            sdi_arr.append(sdi(cc)/np.log(len(l)))
        else:
            sdi_arr.append(sdi(cc))

    iarr = []
    sarr = []
    for idx,(l,s) in enumerate(zip(larr,sdi_arr)):
        iarr.append(np.zeros(len(l),dtype=int)+idx)
        sarr.append(np.zeros(len(l),dtype=float)+s)

    iarr = list_flatten(iarr)
    larr = list_flatten(larr)
    varr = list_flatten(varr)
    sarr = list_flatten(sarr)

    zipped = zip(sarr,larr,iarr,varr)
    zipped.sort()
    sarr,larr,iarr,varr = zip(*zipped)

    ikeys = list(OrderedDict.fromkeys(iarr))
    inew = [ikeys.index(ii) for ii in iarr]

    sc = ax.scatter(inew,larr,c=sarr,s=((np.array(varr)+5)**2),cmap = cm.RdBu_r,vmin=0.,vmax =1.0)
    cbar = plt.colorbar(sc)
    cbar.set_label('Normalized Shannon entropy')

    ax.set_xlim(-1,101)
    ax.set_xlabel('Galaxy')
    ax.set_ylabel('Number of IR sources')
    ax.set_aspect('auto')

    lnp = np.array(larr)
    inp = np.array(inew)
    for i in inp:
        if (inp == i).sum() > 1:
            lmin = np.min(lnp[inp == i])
            lmax = np.max(lnp[inp == i])

            ax.plot([i,i],[lmin,lmax],color='black')

    fig.show()
    fig.tight_layout()

    # Now analyzing the radio components. Order can definitely change, even when associated with a single IR source. 

    # Sort first? Sort first on second column, then radio? That would make sure everything agrees...

    
    '''
    OK. So, Kevin doesn't have any classifications in the set recorded. There are 134 classifications by his IP address in the timeframe
    in question, but they're very short (only 18 minutes), and none of them are in the set of 100. Lots of duplicates, too.

    Eight users classified every galaxy in the sample of 100. 
    42jkb, ivywong, enno.middelberg, xDocR, KWillett, stasmanian, akpinska, vrooje

    klmasters only shows up for 25 of them (as expected)        Galaxies done in timeframe:
    Kevin appears in 0                                          139                                   - Claims he was logged in, but zero classifications under username Kevin
                                                                                                        There are 139 galaxies done by someone matching his IP address,
                                                                                                            but none are in the expert sample of 100 (no idea why). 
                                                                                                            Assume we got no useful classification data from him. 


    '''

    return ir_array

def compare_volunteer_consensus(subjects,classifications,users):

    # Just looks at the total number of IR sources per subject as measured by volunteers. 
    # Should be able to get this by querying MongoDB directly.

    with open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir,'rb') as f:
        zooniverse_ids = [line.rstrip() for line in f]

    ir_array = []
    zid_array = []
    usernames_array = []

    experts = load_expert_parameters()
    expert_usernames = [x['expert_user'] for x in experts]

    for zid in zooniverse_ids:
        ir_temp = []
        username_temp = []
        subject_id = subjects.find_one({'zooniverse_id':zid})['_id']
        clist = classifications.find({'subject_ids.0':subject_id})

        zid_array.append(zid)

        for c in clist:

            usernames_bad = expert_usernames

            # Need to remove both duplicates and experts
            try:
                user_name = c['user_name']
                username_temp.append(user_name)

                if (user_name not in usernames_bad) and \
                    (users.find_one({'name':user_name})['projects']['52afdb804d69636532000001']['classification_count'] > 10):

                    # Hate the way annotations are stored. 

                    # Need to remove counts both for metadata and where there's a glitch recording; 
                    # lots of records have a spurious "No Sources" and "No Contours", possibly due to reset

                    annotations = c['annotations']
                    na = len(annotations)
                    for a in annotations:
                        if a.keys()[0] in bad_keys:
                            na -= 1
                        if 'ir' in a.keys():
                            if a['ir'] is 'No Sources' and a['radio'] is 'No Contours':
                                na -= 1

                    ir_temp.append(na)

                    # Records only the FIRST time someone with a particular username classifies it. Not ideal. 
                    usernames_bad.append(user_name)

            # Anonymous users
            except KeyError:
                pass
                '''
                username_temp.append('Anonymous')

                annotations = c['annotations']
                na = len(annotations)
                for a in annotations:
                    if a.keys()[0] in bad_keys:
                        na -= 1
                    if 'ir' in a.keys():
                        if a['ir'] is 'No Sources' and a['radio'] is 'No Contours':
                            na -= 1

                ir_temp.append(na)
                '''

        #print zid,ir_temp
        ir_array.append(ir_temp)
        usernames_array.append(username_temp)
        
    i_nozeros = []
    for ii,zz,uu in zip(ir_array,zid_array,usernames_array):
        if len(ii) > 0:
            i_nozeros.append(ii)
        '''
        else:
            print 'No non-expert classifications for %s' % zz,uu
        '''

    c = [Counter(i) for i in i_nozeros]

    fig = plt.figure(2,(15,4))
    fig.clf()
    ax = fig.add_subplot(111)

    larr = []
    varr = []
    sdi_arr = []
    for idx,cc in enumerate(c):
        l,v = zip(*cc.items())
        larr.append(list(l))
        varr.append(list(v))
        if len(l) > 1:
            sdi_arr.append(sdi(cc)/np.log(len(l)))
        else:
            sdi_arr.append(sdi(cc))

    iarr = []
    sarr = []
    for idx,(l,s) in enumerate(zip(larr,sdi_arr)):
        iarr.append(np.zeros(len(l),dtype=int)+idx)
        sarr.append(np.zeros(len(l),dtype=float)+s)

    iarr = list_flatten(iarr)       # Index of galaxy image
    larr = list_flatten(larr)       # Number of IR sources
    varr = list_flatten(varr)       # Number of users who selected the given number of IR sources
    sarr = list_flatten(sarr)       # Shannon diversity index

    zipped = zip(sarr,larr,iarr,varr)
    zipped.sort()
    sarr,larr,iarr,varr = zip(*zipped)

    ikeys = list(OrderedDict.fromkeys(iarr))
    inew = [ikeys.index(ii) for ii in iarr]

    sc = ax.scatter(inew,larr,c=sarr,s=((np.array(varr)+5)**2),cmap = cm.RdBu_r,vmin=0.,vmax =1.0)
    cbar = plt.colorbar(sc)
    cbar.set_label('Normalized Shannon entropy')

    ax.set_xlim(-1,len(c)+1)
    ax.set_title('RGZ volunteer classifications')
    ax.set_xlabel('Galaxy')
    ax.set_ylabel('Number of IR sources')
    ax.set_aspect('auto')

    lnp = np.array(larr)
    inp = np.array(inew)

    # Draw black lines between dots
    for i in inp:
        if (inp == i).sum() > 1:
            lmin = np.min(lnp[inp == i])
            lmax = np.max(lnp[inp == i])

            ax.plot([i,i],[lmin,lmax],color='black')

    fig.show()
    #fig.tight_layout()

    return None

def expert_vs_volunteer():

    # Direct comparison of the expert vs. volunteer classifications for all galaxies?

    return None

def histogram_experts(classifications,users):

    # Goal: find the distribution and average number of IR sources per image for the science team

    do_experts = True

    if do_experts:
        experts = load_expert_parameters()

        # Add Larry's regular username as well as Ray Norris

        experts.append({'expert_user':'DocR'})
        experts.append({'expert_user':'raynorris'})

        expert_avg = []
        for ex in experts:

            username = ex['expert_user']
            classcount = classifications.find({'user_name':username}).count()

            if classcount > 0:
                c = classifications.find({'user_name':username})
                nir = []
                for cc in list(c):
                    annotations = cc['annotations']
                    na = len(annotations)
                    for a in annotations:
                        if a.keys()[0] in bad_keys:
                            na -= 1
                        if 'ir' in a.keys():
                            if a['ir'] is 'No Sources' and a['radio'] is 'No Contours':
                                na -= 1

                    nir.append(na)

                print '%20s averages %.2f IR sources per image over %i classifications' % (username,np.mean(nir),classcount)

                expert_avg.append(np.mean(nir))

        print '-----------------------------------'

    # Now look at the volunteers

    rgz_id = classifications.find_one()['project_id']

    # Odd that about half of the users don't seem to have a classification count for RGZ. Is that actually true?

    '''
    In [80]: users.find({'projects.%s.classification_count' % rgz_id:{'$exists':True}}).count()
    Out[80]: 4907

    In [79]: users.find({'projects.%s.classification_count' % rgz_id:{'$exists':False}}).count()
    Out[79]: 3312
    
    All users with a classification count do have at least one classification. 
    In the second group, though, most have zero, but some have a couple classifications (maximum of 6)

        742 have at least one classification
        2570 have no classifications

    So we actually have only 4907+742 = 5,649 contributing users, rather than the 8,219 people in the users db and the 4,955 listed on the API

    '''

    # Concatenate the two groups

    users_good = list(users.find({'projects.%s.classification_count' % rgz_id:{'$exists':True}}))
    users_unsure = users.find({'projects.%s.classification_count' % rgz_id:{'$exists':False}})
    for u in list(users_unsure):
        if classifications.find({'user_id':u['_id']}).count() > 0:
            users_good.append(u)

    nir_username_volunteers = []
    nir_avg_volunteers = []
    nir_count_volunteers = []
    for u in users_good:

        classcount = classifications.find({'user_id':u['_id']}).count()

        if classcount > 0:
            c = classifications.find({'user_id':u['_id']})
            nir = []
            for cc in list(c):
                annotations = cc['annotations']
                na = len(annotations)
                for a in annotations:
                    if a.keys()[0] in bad_keys:
                        na -= 1
                    if 'ir' in a.keys():
                        if a['ir'] is 'No Sources' and a['radio'] is 'No Contours':
                            na -= 1

                nir.append(na)

            
            avg = np.mean(nir)

            #print '%20s averages %.2f IR sources per image over %i classifications' % (u['name'],avg,classcount)

        else:   # Shouldn't happen with this list
            print 'No classifications found for %s' % u['name']
            avg = 0.


        nir_username_volunteers.append(u['name'])
        nir_avg_volunteers.append(avg)
        nir_count_volunteers.append(classcount)

    # If we eliminate users who average more than two IR sources per image, how much of the data would that reject?

    # counts_lost = np.sum([c for a,b,c in zip(nir_username_volunteers, nir_avg_volunteers, nir_count_volunteers) if b > 2.0])
    # Only 413 classifications. Negligible.

    return nir_username_volunteers, nir_avg_volunteers, nir_count_volunteers, expert_avg

def plot_histogram_experts(names, avgs, counts, expert_avg):

    xpairs = [[x,x] for x in expert_avg]
    xlist = []
    for xends in xpairs:
        xlist.extend(xends)
        xlist.append(None)

    avg_cutoff = np.linspace(0,4,50)
    frac_lost = []
    for ac in avg_cutoff:
        frac_lost.append(np.sum([c for a,c in zip(avgs,counts) if a > ac])/float(sum(counts)))

    # Plot results

    fig = plt.figure(1)
    fig.clf()

    ax1 = fig.add_subplot(221)
    ax1.scatter(avgs,counts,color='black',marker='.',s=1,alpha=0.5)
    ax1.set_xlabel('Mean IR sources per image')
    ax1.set_ylabel('Total number of classifications')
    ax1.set_yscale('log')

    ax2 = fig.add_subplot(222)
    ax2.plot(avg_cutoff,frac_lost,color='green',lw=3)
    ax2.set_ylim(-0.02,1.02)
    ax2.set_xlabel('Cutoff for IR sources/image')
    ax2.set_ylabel('Fraction of data affected')

    ax3 = fig.add_subplot(223)
    ax3.hist(avgs,bins=np.linspace(0,4,100))
    ax3.text(2.5,700,'All',fontsize=16)
    ax3.set_xlabel('Mean IR sources per image')
    ax3.set_ylabel('Count (users)')

    ax4 = fig.add_subplot(224)
    ax4.hist(np.array(avgs)[np.array(counts) > 10],bins=np.linspace(0,4,100),color='cyan')
    ax4.text(2.5,250,r'$N_{class}>10$',fontsize=16)
    ax4.set_xlabel('Mean IR sources per image')
    ax4.set_ylabel('Count (users)')
    ax4.set_xlim(ax3.get_xlim())

    for ax in (ax2,ax3,ax4):
        ypairs = [ax.get_ylim() for x in range(len(xpairs))]
        ylist = []
        for yends in ypairs:
            ylist.extend(yends)
            ylist.append(None)

        ax.plot(xlist,ylist,color='red',alpha=0.5)

    fig.show()

    # Save hard copy of the figure
    fig.savefig('%s/plots/histogram_avg_ir_sources.png' % rgz_dir)

    return None

def update_experts(classifications,experts): 

    for ex in experts:
        expert_dates = (dateutil.parser.parse(ex['started_at']),dateutil.parser.parse(ex['ended_at']))
        classifications.update({"updated_at": {"$gt": expert_dates[0],"$lt":expert_dates[1]},"user_name":ex['expert_user']},{'$set':{'expert':True}},multi=True)

    return None


########################################
# Call program from the command line
########################################

if __name__ == '__main__':

    plt.ioff()

    subjects,classifications,users = load_rgz_data()
    experts = load_expert_parameters()
    update_experts(classifications,experts)

    # Experts on sample of 100
    run_expert_sample(subjects,classifications)
    #plot_npeaks()
    # Volunteers on sample of 100
    run_volunteer_sample(subjects,classifications)
    #plot_npeaks(volunteers=True)

