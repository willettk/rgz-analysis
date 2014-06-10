# import necessary python packages
import numpy as np
import pandas as pd
import datetime
import os
import urllib
import cStringIO

from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from collections import Counter
from matplotlib import pyplot as plt
from pymongo import MongoClient
from astropy.io import fits
from astropy import wcs
from scipy import stats
from scipy import optimize
from scipy.linalg.basic import LinAlgError
from PIL import Image
#from astropy import coordinates as coord
#from astropy.io import votable

#------------------------------------------------------------------------------------------------------------

# Setup path locations

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
csv_dir = '%s/csv' % rgz_dir

plot_dir = '../plots/expert' % expert_user
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

ir_peaks_plot_dir = '%s/ir_peaks' % plot_dir
if not os.path.isdir(ir_peaks_plot_dir):
    os.mkdir(ir_peaks_plot_dir)

dat_dir = '../datfiles/expert' % expert_user
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

def plot_npeaks(expert_user):

    # Read in data

    with open('%s/npeaks_ir_expert_%s.csv' % (csv_dir,expert_user),'rb') as f:
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
    fig.savefig('%s/%s/ir_peaks_histogram_expert_%s.png' % (plot_dir,expert_user,expert_user))

    plt.close()

    return None

def find_ir_peak(x,y,srcid):

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

def plot_image(x,y,npeaks,sub,all_radio,radio_unique,no_radio=False):

    '''
    # Find the peak

    xpeak = X[Z==Z.max()][0]
    ypeak = Y[Z==Z.max()][0]
    '''

    # Plot the infrared results
    
    fig = plt.figure(1,(12,4))
    ax = fig.add_subplot(133,aspect='equal')

    '''
    # Plot the KDE map
    ax.imshow(np.rot90(Z), cmap=plt.cm.hot_r,extent=[xmin, xmax, ymin, ymax])
    # Plot the individual sources
    ax.plot(x, y, 'go', markersize=4)
    '''

    ax.text(270,40,r'IR peak: $(%i,%i)$'%(x,y),color='k',fontsize=14)
    ax.text(270,70,r'$N_{peaks}$ = %i' % npeaks,color='k',fontsize=14)
    if npeaks > 0:
        ax.plot(x,y,'r*',markersize=12)

    # Plot the radio counts

    if not no_radio:
        radio_flattened = [item for sublist in all_radio for item in sublist]
        uniques = set(radio_flattened)
        d = dict(zip(uniques,np.arange(len(uniques))))
        c = Counter(all_radio)

        for idx,ckeys in enumerate(c.keys()):
            if len(ckeys) > 1:
                t = ' and R'.join([str(d[x]) for x in ckeys])
            else:
                t = d[ckeys[0]]
            singular = 's' if c[ckeys] != 1 else ''
            ax.text(150,400-idx*20,'%3i vote%s: R%s' % (c[ckeys],singular,t))

        # Rectangle showing the radio box size

        radio_ir_scaling_factor = 435./132

        box_counts = Counter(radio_flattened)
        for ru in radio_unique:

            x0,x1,y0,y1 = [float(ru_) * radio_ir_scaling_factor for ru_ in ru]

            # Assume xmax matching is still good

            xmax_index = '%.6f' % float(ru[1])
            component_number = d[xmax_index]
            number_votes = box_counts[xmax_index]
            rectangle = plt.Rectangle((x0,y0), x1-x0, y1-y0, fill=False, linewidth=number_votes/2., edgecolor = 'b')
            ax.add_patch(rectangle)
            ax.text(x0-15,y0-15,'R%s' % component_number)
    else:
        ax.text(270,100,'No radio components',color='k',fontsize=14)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymax, ymin])
    ax.set_title('%s (%s)' % (sub['metadata']['source'],sub['zooniverse_id']))
    #fig.show()

    # Display IR and radio images

    url_standard = sub['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax2 = fig.add_subplot(131)
    ax2.imshow(im_standard,origin='upper')

    url_radio = sub['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax3 = fig.add_subplot(132)
    ax3.imshow(im_radio,origin='upper')

    # Save hard copy of the figure
    fig.savefig('%s/ir_peaks/%s_ir_peak.png' % (plot_dir,sub['zooniverse_id']))

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    plt.close()

    return None
    
def find_consensus(sub,classifications,expert_user,expert_dates,verbose=False):
    
    Nclass = sub["classification_count"]	# number of classifications made per image
    srcid = sub["metadata"]["source"]	# determine the image source id

    imgid = sub["_id"]			# grab the ObjectId corresponding for this image 

    # locate all the classifications of this image by user
    user_classifications = classifications.find({"subject_ids": imgid, "updated_at": {"$gt": expert_dates[0],"$lt":expert_dates[1]},"user_name":expert_user})
    # count the number of users who classified this object
    Nusers = classifications.find({"subject_ids": imgid, "updated_at": {"$gt": expert_dates[0],"$lt":expert_dates[1]},"user_name":expert_user}).count()

    # loop over the number of classifications
    if Nusers is not None:			# the number of classifications should equal the number of users who classified

        classfile2 = open('%s/RGZBETA2-%s-classifications.txt' % (dat_dir,srcid), 'w')
    
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

        # Plot analysis images

        npeaks = len(ir_x)
        if npeaks == 0:
            ir_x,ir_y = [0],[0]
        plot_image(ir_x[0],ir_y[0],npeaks,sub,all_radio,radio_unique,no_radio = nr)
        
        # calculate the median number of components for both IR and radio for each object in image.
        radio_med = np.median(radio_comp)				# median number of radio components
        Ncomp_radio = np.size(np.where(radio_comp == radio_med))	# number of classifications = median number
        ir_med = np.median(ir_comp)					# median number of infrared components
        Ncomp_ir = np.size(np.where(ir_comp == ir_med))		# number of classifications = median number
        
        if verbose:
            print ' '
            print 'Source.....................................................................................: %s' % srcid
            print 'Number of users who classified the object..................................................: %d' % len(radio_comp)
            print '................'
            print 'Number of users who classified the radio source with the median value of radio components..: %d' % Ncomp_radio
            print 'Median number of radio components per user.................................................: %f' % radio_med
            print 'Number of users who classified the IR source with the median value of IR components........: %d' % Ncomp_ir
            print 'Median number of IR components per user....................................................: %f' % ir_med
            print ' '

        classfile2.close()

    else:
        print 'Expert (s) did not classify subject %s.' % (expert_user,sub['zooniverse_id'])
      
    return npeaks

def load_rgz_data():

    # Connect to Mongo database
    # Make sure to run mongorestore /path/to/database to restore the updated files
    # mongod client must be running locally
    
    client = MongoClient('localhost', 27017)
    db = client['ouroboros'] 
    
    subjects = db['radio_subjects'] 		# subjects = images
    classifications = db['radio_classifications']	# classifications = classifications of each subject per user
    users = db['radio_users']	# volunteers doing each classification (can be anonymous)

    return subjects,classifications,users

def run_expert_sample(subjects,classifications,expert_user,expert_dates):

    expert_zid = open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir).read().splitlines()

    N = 0
    with open('%s/npeaks_ir_expert_%s.csv' % (csv_dir,expert_user),'wb') as f:
        for sub in list(subjects.find({'zooniverse_id':{'$in':expert_zid},'classification_count':{'$gt':0}})):
        
            Nclass = sub["classification_count"]	# number of classifications made per image
            if Nclass > 0:			# if no classifications move to next image
                npeak = find_consensus(sub,classifications,expert_user,expert_dates)
                print >> f, npeak

            N += 1
            
            # Check progress by printing to screen every 10 classifications
            if not N % 10:
                print N, datetime.datetime.now().strftime('%H:%M:%S.%f')
    
    
    return None

# If program is called from the command line, process the full dataset

if __name__ == '__main__':

    # Larry Rudnick
    expert_user = 'xDocR'
    expert_dates = (datetime.datetime(2014, 06, 07, 17, 43, 0, 0),datetime.datetime(2014, 06, 07, 18, 30, 0, 0))

    # Kyle Willett
    # expert_user = 'KWillett'
    # expert_dates = (datetime.datetime(2014, 06, 9, 04, 00, 0, 0),datetime.datetime(2014, 06, 9, 04, 45, 0, 0))

    subjects,classifications,users = load_rgz_data()
    run_expert_sample(subjects,classifications,expert_user,expert_dates)
    plot_npeaks(expert_user)

