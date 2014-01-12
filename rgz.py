# import necessary python packages
import numpy as np
from matplotlib import pyplot as plt
import datetime
import pandas as pd
import os

from pymongo import MongoClient

from astropy.io import fits
from astropy import wcs
#from astropy import coordinates as coord
#from astropy.io import votable

#------------------------------------------------------------------------------------------------------------
# Define some functions

#http://radio.galaxyzoo.org.s3.amazonaws.com/beta/subjects/radio/S334.jpg
#http://radio.galaxyzoo.org.s3.amazonaws.com/beta/subjects/standard/S334.jpg
#http://radio.galaxyzoo.org.s3.amazonaws.com/beta/subjects/raw/S334.fits.gz

# (1) Determine the WCS object based on subject
def getWCSObj(subject):
  src = subject["metadata"]["source"]
  path = "./IMGS/%s.fits" % src
  hdulist = fits.open(path)
  w = wcs.WCS(hdulist[0].header)
  return w
  

def plot_user_counts(df):
    
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(211)

    volunteers = pd.value_counts(df.user_name)

    # Calculate number of anonymous users and include in data

    anonymous_count = df._id.count() - df.user_name.count()
    volunteers = volunteers.set_value("anonymous", anonymous_count);

    vcplot = volunteers.plot(ax=ax1,use_index=True)

    vcplot.set_title('RGZ volunteer distribution')
    vcplot.set_ylabel('Number of classifications')
    vcplot.set_ylim((1,1e5))
    vcplot.set_yscale('log')

    ax2 = fig.add_subplot(212)

    vchist = volunteers.hist(ax=ax2,bins=100,bottom=0.1)

    vchist.set_ylabel('Classifications per volunteer')
    vchist.set_xlabel('Number of classifications')
    vchist.set_yscale('log')

    fig.show()
    fig.tight_layout()

    # Save hard copy of the figure
    fig.savefig('%s/plots/classifications_per_user.png')

    return None



#------------------------------------------------------------------------------------------------------------

# Set constants
beta_release_date = datetime.datetime(2013, 10, 20, 12, 0, 0, 0)	# date of beta release (YYY,MM,DD,HH,MM,SS,MS)
main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)

IMG_HEIGHT = 424.0			# number of pixels in the JPG image along the y axis
IMG_WIDTH = 424.0			# number of pixels in the JPG image along the x axis
FITS_HEIGHT = 301.0			# number of pixels in the FITS image along the y axis
FITS_WIDTH = 301.0			# number of pixels in the FITS image along the x axis
PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image

xjpg2fits = float(IMG_WIDTH/FITS_WIDTH)		# map the JPG pixels to the FITS pixels in x
yjpg2fits = float(IMG_HEIGHT/FITS_HEIGHT)	# map the JPG pixels to the FITS pixels in y

# Connect to Mongo database
# Make sure to run mongorestore /path/to/database to restore the updated files
# mongod client must be running locally

client = MongoClient('localhost', 27017)
db = client['ouroboros'] 

subjects = db['radio_subjects'] 		# subjects = images
classifications = db['radio_classifications']	# classifications = classifications of each subject per user
users = db['radio_users']	# volunteers doing each classification (can be anonymous)

# Retrieve RGZ data since the launch data
batch = classifications.find({"updated_at": {"$gt": main_release_date}})
df = pd.DataFrame( list(batch) )

# Plot number of classifications per volunteer
plot_user_counts(df)

# Get some quick statistics on the dataset so far
n_subjects = subjects.count()		# determine the number of images in the data set
n_classifications = classifications.find({"updated_at": {"$gt": main_release_date}}).count() # total number of classifications
n_users = users.count()

# Find the most recent classification in this data dump
mrc = db.radio_classifications.find().sort([("updated_at", -1)]).limit(1)
most_recent_date = [x for x in mrc][0]['updated_at']

# Find number of anonymous classifications
anonymous_count = df._id.count() - df.user_name.count()
total_count = df.user_name.count()
anonymous_percent = float(anonymous_count)/total_count * 100

print ' '
print 'RGZ data as of %s' % most_recent_date.strftime("%H:%M:%S%Z %b %d, %Y")
print '---------------------------------'
print 'Total classifications   : %i' % n_classifications
print 'Total distinct subjects : %i' % n_subjects
print 'Total distinct users    : %i' % n_users
print ' '
print 'Percent of classifications by anonymous users: %.1f (%i/%i)' % (anonymous_percent,anonymous_count,total_count)
print ' '

'''


classfile = open('RGZBETA2-class.dat', 'w')
#compfile = open('RGZBETA-components.dat','w')
#print >> compfile,'{0:10} {1:7} {2:7} {3:11} {4:7} {5:11}'.format('#SRC_ID   ','N_TOTAL','N_RADIO','N_RADIO_MED','   N_IR','   N_IR_MED')

#------------------------------------------------------------------------------------------------------------------------
#---START: loop through images/subjects
for item in list(subjects.find().limit(100)):
  #print item
  #print '------------------------'
  Nclass = item["classification_count"]	# number of classifications made per image

  if Nclass <= 0:			# if no classifications move to next image
    continue

  srcid = item["metadata"]["source"]	# determine the image source id

  classfile2 = open('datfiles/RGZBETA2-%s-classifications.dat' % srcid, 'w')

  #wgetcmd = 'wget http://radio.galaxyzoo.org.s3.amazonaws.com/beta/subjects/raw/' + srcid + '.fits.gz'
  #os.system(wgetcmd)
  #mvcmd = 'mv '+ srcid + '.fits.gz ./IMGS/'
  #os.system(mvcmd)
  #zipcmd = 'gunzip ./IMGS/%s.fits.gz' % srcid
  #os.system(zipcmd)

  #wcsObj = getWCSObj(item)		# get WCS object for pixel to sky transformation
 
  #zipcmd = 'gzip ./IMGS/%s.fits' % srcid
  #os.system(zipcmd)
 
  print >> classfile, srcid, Nclass
  print >> classfile, '-----------'
  annfile = open('annfiles/%s.ann' % srcid, 'w')		# kvis annotation file
  print >> annfile, '#KARMA annotation file'
  print >> annfile, 'COORD W'
  print >> annfile, 'PA STANDARD'
  print >> annfile, 'FONT hershey14'
  
  imgid = item["_id"]			# grab the ObjectId corresponding for this image 
  # locat all the classifications of this image by user
  user_classifications = classifications.find({"subject_ids": imgid, "updated_at": {"$gt": beta_release_date}})
  # count the number of users who classified this object
  Nusers = classifications.find({"subject_ids": imgid, "updated_at": {"$gt": beta_release_date}}).count()
  # loop over the number of classifications
  if Nclass == Nusers:			# the number of classifications should equal the number 
                                        # of users who classified
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
    # loop over users who classified
    print 'Working on source: %s' % srcid
    Nuser_id = 0			# User id number


    #---------------------------------------------------------------------------------------------------------------------
    #---START: loop through the users who classified the image
    for k in list(user_classifications):
      compid = 0				# Component id per image
      c = k["annotations"][0]			# get the user annotations
      print >> classfile, c
      Nuser_id = Nuser_id + 1			# increase the number of user who classified by 1.

      #-------------------------------------------------------------------------------------------------------------------
      #---START: loop through the keys in the annotation array, making sure that a classification has been made
      for key in c:
        if key in ["started_at", "finished_at"]:	# no user classifications, move to next classification
          continue
   
        markings = c[key]				# get keys related to the user annotations

        print 'MARKINGS:'
        print markings

        Nradio = 0					# counter for the number of radio componets per classification
        Nir = 0						# counter for the number of IR components per classification

        try:
          keys = markings.has_key("radio")
          if markings.has_key("radio"):			# get the radio annotations
            radio = markings["radio"]
            Nradio = len(radio)				# count the number of radio components per classification
            print 'RADIO:'
            print radio
            compid = compid + 1				# we have a radio source - all componets will be id with this number

            #---------------------------------------------------------------------------------------------------------------
            #---START: loop through number of radio components in user classification
            for rr in radio:
              radio_marking = radio[rr]
              
              # NOTE: JPG images begin at (0,0) in the upper right of the image.  FITS images begin at (1,1) in 
              #       the lower left of the image.  I have added 1.0 to the x and y values of the annotation.
              #       The pixels in the JPG image increaser to the right for x and down for y.
              xULC = 1.0 + float(radio_marking["x"])	# x pixel position in the JPG file for the ULC of the box
                                                  	# surrounding the radio blob.  Usually the 3 sigma contour.
              yULC = 1.0 + float(radio_marking["y"])	# y pixel position in the JPG file for the ULC of the box
                                                  	# surrounding the radio blob.  Usually the 3 sigma contour.
              yULC = (IMG_HEIGHT - yULC)			# transfer the pixel value to increase upwards in y in JPG pixels.
              width = float(radio_marking["width"])	# width of the box surrounding the radio blob in JPG pixels.
              height = float(radio_marking["height"])	# height of the box surrounding the radio blob in JPG pixels.

              # convert the pixels into FITS image pixels
              xULC = xULC/xjpg2fits			# x pixel position in the FITS image for the ULC of the box
		  					# surrounding the radio blob.
              yULC = yULC/yjpg2fits			# y pixel position in the FITS image for the ULC of the box
		  					# surrounding the radio blob.
              width = width/xjpg2fits			# width of the box surround the radio blob in FITS pixels.
              height = height/yjpg2fits			# height of the box surrounding the radio blob in FITS pixels.

              # Want to create a kvis annotation file to overlay on source
              # the radio will be given by the box and the IR will be shown in a circle
              #
              # BOX [coord type] <x1> <y1> <x2> <y2>
              # CIRCLE [coord type] <x_cent> <y_cent> <radius> 
#              coords1 = wcsObj.wcs_pix2world([[xULC,yULC]], 0)[0]	# convert ULC pixels to RA and DEC (degrees)
#              xBRC = xULC + width				# x pixel position in the FITS image for the BRC of the box
#		  					# surrounding the radio blobe.
#              yBRC = yULC - height			# y pixel position in the FITS image for the BRC of the box
#		  					# surrounding the radio blobe.
#              coords2 = wcsObj.wcs_pix2world([[xBRC,yBRC]], 0)[0] # convert BRC pixels to RA and DEC (degrees)
#              # write to annotation file
#              print >> annfile, 'COLOUR BLUE'
#              print >> annfile, 'BOX',coords1[0],coords1[1],coords2[0],coords2[1]
#              print >> classfile2, Nuser_id, compid,'RADIO', xBRC, yBRC, width, height, coords1[0],coords1[1]

            #---END: loop through number of radio components in user classification
            #---------------------------------------------------------------------------------------------------------------

            # get IR counterpart
            if markings.has_key("infrared"):		# get the infrared annotation for the radio classification.
              ir = markings["infrared"]
              Nir = 1 #len(ir)				# number of IR counterparts.  NOTE: for beta this is 1 but
		  					# future Galaxy Zoo: Radio releases this will be changed.
              print 'IR:'
              print ir
              #exit()
              #jj = 0
              for ii in ir:
                ir_marking = ir[ii]
                
                #if jj == 0: yir = 1.0 + float(ir_marking)	# y pixel location in the JPG image for the IR source.
                #if jj == 1: xir = 1.0 + float(ir_marking)	# x pixel location in the JPG image for the IR source.
                #if jj == 2: radiusir = float(ir_marking)	# radius of circle centred at xir,yir.
                #jj = jj+1
             
              # convert the pixels into FITS image pixels
              #xir = xir/xjpg2fits				# x pixel location in FITS image for IR source.
              #yir = (IMG_HEIGHT - yir)/xjpg2fits		# y pixel location in FITS image for IR source.
              #radiusir = (radiusir/xjpg2fits) * PIXEL_SIZE # radius of circle around IR source in degrees in FITS image.
		  					# NOTE: radius removed in following releases of project.

              #coords3 = wcsObj.wcs_pix2world([[xir,yir]], 0)[0] # convert IR pixel location to RA and DEC (degrees)
              # write to annotation file
              #print >> annfile, 'COLOUR RED'
              #print >> annfile, 'CIRCLE',coords3[0],coords3[1],radiusir
              #print >> annfile, 'CROSS',coords3[0],coords3[1],0.0005,0.0005
              #print >> classfile2, Nuser_id, compid, 'IR', xir, yir, radiusir, coords3[0],coords3[1]

            else:	# user did not classify an infrared source
              Nir = 0
              xir = -99.
              yir = -99.
              radiusir = -99.
              #print >> classfile2, Nuser_id, compid, 'IR', xir, yir, radiusir, -99.0, -99.0
            
        except AttributeError:	# user did not classifiy a radio source
          Nradio = 0
          Nir = 0
          # there should always be a radio source, bug in program if we reach this part.
          print ''
          print >> classfile2,'No radio source?...error on image %s' % srcid
          #print 'exiting.....'
          #exit()
        
        radio_comp.append( Nradio )			# add the number of radio components per user source to array.
        ir_comp.append( Nir )				# add the number of IR counterparts per user soruce to array.
    #---END: loop through the users who classified the image
    #---------------------------------------------------------------------------------------------------------------------

  else:
    print 'number of users who classified does not equal classification count?'
#    print 'exiting....'
    


  # calculate the median number of components for both IR and radio for each object in image.
  #radio_med = np.median(radio_comp)				# median number of radio components
  #Ncomp_radio = np.size(np.where(radio_comp == radio_med))	# number of classifications = median number
  #ir_med = np.median(ir_comp)					# median number of infrared components
  #Ncomp_ir = np.size(np.where(ir_comp == ir_med))		# number of classifications = median number

#  print >> compfile,'{0:10} {1:7d} {2:7d} {3:11.1f} {4:7d} {5:11.1f}'.format(srcid,len(radio_comp),Ncomp_radio,radio_med,Ncomp_ir,ir_med)
  annfile.close()
  classfile2.close()
#  print >> classfile2,' '
#  print >> classfile2,'Source.....................................................................................: %s' % srcid
#  print >> classfile2,'Number of users who classified the object..................................................: %d' % len(radio_comp)
#  print >> classfile2,'Number of users who classified the radio source with the median value of radio components..: %d' % Ncomp_radio
#  print >> classfile2,'Median value of radio components from the users............................................: %f' % radio_med
#  print >> classfile2,'Number of users who classified the IR source with the median value of IR components........: %d' % Ncomp_ir
#  print >> classfile2,'Median value of IFcomponents from the users................................................: %f' % ir_med

#---END: loop through images/subjects
#------------------------------------------------------------------------------------------------------------------------

classfile.close()
#compfile.close()
'''
