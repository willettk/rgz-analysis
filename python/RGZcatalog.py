'''
RGZcatalog is a pipeline that takes all of the completed RGZ subjects and creates a Mongo database containing
consensus matching information, radio morphology, IR counterpart location, and data from corresponding
AllWISE and SDSS catalogs.
'''

import logging, urllib2, time, argparse, json, os, datetime
import pymongo
import numpy as np
import StringIO, gzip
from astropy.io import fits
from astropy import wcs, coordinates as coord, units as u

#custom modules for the RGZ catalog pipeline
import catalog_functions as fn #contains miscellaneous helper functions
import processing as p #contains functions that process the data
from update_consensus_csv import updateConsensus #replaces the current consensus collection with a specified csv
from find_duplicates import find_duplicates #finds and marks any radio components that are duplicated between sources

from consensus import rgz_path, data_path, db
in_progress_file = '%s/subject_in_progress.txt' % rgz_path

def RGZcatalog():
    
    #start timer
    starttime = time.time()
    
    #begin logging even if not run from command line
    logging.basicConfig(filename='RGZcatalog.log', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logging.captureWarnings(True)

    #check if consensus collection needs to be updated; if so, drop entire consensus collection and replace it with entries from the designated CSV file
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--consensus', help='replace the current consensus collection with a specified csv')
    args = parser.parse_args()
    if args.consensus:
        updateConsensus(args.consensus)
    
    #connect to database of subjects
    #logging.info('Connecting to MongoDB')
    #db = pymongo.MongoClient()['radio']
    subjects = db['radio_subjects']
    consensus = db['consensus']
    catalog = db['catalog'] #this is being populated by this program
    if catalog.count():
        logging.info('Catalog contains entries; appending')
    else:
        catalog.create_index([('catalog_id', 1)])
    
    #get dictionary for finding the path to FITS files and WCS headers
    with open('%s/first_fits.txt' % rgz_path) as f:
        lines = f.readlines()
    pathdict = {}
    for l in lines:
        spl = l.split(' ')
        pathdict[spl[1].strip()] = '%s/rgz/raw_images/RGZ-full.%i/FIRST-IMGS/%s.fits' % (data_path, int(spl[0]), spl[1].strip())

    #count the number of entries from this run and how many entries are in the catalog total
    count = 0
    if catalog.count() != 0:
        for entry in catalog.find().sort('catalog_id', -1).limit(1):
            IDnumber = entry['catalog_id']
    else:
        IDnumber = 0

    #find completed catalog entries so they can be skipped
    consensus_set = set()
    for source in consensus.find():
        consensus_set.add(source['zooniverse_id'])
    catalog_set = set()
    for entry in catalog.find():
        catalog_set.add(entry['zooniverse_id'])
    to_be_completed = consensus_set.difference(catalog_set)
    if os.path.exists(in_progress_file):
        with open(in_progress_file, 'r') as f:
            in_progress_zid = f.read()
        to_be_completed.discard(in_progress_zid)
    to_be_completed = list(to_be_completed)
    
    #iterate through all noncompleted subjects
    for subject in subjects.find({'zooniverse_id': {'$in':to_be_completed} }).batch_size(10):
    #for subject in subjects.find({'zooniverse_id': {'$in': ['ARG00000sl', 'ARG0003f9l']} }):
    #for subject in subjects.find({'zooniverse_id':'ARG00000sl'}): #sample subject with distinct sources
    #for subject in subjects.find({'zooniverse_id':'ARG0003f9l'}): #sample subject with multiple-component source
        
        #mark subject as being in-progress
        with open(in_progress_file, 'w') as f:
            f.write(subject['zooniverse_id'])
        
        #iterate through all consensus groupings
        for source in consensus.find({'zooniverse_id':subject['zooniverse_id'], 'first_id':{'$exists':True}}):
            
            #do not process if this object in this source is already in the catalog
            process = True
            for i in catalog.find({'zooniverse_id':subject['zooniverse_id']}):
                if i['consensus']['label'] == source['label']:
                    process = False
            
            if process:
                
                logging.info('Processing consensus object %s within subject field %s', source['label'], subject['zooniverse_id'])
                
                count += 1
                IDnumber += 1
                
                #display which entry is being processed to see how far the program is
                print 'Processing entry %i (consensus %s in subject %s)' % (IDnumber, source['label'], subject['zooniverse_id'])
                entry = {'catalog_id':IDnumber, 'zooniverse_id':str(subject['zooniverse_id'])}
                
                #find location of FITS file; once non-FIRST sources are included, modify this and line 91
                fid = source['first_id']
                #if fid[0] == 'F':
                fits_loc = pathdict[fid]
                entry.update({'first_id':str(fid)})
                #else:
                #    raise RuntimeError('Not expecting non-FIRST data')
                #    fits_loc = '%s/rgz/raw_images/ATLAS/2x2/%s_radio.fits' % (data_path, fid)
                #    entry.update({'atlas_id':str(fid)})
                
                #find IR counterpart from consensus data, if present
                w = wcs.WCS(fits.getheader(fits_loc, 0)) #gets pixel-to-WCS conversion from header
                ir_coords = source['ir_peak']
                if ir_coords[0] == -99:
                    ir_pos = None
                    wise_match = None
                    sdss_match = None
                else:
                    #this only works for FIRST images; will need changing when ATLAS is added
                    p2w = w.wcs_pix2world
                    ir_ra_pixels = ir_coords[0] * 132./500.
                    ir_dec_pixels = 133 - ir_coords[1] * 132./500.
                    ir_peak = p2w( np.array([[ir_ra_pixels, ir_dec_pixels]]), 1)
                    ir_pos = coord.SkyCoord(ir_peak[0][0], ir_peak[0][1], unit=(u.deg,u.deg), frame='icrs')

                entry.update({'consensus':{'n_votes':source['n_votes'], 'n_total':source['n_total'], \
                                           'level':source['consensus_level'], 'label':source['label']}})
                if ir_pos:
                    logging.info('IR counterpart found')
                    entry['consensus'].update({'IR_ra':ir_pos.ra.deg, 'IR_dec':ir_pos.dec.deg})
                else:
                    logging.info('No IR counterpart found')

                #if an IR peak exists, search AllWISE and SDSS for counterparts
                if ir_pos:

                    wise_match = p.getWISE(entry)
                    if wise_match:
                        entry.update({'AllWISE':wise_match})
                    
                    sdss_match = p.getSDSS(entry)
                    if sdss_match:
                        entry.update({'SDSS':sdss_match})

                #try block attempts to read JSON from web; if it exists, calculate data
                try:
                    link = subject['location']['contours'] #gets url as Unicode string

                    # Use local file if available

                    jsonfile = link.split("/")[-1]
                    jsonfile_path = "{0}/rgz/contours/{1}".format(data_path,jsonfile)
                    if os.path.exists(jsonfile_path):
                        with open(jsonfile_path,'r') as jf:
                            data = json.load(jf)

                    # Otherwise, read from web

                    else:

                        # Reform weblink to point to the direct S3 URL, which will work even with older SSLv3
                        
                        link_s3 = "http://zooniverse-static.s3.amazonaws.com/"+link.split('http://')[-1]
                        
                        tryCount = 0
                        while(True): #in case of error, wait 10 sec and try again; give up after 5 tries
                            tryCount += 1
                            try:
                                compressed = urllib2.urlopen(str(link_s3)).read() #reads contents of url to str
                                break
                            except (urllib2.URLError, urllib2.HTTPError) as e:
                                if tryCount>5:
                                    message = 'Unable to connect to Amazon Web Services; trying again in 10 min'
                                    logging.exception(message)
                                    print message
                                    raise fn.DataAccessError(message)
                                logging.exception(e)
                                time.sleep(10)
                        
                        tempfile = StringIO.StringIO(compressed) #temporarily stores contents as file (emptied after unzipping)
                        uncompressed = gzip.GzipFile(fileobj=tempfile, mode='r').read() #unzips contents to str
                        data = json.loads(uncompressed) #loads JSON object
                    
                    radio_data = p.getRadio(data, fits_loc, source)
                    entry.update(radio_data)

                    #create RGZ name from radio position
                    radio_ra = radio_data['radio']['ra']
                    ra_h = int(radio_ra/15.)
                    ra_m = int((radio_ra - ra_h*15)*4)
                    ra_s = (radio_ra - ra_h*15 - ra_m/4.)*240
                    radio_dec = radio_data['radio']['dec']
                    dec_d = int(radio_dec)
                    dec_m = int((radio_dec - dec_d)*60)
                    dec_s = int((radio_dec - dec_d - dec_m/60.)*3600)
                    name = 'RGZJ{:0=2}{:0=2}{:0=4.1f}{:0=+3}{:0=2}{:0=2}'.format(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)
                    entry.update({'rgz_name':name})

                    #calculate physical data using redshift
                    if sdss_match and 'redshift' in sdss_match:
                        z = sdss_match['redshift']
                        lz = np.log10(z)
                        DAkpc = pow(10, -0.0799*pow(lz,3)-0.406*pow(lz,2)+0.3101*lz+3.2239)*1000 #angular size distance approximation in kpc
                        DLkpc = DAkpc*np.square(1+z) #luminosity distance approximation in kpc
                        maxPhysicalExtentKpc = DAkpc*radio_data['radio']['maxAngularExtent']*np.pi/180/3600 #arcseconds to radians
                        totalCrossSectionKpc2 = np.square(DAkpc)*radio_data['radio']['totalSolidAngle']*np.square(np.pi/180/3600) #arcseconds^2 to radians^2
                        totalLuminosityWHz = radio_data['radio']['totalFlux']*1e-29*4*np.pi*np.square(DLkpc*3.09e19) #mJy to W/(m^2 Hz), kpc to m
                        totalLuminosityErrWHz = radio_data['radio']['totalFluxErr']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                        peakLuminosityErrWHz = radio_data['radio']['peakFluxErr']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                        for component in radio_data['radio']['components']:
                            component['physicalExtent'] = DAkpc*component['angularExtent']*np.pi/180/3600
                            component['crossSection'] = np.square(DAkpc)*component['solidAngle']*np.square(np.pi/180/3600)
                            component['luminosity'] = component['flux']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                            component['luminosityErr'] = component['fluxErr']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                        for peak in radio_data['radio']['peaks']:
                            peak['luminosity'] = peak['flux']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                        entry['radio'].update({'maxPhysicalExtent':maxPhysicalExtentKpc, 'totalCrossSection':totalCrossSectionKpc2, \
                                               'totalLuminosity':totalLuminosityWHz, 'totalLuminosityErr':totalLuminosityErrWHz, \
                                               'peakLuminosityErr':peakLuminosityErrWHz})

                    logging.info('Radio data added')
                                       
                #if the link doesn't have a JSON, no data can be determined
                except urllib2.HTTPError as e:
                    if e.code == 404:
                        logging.info('No radio JSON detected')
                    else:
                        logging.exception(e)
                        raise
                
                catalog.insert(entry)
                find_duplicates(entry['zooniverse_id'])
                logging.info('Entry %i added to catalog', IDnumber)
        
        with open(in_progress_file, 'w') as f:
            f.write('')
        
    #end timer
    endtime = time.time()
    output = 'Time taken: %f' % (endtime-starttime)
    logging.info(output)
    print output

    return count

if __name__ == '__main__':
    logging.basicConfig(filename='RGZcatalog.log', level=logging.DEBUG, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logging.captureWarnings(True)
    logging.info('Catalog run from command line')
    done = False
    while not done:
        try:
            output = '%i entries added.' % RGZcatalog()
            logging.info(output)
            print output
            done = True
        except pymongo.errors.CursorNotFound as c:
            time.sleep(10)
            output = 'Cursor timed out; starting again.'
            logging.info(output)
            print output
        except fn.DataAccessError as d:
            resume = datetime.datetime.now() + datetime.timedelta(minutes=10)
            output = "RGZcatalog.py can't connect to external server; will resume at {:%H:%M}".format(resume)
            logging.info(output)
            print output
            time.sleep(600)
        except BaseException as e:
            logging.exception(e)
            raise
