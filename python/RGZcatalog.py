'''
RGZcatalog is a pipeline that takes all of the completed RGZ subjects and creates a Mongo database containing
consensus matching information, radio morphology, IR counterpart location, and data from corresponding
AllWISE and SDSS catalogs.
'''

import logging, urllib2, time, json, os, datetime
import pymongo
import numpy as np
import StringIO, gzip
from astropy.io import fits
from astropy import wcs, coordinates as coord, units as u

#custom modules for the RGZ catalog pipeline
import catalog_functions as fn #contains miscellaneous helper functions
import processing as p #contains functions that process the data
from find_duplicates import find_duplicates #finds and marks any radio components that are duplicated between sources

from consensus import rgz_path, data_path, db, version, logfile
in_progress_file = '%s/subject_in_progress.txt' % rgz_path

def RGZcatalog():
    
    #start timer
    starttime = time.time()
    
    #begin logging even if not run from command line
    logging.basicConfig(filename='{}/{}'.format(rgz_path,logfile), level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logging.captureWarnings(True)
    
    #connect to database of subjects
    subjects = db['radio_subjects']
    consensus = db['consensus{}'.format(version)]
    catalog = db['catalog{}'.format(version)] #this is being populated by this program
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
        to_be_completed = to_be_completed.union(in_progress_zid)
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

                entry.update({'consensus':{'n_radio':source['n_votes'], 'n_total':source['n_total'], 'n_ir':source['n_ir'], 'ir_flag':source['ir_flag'], \
                                           'ir_level':source['ir_level'], 'radio_level':source['consensus_level'], 'label':source['label']}})
                if ir_pos:
                    logging.info('IR counterpart found')
                    entry['consensus'].update({'ir_ra':ir_pos.ra.deg, 'ir_dec':ir_pos.dec.deg})
                else:
                    logging.info('No IR counterpart found')

                #if an IR peak exists, search AllWISE and SDSS for counterparts
                if ir_pos:

                    wise_match = p.getWISE(entry)
                    if wise_match:
                        entry.update({'AllWISE':wise_match})

                    tryCount = 0
                    while(True):
                        tryCount += 1
                        try:
                            sdss_match = p.getSDSS(entry)
                            if sdss_match:
                                entry.update({'SDSS':sdss_match})
                            break
                        except KeyError as e:
                            if tryCount>5:
                                message = 'Bad response from SkyServer; trying again in 10 min'
                                logging.exception(message)
                                print message
                                raise fn.DataAccessError(message)
                            elif e.message == 'ra':
                                #unable to reproduce; no error when I try again, so let's just do that
                                logging.exception(e)
                                time.sleep(10)
                            else:
                                raise e

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

                    #use WISE catalog name if available
                    if wise_match:
                        entry.update({'rgz_name':'RGZ{}{}'.format(wise_match['designation'][5:14], wise_match['designation'][15:22])})
                        
                    else:
                        #if not, try consensus IR position
                        if ir_pos:
                            ra = ir_pos.ra.deg
                            dec = ir_pos.dec.deg
                        #finally, just use radio center
                        else:
                            ra = radio_data['radio']['ra']
                            dec = radio_data['radio']['dec']
                        
                        ra_h = int(ra/15.)
                        ra_m = int((ra - ra_h*15)*4)
                        ra_s = (ra - ra_h*15 - ra_m/4.)*240
                        dec_d = int(dec)
                        dec_m = int((dec - dec_d)*60)
                        dec_s = int((dec - dec_d - dec_m/60.)*3600)
                        entry.update({'rgz_name':'RGZJ{:0=2}{:0=2}{:0=4.1f}{:0=+3}{:0=2}{:0=2}'.format(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)})
                    
                    #calculate physical data using redshift
                    if sdss_match:
                        z = 0
                        if 'spec_redshift' in sdss_match:
                            z = sdss_match['spec_redshift']
                        elif 'photo_redshift' in sdss_match:
                            z = sdss_match['photo_redshift']
                        if z>0:
                            lz = np.log10(z)
                            DAkpc = pow(10, -0.0799*pow(lz,3)-0.406*pow(lz,2)+0.3101*lz+3.2239)*1000 #angular size distance approximation in kpc
                            DLkpc = DAkpc*np.square(1+z) #luminosity distance approximation in kpc
                            maxPhysicalExtentKpc = DAkpc*radio_data['radio']['max_angular_extent']*np.pi/180/3600 #arcseconds to radians
                            totalCrossSectionKpc2 = np.square(DAkpc)*radio_data['radio']['total_solid_angle']*np.square(np.pi/180/3600) #arcseconds^2 to radians^2
                            totalLuminosityWHz = radio_data['radio']['total_flux']*1e-29*4*np.pi*np.square(DLkpc*3.09e19) #mJy to W/(m^2 Hz), kpc to m
                            totalLuminosityErrWHz = radio_data['radio']['total_flux_err']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                            peakLuminosityErrWHz = radio_data['radio']['peak_flux_err']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                            for component in radio_data['radio']['components']:
                                component['physical_extent'] = DAkpc*component['angular_extent']*np.pi/180/3600
                                component['cross_section'] = np.square(DAkpc)*component['solid_angle']*np.square(np.pi/180/3600)
                                component['luminosity'] = component['flux']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                                component['luminosity_err'] = component['flux_err']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                            for peak in radio_data['radio']['peaks']:
                                peak['luminosity'] = peak['flux']*1e-29*4*np.pi*np.square(DLkpc*3.09e19)
                            entry['radio'].update({'max_physical_extent':maxPhysicalExtentKpc, 'total_cross_section':totalCrossSectionKpc2, \
                                                   'total_luminosity':totalLuminosityWHz, 'total_luminosity_err':totalLuminosityErrWHz, \
                                                   'peak_luminosity_err':peakLuminosityErrWHz})

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
    logging.basicConfig(filename='{}/{}'.format(rgz_path,logfile), level=logging.DEBUG, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
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
