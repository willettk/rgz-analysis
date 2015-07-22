#to copy individual FITS files for test runs:
#find /data/extragal/willett/rgz/raw_images/ -name '[whatever].fits' | xargs cp -t .

from pymongo import MongoClient
import numpy as np
import pandas as pd
import urllib2
from StringIO import StringIO
from gzip import GzipFile
import json
from ast import literal_eval
import os
import scandir
from astropy.io import fits
from astropy import wcs
from astropy import coordinates as coord
from astropy import units as u
from astroquery.irsa import Irsa
import mechanize
import time
import contourNode as c #contains Node class

def RGZcatalog():

    #connect to database of subjects
    client = MongoClient('localhost',27017)
    db = client['radio']
    subjects = db['radio_subjects']
    #import consensus database with
    #mongoimport -d radio -c consensus --file consensus_all_75percent.csv --type csv --headerline
    consensus = db['consensus']
    catalog = db['catalog'] #this is being populated by this program

    #for testing
    if catalog.count():
        overwrite = raw_input('Catalog has entries. Overwrite or append? (o/a) ').lower()
        while overwrite!='o' and overwrite!='a':
            overwrite = raw_input('Invalid choice. Overwrite or append? (o/a) ').lower()
        if overwrite=='o':
            db.drop_collection('catalog')

    fits_dir = raw_input('Directory where radio FITS images are stored: ')
    while not os.path.exists(fits_dir):
        fits_dir = raw_input('Invalid path. Please re-enter: ')

    IDnumber = 0 #ID number for catalog
    count = 0 #number of sources added

    #start timer
    starttime = time.time()

    #iterate through all subjects
    for subject in subjects.find().batch_size(30):
    #for subject in subjects.find({'zooniverse_id': {'$in': ['ARG00000sl', 'ARG0003f9l']} }):
    #for subject in subjects.find({'zooniverse_id':'ARG0002qyh'}): #sample multipeaked subject
    #for subject in subjects.find({'zooniverse_id':'ARG00000sl'}): #sample subject with distinct galaxies
    #for subject in subjects.find({'zooniverse_id':'ARG0003f9l'}): #sample subject with multiple components

        IDnumber += 1

        link = subject['location']['contours'] #gets url as Unicode string

        for consensusObject in consensus.find({'zooniverse_id':subject['zooniverse_id']}):

            count += 1
            if count>10:
                print 'Time taken:', time.time()-starttime
                return [count, IDnumber]

            #display which entry is being processed to see how far the program is
            catalog_id = 'RGZ'+str(IDnumber)+str(consensusObject['label'])
            print catalog_id
            
            searchtime = time.time()
            #find location of FITS file
            fid = consensusObject['FIRST_id']
            for root, dirnames, filenames in scandir.walk(fits_dir):
                for filename in filenames:
                    if filename == fid + '.fits':
                        fits_loc = os.path.join(root, filename)
            print 'Time spent finding FITS:', time.time()-searchtime
            
            #find IR counterpart from consensus data, if present
            w = wcs.WCS(fits.open(fits_loc)[0].header) #gets pixel-to-WCS conversion from header
            ir_coords = literal_eval(consensusObject['ir_peak'])
            if ir_coords == (-99, -99):
                ir_peak = np.array([[None, None]])
            else:
                p2w = w.wcs_pix2world
                ir_ra_pixels = ir_coords[0] * 132./500.
                ir_dec_pixels = 133 - ir_coords[1] * 132./500.
                ir_peak = p2w( np.array([[ir_ra_pixels, ir_dec_pixels]]), 1)
            ir_pos = coord.SkyCoord(ir_peak[0][0], ir_peak[0][1], unit=(u.deg,u.deg), frame='icrs')

            #if an IR peak exists, search AllWISE and SDSS for counterparts
            if ir_pos:
                wisetime = time.time()
                #get IR data from AllWISE Source Catalog
                table = Irsa.query_region(ir_pos, catalog='wise_allwise_p3as_psd', radius=3*u.arcsec)
                if len(table):
                    numberMatches = 0
                    if table[0]['w1snr']>5:
                        match = table[0]
                        dist = match['dist']
                        numberMatches += 1
                    else:
                        match = None
                        dist = np.inf
                    if len(table)>1:
                        for entry in table:
                            if entry['dist']<dist and entry['w1snr']>5:
                                match = entry
                                dist = match['dist']
                                numberMatches += 1
                    if match:
                        wise_match = {'designation':'WISEA'+match['designation'], 'ra':match['ra'], 'dec':match['dec'], 'numberMatches':np.int16(numberMatches), \
                                      'w1mpro':match['w1mpro'], 'w1sigmpro':match['w1sigmpro'], 'w1snr':match['w1snr'], \
                                      'w2mpro':match['w2mpro'], 'w2sigmpro':match['w2sigmpro'], 'w2snr':match['w2snr'], \
                                      'w3mpro':match['w3mpro'], 'w3sigmpro':match['w3sigmpro'], 'w3snr':match['w3snr'], \
                                      'w4mpro':match['w4mpro'], 'w4sigmpro':match['w4sigmpro'], 'w4snr':match['w4snr']}
                        for key in wise_match:
                            if wise_match[key] is np.ma.masked:
                                wise_match[key] = None
                    else:
                        wise_match = None
                else:
                    wise_match = None

                if wise_match:
                    matched = 'match'
                else:
                    matched = 'no match'
                print 'WISE:', time.time()-wisetime, matched
                sdsstime = time.time()

                #get optical magnitude data from Galaxy table in SDSS
                query = '''select objID, ra, dec, u, g, r, i, z, err_u, err_g, err_r, err_i, err_z from Galaxy
                           where (ra between ''' + str(ir_pos.ra.deg) + '-1.5/3600 and ' + str(ir_pos.ra.deg) + '''+1.5/3600) and
                                 (dec between ''' + str(ir_pos.dec.deg) + '-1.5/3600 and ' + str(ir_pos.dec.deg) + '+1.5/3600)'
                df = SDSS_select(query)
                if len(df):
                    numberMatches = 0
                    matchPos = coord.SkyCoord(df.iloc[0]['ra'], df.iloc[0]['dec'], unit=(u.deg, u.deg))
                    tempDist = ir_pos.separation(matchPos).arcsecond
                    if tempDist<3:
                        match = df.iloc[0]
                        dist = tempDist
                        numberMatches += 1
                    else:
                        match = None
                        dist = np.inf
                    if len(df)>1:
                        for i in range(len(df)):
                            matchPos = coord.SkyCoord(df.iloc[i]['ra'], df.iloc[i]['dec'], unit=(u.deg, u.deg))
                            tempDist = ir_pos.separation(matchPos).arcsecond
                            if tempDist<3 and tempDist<dist:
                                match = df.iloc[i]
                                dist = tempDist
                                numberMatches += 1
                    if match is not None:
                        sdss_match = {'objID':df['objID'][match.name], 'ra':match['ra'], 'dec':match['dec'], 'numberMatches':np.int16(numberMatches), \
                                      'u':match['u'], 'g':match['g'], 'r':match['r'], 'i':match['i'], 'z':match['z'], \
                                      'u_err':match['err_u'], 'g_err':match['err_g'], 'r_err':match['err_r'], 'i_err':match['err_i'], 'z_err':match['err_z']}
                    else:
                        sdss_match = None
                else:
                    sdss_match = None

                #only get more data from SDSS if a postional match exists
                if sdss_match:

                    #get photo redshift and uncertainty from Photoz table
                    query = 'select z, zErr from Photoz where objID=' + str(sdss_match['objID'])
                    df = SDSS_select(query)
                    if len(df):
                        photoZ = df['z'][0]
                        photoZErr = df['zErr'][0]
                    else:
                        photoZ = None
                        photoZErr = None

                    #get spectral lines from GalSpecLine tables
                    query = '''select oiii_5007_flux, oiii_5007_flux_err, h_beta_flux, h_beta_flux_err,
                                      nii_6584_flux, nii_6584_flux_err, h_alpha_flux, h_alpha_flux_err
                               from GalSpecLine AS g
                                  join SpecObj AS s ON s.specobjid = g.specobjid
                               where s.bestObjID = ''' + str(sdss_match['objID'])
                    df = SDSS_select(query)
                    if len(df):
                        sdss_match.update({'oiii_5007_flux':df['oiii_5007_flux'][0], 'oiii_5007_flux_err':df['oiii_5007_flux_err'][0], \
                                           'h_beta_flux':df['h_beta_flux'][0], 'h_beta_flux_err':df['h_beta_flux_err'][0], \
                                           'nii_6584_flux':df['nii_6584_flux'][0], 'nii_6584_flux_err':df['nii_6584_flux_err'][0], \
                                           'h_alpha_flux':df['h_alpha_flux'][0], 'h_alpha_flux_err':df['h_alpha_flux_err'][0]})
                    else:
                        sdss_match.update({'oiii_5007_flux':None, 'oiii_5007_flux_err':None, 'h_beta_flux':None, 'h_beta_flux_err':None, \
                                           'nii_6584_flux':None, 'nii_6584_flux_err':None, 'h_alpha_flux':None, 'h_alpha_flux_err':None})

                    #get spectral class and redshiftfrom SpecPhoto table
                    query = '''select z, zErr, case when class like 'GALAXY' then 0
                                                    when class like 'QSO' then 1
                                                    when class like 'STAR' then 2 end as classNum
                               from SpecObj
                               where bestObjID = ''' + str(sdss_match['objID'])
                    df = SDSS_select(query)
                    if len(df):
                        sdss_match.update({'spectralClass':np.int16(df['classNum'][0])})
                        specZ = df['z'][0]
                        specZErr = df['zErr'][0]
                    else:
                        sdss_match.update({'spectralClass':None})
                        specZ = None
                        specZErr = None

                    #use specZ is present, otherwise use photoZ
                    if specZ:
                        sdss_match.update({'redshift':specZ, 'redshift_err':specZErr, 'redshift_type':np.int16(1)})
                    else:
                        sdss_match.update({'redshift':photoZ, 'redshift_err':photoZErr, 'redshift_type':np.int16(0)})

                #end of 'if SDSS match'

            #end of 'if IR peak'
            if sdss_match:
                matched = 'match'
            else:
                matched = 'no match'
            print 'SDSS:', time.time()-sdsstime, matched

            #convert match data from numpy types to native Python types for JSON compatibility
            if wise_match:
                for key in wise_match:
                    if wise_match[key] and type(wise_match[key]) is not str:
                        wise_match[key] = wise_match[key].item()
                    elif wise_match[key] == 0:
                        wise_match[key] = 0
            if sdss_match:
                for key in sdss_match:
                    if sdss_match[key] and type(sdss_match[key]) is not str:
                        sdss_match[key] = sdss_match[key].item()
                    elif sdss_match[key] == 0:
                        sdss_match[key] = 0

            #save consensus data as dict for printing to JSON
            outputDict = { 'catalog_id':catalog_id, 'Zooniverse_id':str(subject['zooniverse_id']), 'FIRST_id':str(fid), 'AllWISE':wise_match, 'SDSS':sdss_match, \
                           'consensus':{'n_users':consensusObject['n_users'], 'n_total':consensusObject['n_total'], \
                                        'level':consensusObject['consensus_level'], 'IR_ra':ir_pos.ra.deg, 'IR_dec':ir_pos.dec.deg} }

            radiotime = time.time()
            #try block attempts to read JSON from web; if it exists, calculate data
            try:
                compressed = urllib2.urlopen(str(link)).read() #reads contents of url to str
                tempfile = StringIO(compressed) #temporarily stores contents as file (emptied after unzipping)
                uncompressed = GzipFile(fileobj=tempfile, mode='r').read() #unzips contents to str
                data = json.loads(uncompressed) #loads JSON object

                #scaling factors, pixels to arcminutes
                #w = 3./data['width']
                #h = 3./data['height']

                #create list of trees, each containing a contour and its contents
                contourTrees = []
                consensusBboxes = literal_eval(consensusObject['bbox'])
                for contour in data['contours']:
                    for bbox in consensusBboxes:
                        if approx(contour[0]['bbox'][0], bbox[0]) and approx(contour[0]['bbox'][1], bbox[1]) and \
                           approx(contour[0]['bbox'][2], bbox[2]) and approx(contour[0]['bbox'][3], bbox[3]):
                            tree = c.Node(contour=contour, fits_loc=fits_loc)
                            contourTrees.append(tree)

                #get component fluxes and sizes
                components = []
                for tree in contourTrees:
                    bboxP = bboxToDS9(findBox(tree.value['arr']))[0] #bbox in DS9 coordinate pixels
                    bboxCornersRD = tree.w.wcs_pix2world( np.array( [[bboxP[0],bboxP[1]], [bboxP[2],bboxP[3]] ]), 1) #two opposite corners of bbox in ra and dec
                    raRange = [ min(bboxCornersRD[0][0], bboxCornersRD[1][0]), max(bboxCornersRD[0][0], bboxCornersRD[1][0]) ]
                    decRange = [ min(bboxCornersRD[0][1], bboxCornersRD[1][1]), max(bboxCornersRD[0][1], bboxCornersRD[1][1]) ]
                    pos1 = coord.SkyCoord(raRange[0], decRange[0], unit=(u.deg, u.deg))
                    pos2 = coord.SkyCoord(raRange[1], decRange[1], unit=(u.deg, u.deg))
                    extent = pos1.separation(pos2).arcminute
                    solidAngle = tree.area #square arcsec
                    components.append({'flux':tree.flux, 'fluxErr':tree.fluxErr, 'angularExtent':extent, 'solidAngle':solidAngle, \
                                       'raRange':raRange, 'decRange':decRange, 'physicalExtent':None, 'crossSection':None})

                #adds up total flux of all components
                totalFlux = 0
                totalFluxErr = 0
                for component in components:
                    totalFlux += component['flux']
                    totalFluxErr += pow(component['fluxErr'], 2)
                totalFluxErr = np.sqrt(totalFluxErr)

                #finds total area enclosed by contours in arcminutes
                totalSolidAngle = 0
                for component in components:
                    totalSolidAngle += component['solidAngle']

                #find maximum extent of component bboxes in arcseconds
                maxAngularExtent = 0
                if len(components)==1:
                    maxAngularExtent = components[0]['angularExtent']
                else:
                    for i in range(len(components)-1):
                        for j in range(1,len(components)-i):
                            corners1 = [ [components[i]['raRange'][0],components[i]['decRange'][0]], \
                                         [components[i]['raRange'][0],components[i]['decRange'][1]], \
                                         [components[i]['raRange'][1],components[i]['decRange'][0]], \
                                         [components[i]['raRange'][1],components[i]['decRange'][1]] ]
                            corners2 = [ [components[i+j]['raRange'][0],components[i+j]['decRange'][0]], \
                                         [components[i+j]['raRange'][0],components[i+j]['decRange'][1]], \
                                         [components[i+j]['raRange'][1],components[i+j]['decRange'][0]], \
                                         [components[i+j]['raRange'][1],components[i+j]['decRange'][1]] ]
                            for corner1 in corners1:
                                for corner2 in corners2:
                                    pos1 = coord.SkyCoord(corner1[0], corner1[1], unit=(u.deg, u.deg))
                                    pos2 = coord.SkyCoord(corner2[0], corner2[1], unit=(u.deg, u.deg))
                                    angularExtent = pos1.separation(pos2).arcminute
                                    if angularExtent>maxAngularExtent:
                                        maxAngularExtent = angularExtent

                #add all peaks up into single list
                peakList = []
                for tree in contourTrees:
                    for peak in tree.peaks:
                        peakList.append(peak)

                #calculate physical data using redshift
                if sdss_match:
                    if sdss_match['redshift']:
                        z = sdss_match['redshift']
                        lz = np.log10(z)
                        D_A = pow(10, -0.0799*pow(lz,3)-0.406*pow(lz,2)+0.3101*lz+3.2239) #angular size distance approximation
                        D_L = D_A*pow(1+z, 2) #luminosity distance approximation
                        maxPhysicalExtent = D_A*maxAngularExtent*np.pi/180/3600 #arcseconds to radians
                        totalCrossSection = pow(D_A,2)*totalSolidAngle*pow(np.pi/180/60,2) #arcminutes^2 to radians^2
                        for component in components:
                            component['physicalExtent'] = D_A*component['angularExtent']*np.pi/180/3600
                            component['crossSection'] = pow(D_A,2)*component['solidAngle']*pow(np.pi/180/60,2)
                else:
                    maxPhysicalExtent = None
                    totalCrossSection = None

                #save radio data as dict for printing
                outputDict.update({'totalFlux':totalFlux, 'totalFluxErr':totalFluxErr, \
                                   'outermostLevel':data['contours'][0][0]['level']*1000, 'numberComponents':len(contourTrees), 'numberPeaks':len(peakList), \
                                   'maxAngularExtent':maxAngularExtent, 'maxPhysicalExtent':maxPhysicalExtent, 'totalSolidAngle':totalSolidAngle, \
                                   'totalCrossSection':totalCrossSection, 'peakFluxErr':contourTrees[0].sigma*1000, 'peaks':peakList, 'components':components})
                                   
            #if the link doesn't have a JSON, no data can be determined
            except urllib2.HTTPError, err:
                if err.code == 404:
                    outputDict.update({'totalFlux':None, 'totalFluxErr':None, 'outermostLevel':None, \
                                       'numberComponents':None, 'numberPeaks':None, 'maxExtent':None, 'totalArea':None, \
                                       'peakFluxErr':None, 'peaks':None, 'components':None})
                    pass
                else:
                    raise

            print 'radio:', time.time()-radiotime

            catalog.insert(outputDict)

    #end timer
    endtime = time.time()
    print 'Time taken:', endtime-starttime

    return [count, IDnumber]

#pass an SQL query to SDSS and return a pandas dataframe
def SDSS_select(sql):
    br = mechanize.Browser()
    br.open('http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx')
    br.select_form(name='sql')
    br['cmd'] = sql
    br['format'] = ['csv']
    response = br.submit()
    file_like = StringIO(response.get_data())
    return pd.read_csv(file_like, skiprows=1)

#determines if two floats are approximately equal
def approx(a, b, uncertainty=1e-5):
   return np.abs(a-b) < uncertainty

#creates a bounding box for a given contour path
#loop = data['contours'][0][0]['arr'] #outermost contour (for testing)
def findBox(loop):
   xmax = loop[0]['x']
   ymax = loop[0]['y']
   xmin = loop[0]['x']
   ymin = loop[0]['y']
   for i in loop:
      if i['x']>xmax:
         xmax = i['x']
      elif i['x']<xmin:
         xmin = i['x']
      if i['y']>ymax:
         ymax = i['y']
      elif i['y']<ymin:
         ymin = i['y']
   return [xmax, ymax, xmin, ymin]

#finds the coordinates of the bbox in DS9's system and the input values for drawing a box in DS9
#bbox = tree.value['bbox'] #outermost bbox (for testing)
def bboxToDS9(bbox):
    xmax = bbox[0]
    ymax = bbox[1]
    xmin = bbox[2]
    ymin = bbox[3]
    temp = 133-ymax
    ymax = 133-ymin
    ymin = temp
    newBbox = [xmax, ymax, xmin, ymin]
    ds9Box = [ (xmax+xmin)/2., (ymax+ymin)/2., xmax-xmin, ymax-ymin ]
    return [newBbox, ds9Box]

if __name__ == '__main__':
    counts = RGZcatalog()
    print counts[0], 'entries added from', counts[1], 'subject fields'
