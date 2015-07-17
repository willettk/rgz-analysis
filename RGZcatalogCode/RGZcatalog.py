#to copy individual FITS files for test runs:
#find /data/extragal/willett/rgz/raw_images/ -name '[whatever].fits' | xargs cp -t .

from pymongo import MongoClient
import numpy as np
import urllib2
from StringIO import StringIO
from gzip import GzipFile
import json
from ast import literal_eval
import os
from astropy.io import fits
from astropy import wcs
from astroquery.irsa import Irsa
from astroquery.sdss import SDSS
import astropy.coordinates as coord
import astropy.units as u
import catalogFunctions as fn #contains custom functions
import contourNode as c #contains Node class

def RGZcatalog():

    #connect to database of subjects
    client = MongoClient("localhost",27017)
    db = client['radio']
    subjects = db['radio_subjects']
    #import consensus database with
    #mongoimport -d radio -c consensus --file consensus_all_75percent.csv --type csv --headerline
    consensus = db['consensus']

    fits_dir = raw_input('Directory where radio FITS images are stored: ')
    while not os.path.exists(fits_dir):
        fits_dir = raw_input('Invalid path. Please re-enter: ')

    #initilize file for output
    with open('RGZcatalog.json', 'w') as f:
##        #output = 'zid,totalFlux,totalFluxUncertainty,outerLevel,numberComponents,numberPeaks,maxExtent,totalArea,peakFluxDensityUncertainty,peaks\n' #csv
##        output = 'catalog id,zooniverse id,FIRST id,n_users,n_total,consensus level,IR counterpart ra,IR counterpart dec,total flux (mJy),' + \
##                 'total flux uncertainty (mJy),outermost level (mJy/beam),number of components,number of peaks,max extent of bboxes (arcminutes),' + \
##                 'total area of components (arcminutes^2),peak flux density uncertainty (mJy/beam),peaks (flux density in mJy/beam; ra and dec in degrees),,,' + \
##                 'components (flux and uncertainty in mJy; extent in arcminutes; area in arcminutes^2; ra and dec in degrees),AllWISE match\n' #human readable
##        f.write(output)
        f.write('{')

        IDnumber = 0 #ID number for catalog
        count = 0 #number of sources added

        #iterate through all subjects
        #for subject in subjects.find().batch_size(30):
        for subject in subjects.find({'zooniverse_id': {'$in': ['ARG00000sl', 'ARG0003f9l']} }):
        #for subject in subjects.find({'zooniverse_id':'ARG0002qyh'}): #sample multipeaked subject
        #for subject in subjects.find({'zooniverse_id':'ARG00000sl'}): #sample subject with distinct galaxies
        #for subject in subjects.find({'zooniverse_id':'ARG0003f9l'}): #sample subject with multiple components

            IDnumber += 1

            link = subject['location']['contours'] #gets url as Unicode string

            for consensusObject in consensus.find({'zooniverse_id':subject['zooniverse_id']}):

                count += 1

                #find location of FITS file
                fid = consensusObject['FIRST_id']
                for root, dirnames, filenames in os.walk(fits_dir):
                    for filename in filenames:
                        if filename == fid + '.fits':
                            fits_loc = os.path.join(root, filename)

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
                ir_ra = ir_peak[0][0]
                ir_dec = ir_peak[0][1]
                pos = coord.SkyCoord(ir_ra, ir_dec, unit=(u.deg,u.deg), frame='icrs')

                #get IR data from AllWISE Source Catalog
                table = Irsa.query_region(pos, catalog='wise_allwise_p3as_psd', radius=3*u.arcsec)
                if len(table):
                    match = table[0]
                    if len(table)>1:
                        dist = match['dist']
                        for entry in table:
                            if entry['dist'] < dist:
                                match = entry
                                dist = match['dist']
                    wise_match = {'AllWISE_id':'WISEA'+match['designation'], 'ra':match['ra'], 'dec':match['dec'], 'numberMatches':len(table), \
                                  'w1mpro':match['w1mpro'], 'w1sigmpro':match['w1sigmpro'], 'w1snr':match['w1snr'], \
                                  'w2mpro':match['w2mpro'], 'w2sigmpro':match['w2sigmpro'], 'w2snr':match['w2snr'], \
                                  'w3mpro':match['w3mpro'], 'w3sigmpro':match['w3sigmpro'], 'w3snr':match['w3snr'], \
                                  'w4mpro':match['w4mpro'], 'w4sigmpro':match['w4sigmpro'], 'w4snr':match['w4snr']}
                    for key in wise_match:
                        if wise_match[key] is np.ma.masked:
                            wise_match[key] = None
                else:
                    wise_match = None

                #get optical data from SDSS; requires numpy 0.10
                #table = SDSS.query_region(pos, radius=3*u.arcsec)

                #save consensus data as dict for printing to JSON
                catalog_id = 'RGZ'+str(IDnumber)+str(consensusObject['label'])
##                consensusOutputStr = catalog_id + ',' + str(subject['zooniverse_id']) + ',' + \
##                                   str(fid) + ',' + str(consensusObject['n_users']) + ',' + str(consensusObject['n_total']) + ',' + \
##                                   str(consensusObject['consensus_level'])
                outputDict = { 'id':catalog_id, 'Zooniverse_id':str(subject['zooniverse_id']), 'FIRST_id':str(fid), 'AllWISE_match':wise_match, \
                               'consensus':{'n_users':consensusObject['n_users'], 'n_total':consensusObject['n_total'], \
                                            'consensus_level':consensusObject['consensus_level'], 'IR_ra':ir_ra, 'IR_dec':ir_dec} }

                #try block attempts to read JSON from web; if it exists, calculate data
                try:
                    compressed = urllib2.urlopen(str(link)).read() #reads contents of url to str
                    tempfile = StringIO(compressed) #temporarily stores contents as file (emptied after unzipping)
                    uncompressed = GzipFile(fileobj=tempfile, mode='r').read() #unzips contents to str
                    data = json.loads(uncompressed) #loads JSON object

                    #scaling factors, pixels to arcminutes
                    w = 3./data['width']
                    h = 3./data['height']

                    #create list of trees, each containing a contour and its contents
                    contourTrees = []
                    consensusBboxes = literal_eval(consensusObject['bbox'])
                    for contour in data['contours']:
                        for bbox in consensusBboxes:
                            if fn.approx(contour[0]['bbox'][0], bbox[0]) and fn.approx(contour[0]['bbox'][1], bbox[1]) and \
                               fn.approx(contour[0]['bbox'][2], bbox[2]) and fn.approx(contour[0]['bbox'][3], bbox[3]):
                                tree = c.Node(contour=contour, fits_loc=fits_loc)
                                contourTrees.append(tree)

                    #get component fluxes and sizes
                    components = []
                    for tree in contourTrees:
                        bboxP = fn.bboxToDS9(fn.findBox(tree.value['arr']))[0] #bbox in DS9 coordinate pixels
                        bboxCornersRD = tree.w.wcs_pix2world( np.array( [[bboxP[0],bboxP[1]], [bboxP[2],bboxP[3]] ]), 1) #two opposite corners of bbox in ra and dec
                        raRange = [ min(bboxCornersRD[0][0], bboxCornersRD[1][0]), max(bboxCornersRD[0][0], bboxCornersRD[1][0]) ]
                        decRange = [ min(bboxCornersRD[0][1], bboxCornersRD[1][1]), max(bboxCornersRD[0][1], bboxCornersRD[1][1]) ]
                        extent = fn.distance( [raRange[0],decRange[0]], [raRange[1],decRange[1]] )*60 #arcminutes
                        area = fn.area(tree.value['arr'])*w*h*3600 #square arcminutes
                        components.append({'flux':tree.flux, 'fluxUncertainty':tree.fluxUncertainty, 'extent':extent, 'area':area, \
                                           'raRange':raRange, 'decRange':decRange})

                    #adds up total flux of all components
                    totalFlux = 0
                    totalFluxUncertainty = 0
                    for component in components:
                        totalFlux += component['flux']
                        totalFluxUncertainty += pow(component['fluxUncertainty'], 2)
                    totalFluxUncertainty = np.sqrt(totalFluxUncertainty)

                    #finds total area enclosed by contours in arcminutes
                    totalArea = 0
                    for component in components:
                        totalArea += component['area']

                    #find maximum extent of component bboxes in degrees
                    maxExtent = 0
                    if len(components)==1:
                        maxExtent = components[0]['extent']
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
                                        extent = fn.distance(corner1, corner2)
                                        if extent>maxExtent:
                                            maxExtent = extent

                    #add all peaks up into single list
                    peakList = []
                    for tree in contourTrees:
                        for peak in tree.peaks:
                            peakList.append(peak)

                    #save radio data as dict for printing
##                    radioOutputStr = ',' + str(ir_ra) + ',' + str(ir_dec)+ ',' + str(totalFlux)+ ',' + \
##                             str(totalFluxUncertainty) + ',' + str(data['contours'][0][0]['level']*1000) + ',' + \
##                             str(len(contourTrees)) + ',' + str(len(peakList))+ ',' + str(maxExtent*60) + ',' + \
##                             str(totalArea) + ',' + str(contourTrees[0].sigma*1000) + ',' + str(peakList) + ',' + \
##                             str(components) + ',' + str(wise_match)
                    outputDict.update({'totalFlux':totalFlux, 'totalFluxUncertainty':totalFluxUncertainty, \
                                       'outermostLevel':data['contours'][0][0]['level']*1000, 'numberComponents':len(contourTrees), \
                                       'numberPeaks':len(peakList), 'maxExtent':maxExtent, 'totalArea':totalArea, \
                                       'peakFluxUncertainty':contourTrees[0].sigma*1000, 'peaks':peakList, 'components':components})
                                       
                #if the link doesn't have a JSON, no data can be determined
                except urllib2.HTTPError, err:
                    if err.code == 404:
##                        radioOutputStr = ',None,None,None,None,None,None,None,None,None,None,None,None\n'
                        outputDict.update({'totalFlux':None, 'totalFluxUncertainty':None, 'outermostLevel':None, \
                                           'numberComponents':None, 'numberPeaks':None, 'maxExtent':None, 'totalArea':None, \
                                           'peakFluxUncertainty':None, 'peaks':None, 'components':None})
                        pass
                    else:
                        raise

                #output data to file one entry at a time
##                outputStr = consensusOutputStr + radioOutputStr + '\n'
##                f.write(outputStr)
                if count>1:
                    f.write(', ')
                f.write('"' + catalog_id + '": ' + json.dumps(outputDict))

        f.write('}') #close the JSON

    return [count, IDnumber]

if __name__ == '__main__':
    counts = RGZcatalog()
    print counts[0], 'entries added from', counts[1], 'subject fields'
