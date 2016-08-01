import logging, time
from astropy import coordinates as coord, units as u
import mechanize, httplib, StringIO
from astroquery.exceptions import TimeoutError, TableParseError
from astroquery.irsa import Irsa
import numpy as np
import pandas as pd

#custom modules for the RGZ catalog
import catalog_functions as fn #contains miscellaneous helper functions
import contour_node as c #contains Node class

def getWISE(entry):
    '''
    get IR data from AllWISE Source Catalog
    attempts to query Irsa 5 times; if they keep failing, abort
    returns updated entry
    '''
    
    ir_pos = coord.SkyCoord(entry['consensus']['ir_ra'], entry['consensus']['ir_dec'], unit=(u.deg,u.deg), frame='icrs')
    
    tryCount = 0
    while(True): #in case of error, wait 10 sec and try again; give up after 5 tries
        tryCount += 1
        try:
            table = Irsa.query_region(ir_pos, catalog='wise_allwise_p3as_psd', radius=3.*u.arcsec)
            break
        except (TimeoutError, TableParseError) as e:
            if tryCount>5:
                message = 'Unable to connect to IRSA; trying again in 10 min'
                logging.exception(message)
                print message
                raise fn.DataAccessError(message)
            logging.exception(e)
            time.sleep(10)
        except Exception as e:
            if str(e) == 'Query failed\n':
                if tryCount>5:
                    message = 'Unable to connect to IRSA; trying again in 10 min'
                    logging.exception(message)
                    print message
                    raise fn.DataAccessError(message)
                logging.exception(e)
                time.sleep(10)
            else:
                raise
    
    if len(table):
        number_matches = 0
        if table[0]['w1snr']>5:
            match = table[0]
            dist = match['dist']
            number_matches += 1
        else:
            match = None
            dist = np.inf
        if len(table)>1:
            for row in table:
                if row['dist']<dist and row['w1snr']>5:
                    match = row
                    dist = match['dist']
                    number_matches += 1
        if match:
            wise_match = {'designation':'WISEA'+match['designation'], 'ra':match['ra'], 'dec':match['dec'], 'number_matches':np.int16(number_matches), \
                          'w1mpro':match['w1mpro'], 'w1sigmpro':match['w1sigmpro'], 'w1snr':match['w1snr'], \
                          'w2mpro':match['w2mpro'], 'w2sigmpro':match['w2sigmpro'], 'w2snr':match['w2snr'], \
                          'w3mpro':match['w3mpro'], 'w3sigmpro':match['w3sigmpro'], 'w3snr':match['w3snr'], \
                          'w4mpro':match['w4mpro'], 'w4sigmpro':match['w4sigmpro'], 'w4snr':match['w4snr']}     
        else:
            wise_match = None
    else:
        wise_match = None
    
    if wise_match:
        logging.info('AllWISE match found')
        for key in wise_match.keys():
            if wise_match[key] is np.ma.masked:
                wise_match.pop(key)
            elif wise_match[key] and type(wise_match[key]) is not str:
                wise_match[key] = wise_match[key].item()
            elif wise_match[key] == 0:
                wise_match[key] = 0
    else:
        logging.info('No AllWISE match found')
    
    return wise_match

def SDSS_select(sql):
   '''pass an SQL query to SDSS and return a pandas dataframe
   in case of error, wait 10 seconds and try again; give up after 5 tries'''
   logging.basicConfig(filename='RGZcatalog.log', level=logging.DEBUG, format='%(asctime)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
   br = mechanize.Browser()
   br.set_handle_robots(False)
   tryCount = 0
   while(True):
      tryCount += 1
      try:
         br.open('http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx', timeout=4)
         br.select_form(name='sql')
         br['cmd'] = sql
         br['format'] = ['csv']
         response = br.submit()
         file_like = StringIO.StringIO(response.get_data())
         df = pd.read_csv(file_like, skiprows=1)
         break
      except (mechanize.URLError, mechanize.HTTPError, httplib.BadStatusLine, pd.parser.CParserError) as e:
         if tryCount>5:
            message = 'Unable to connect to SkyServer; trying again in 10 min'
            logging.exception(message)
            print message
            raise fn.DataAccessError(message)
         logging.exception(e)
         time.sleep(10)
   return df

def getSDSS(entry):
    '''
    get optical magnitude data from Galaxy table in SDSS
    if a positional match exists, also get photo redshift and uncertainty from Photoz table, spectral lines from GalSpecLine table, and
    spectral class and spec redshift and uncertainty from SpecPhoto table
    '''
    
    ir_pos = coord.SkyCoord(entry['consensus']['ir_ra'], entry['consensus']['ir_dec'], unit=(u.deg,u.deg), frame='icrs')
    
    query = '''select objID, ra, dec, u, r, g, i, z, err_u, err_r, err_g, err_i, err_z,
                 case type when 3 then 'G'
                           when 6 then 'S'
                           else 'U' end as class
               from PhotoPrimary
               where (ra between %f-3./3600 and %f+3./3600) and (dec between %f-3./3600 and %f+3./3600)''' \
               % (ir_pos.ra.deg, ir_pos.ra.deg, ir_pos.dec.deg, ir_pos.dec.deg)
    df = SDSS_select(query)
    if len(df):
        number_matches = 0
        match_pos = coord.SkyCoord(df.iloc[0]['ra'], df.iloc[0]['dec'], unit=(u.deg, u.deg))
        temp_dist = ir_pos.separation(match_pos).arcsecond
        if temp_dist<3.:
            match = df.iloc[0]
            dist = temp_dist
            number_matches += 1
        else:
            match = None
            dist = np.inf
        if len(df)>1:
            for i in range(len(df)):
                match_pos = coord.SkyCoord(df.iloc[i]['ra'], df.iloc[i]['dec'], unit=(u.deg, u.deg))
                temp_dist = ir_pos.separation(match_pos).arcsecond
                if temp_dist<3. and temp_dist<dist:
                    match = df.iloc[i]
                    dist = temp_dist
                    number_matches += 1
        if match is not None:
            sdss_match = {'objID':df['objID'][match.name], 'ra':match['ra'], 'dec':match['dec'], 'number_matches':np.int16(number_matches), \
                          'morphological_class':match['class'], 'u':match['u'], 'r':match['r'], 'g':match['g'], 'i':match['i'], 'z':match['z'], \
                          'u_err':match['err_u'], 'r_err':match['err_r'], 'g_err':match['err_g'], 'i_err':match['err_i'], 'z_err':match['err_z']}
        else:
            sdss_match = None
    else:
        sdss_match = None

    if sdss_match and sdss_match['morphological_class'] == 'G': #query the galaxy tables

        query = '''select p.z as photo_redshift, p.zErr as photo_redshift_err, s.z as spec_redshift, s.zErr as spec_redshift_err,
                     oiii_5007_flux, oiii_5007_flux_err, h_beta_flux, h_beta_flux_err,
                     nii_6584_flux, nii_6584_flux_err, h_alpha_flux, h_alpha_flux_err, class
                   from Photoz as p
                     full outer join SpecObj as s on p.objID = s.bestObjID
                     full outer join GalSpecLine as g on s.specobjid = g.specobjid
                   where p.objID = %i''' % sdss_match['objID']
        df = SDSS_select(query)
        if len(df):
            more_data = {}
            if not np.isnan(df['spec_redshift'][0]):
                more_data['spec_redshift'] = df['spec_redshift'][0]
                more_data['spec_redshift_err'] = df['spec_redshift_err'][0]
            if df['photo_redshift'][0] != -9999:
                more_data['photo_redshift'] = df['photo_redshift'][0]
                more_data['photo_redshift_err'] = df['photo_redshift_err'][0]
            if type(df['class'][0]) is not np.float64 or not np.isnan(df['class'][0]):
                more_data['spectral_class'] = df['class'][0][0]
            for key in ['oiii_5007_flux', 'oiii_5007_flux_err', 'h_beta_flux', 'h_beta_flux_err', \
                        'nii_6584_flux', 'nii_6584_flux_err', 'h_alpha_flux', 'h_alpha_flux_err']:
                if not np.isnan(df[key][0]):
                    more_data[key] = df[key][0]
            sdss_match.update(more_data)

    elif sdss_match and sdss_match['morphological_class'] == 'S': #query the star tables

        query = '''select so.z as spec_redshift, so.zErr as spec_redshift_err, class
                     from Star as s
                       full outer join SpecObj as so on s.objID=so.bestObjID
                   where s.objID = %i''' % sdss_match['objID']
        df = SDSS_select(query)
        if len(df):
            more_data = {}
            if not np.isnan(df['spec_redshift'][0]):
                more_data['spec_redshift'] = df['spec_redshift'][0]
                more_data['spec_redshift_err'] = df['spec_redshift_err'][0]
            if type(df['class'][0]) is not np.float64 or not np.isnan(df['class'][0]):
                more_data['spectral_class'] = df['class'][0][0]
            sdss_match.update(more_data)

    if sdss_match:
        logging.info('SDSS match found')
        for key in sdss_match.keys():
            if sdss_match[key] is None:
                sdss_match.pop(key)
            elif sdss_match[key] and type(sdss_match[key]) is not str:
                sdss_match[key] = sdss_match[key].item()
            elif sdss_match[key] == 0:
                sdss_match[key] = 0
    else:
        logging.info('No SDSS match found')
    
    return sdss_match

def getRadio(data, fits_loc, consensusObject):
    '''
    calculates all of the radio parameters from the fits file
    data is a JSON object downloaded from the online RGZ interface
    fits_loc is the fits file on the physical drive
    '''
    
    #create list of trees, each containing a contour and its contents
    contourTrees = []
    consensusBboxes = consensusObject['bbox']
    for contour in data['contours']:
        for bbox in consensusBboxes:
            if fn.approx(contour[0]['bbox'][0], bbox[0]) and fn.approx(contour[0]['bbox'][1], bbox[1]) and \
               fn.approx(contour[0]['bbox'][2], bbox[2]) and fn.approx(contour[0]['bbox'][3], bbox[3]):
                tree = c.Node(contour=contour, fits_loc=fits_loc)
                contourTrees.append(tree)
    
    #get component fluxes and sizes
    components = []
    for tree in contourTrees:
        bboxP = fn.bboxToDS9(fn.findBox(tree.value['arr']), tree.imgSize)[0] #bbox in DS9 coordinate pixels
        bboxCornersRD = tree.w.wcs_pix2world( np.array( [[bboxP[0],bboxP[1]], [bboxP[2],bboxP[3]] ]), 1) #two opposite corners of bbox in ra and dec
        raRange = [ min(bboxCornersRD[0][0], bboxCornersRD[1][0]), max(bboxCornersRD[0][0], bboxCornersRD[1][0]) ]
        decRange = [ min(bboxCornersRD[0][1], bboxCornersRD[1][1]), max(bboxCornersRD[0][1], bboxCornersRD[1][1]) ]
        pos1 = coord.SkyCoord(raRange[0], decRange[0], unit=(u.deg, u.deg))
        pos2 = coord.SkyCoord(raRange[1], decRange[1], unit=(u.deg, u.deg))
        extentArcsec = pos1.separation(pos2).arcsecond
        solidAngleArcsec2 = tree.areaArcsec2
        components.append({'flux':tree.fluxmJy, 'flux_err':tree.fluxErrmJy, 'angular_extent':extentArcsec, 'solid_angle':solidAngleArcsec2, \
                           'ra_range':raRange, 'dec_range':decRange})
    
    #adds up total flux of all components
    totalFluxmJy = 0
    totalFluxErrmJy2 = 0
    for component in components:
        totalFluxmJy += component['flux']
        totalFluxErrmJy2 += np.square(component['flux_err'])
    totalFluxErrmJy = np.sqrt(totalFluxErrmJy2)
    
    #finds total area enclosed by contours in arcminutes
    totalSolidAngleArcsec2 = 0
    for component in components:
        totalSolidAngleArcsec2 += component['solid_angle']
    
    #find maximum extent of component bboxes in arcseconds
    maxAngularExtentArcsec = 0
    if len(components)==1:
        maxAngularExtentArcsec = components[0]['angular_extent']
    else:
        for i in range(len(components)-1):
            for j in range(1,len(components)-i):
                corners1 = np.array([ [components[i]['ra_range'][0], components[i]['dec_range'][0]], \
                                      [components[i]['ra_range'][0], components[i]['dec_range'][1]], \
                                      [components[i]['ra_range'][1], components[i]['dec_range'][0]], \
                                      [components[i]['ra_range'][1], components[i]['dec_range'][1]] ])
                corners2 = np.array([ [components[i+j]['ra_range'][0], components[i+j]['dec_range'][0]], \
                                      [components[i+j]['ra_range'][0], components[i+j]['dec_range'][1]], \
                                      [components[i+j]['ra_range'][1], components[i+j]['dec_range'][0]], \
                                      [components[i+j]['ra_range'][1], components[i+j]['dec_range'][1]] ])
                pos1 = coord.SkyCoord(corners1.T[0], corners1.T[1], unit=(u.deg, u.deg))
                pos2 = coord.SkyCoord(corners2.T[0], corners2.T[1], unit=(u.deg, u.deg))
                angularExtentArcsec = pos1.separation(pos2).arcsecond
                maxAngularExtentArcsec = max(np.append(angularExtentArcsec, maxAngularExtentArcsec))
    
    #add all peaks up into single list
    peakList = []
    for tree in contourTrees:
        for peak in tree.peaks:
            peakList.append(peak)
    peakFluxErrmJy = contourTrees[0].sigmamJy

    #find center of radio source
    raMin, raMax, decMin, decMax = np.inf, 0, np.inf, 0
    for comp in components:
        if comp['ra_range'][0] < raMin:
            raMin = comp['ra_range'][0]
        if comp['ra_range'][1] > raMax:
            raMax = comp['ra_range'][1]
        if comp['dec_range'][0] < decMin:
            decMin = comp['dec_range'][0]
        if comp['dec_range'][1] > decMax:
            decMax = comp['dec_range'][1]
    meanRa = (raMax+raMin)/2.
    meanDec = (decMax+decMin)/2.
    
    radio_data = {'radio':{'total_flux':totalFluxmJy, 'total_flux_err':totalFluxErrmJy, 'outermost_level':data['contours'][0][0]['level']*1000, \
                           'number_components':len(contourTrees), 'number_peaks':len(peakList), 'max_angular_extent':maxAngularExtentArcsec, \
                           'total_solid_angle':totalSolidAngleArcsec2, 'peak_flux_err':peakFluxErrmJy, 'peaks':peakList, 'components':components, \
                           'ra':meanRa, 'dec':meanDec}}
    
    return radio_data
