from pymongo import MongoClient

# Make a version of the Radio Galaxy Zoo catalog that's perusable as a flat FITS or CSV table. 
# This is based on the output of:
#   consensus.py
#   RGZcatalog.py

'''
Set the desired consensus level of the RGZ classifications for this output (between 0 and 1).
0 will include every source in the catalog; 1 would include only those with 100% consensus. 
Default right now is (arbitrarily) set at 50%. One could also output the entire catalog
with consensus_level = 0 and then perform cuts later in their analysis, since consensus_level
is included as an output parameter.
'''

path = '.'
consensus_level = 0.0

# Define a suffix that will be appended to the filename of the new catalog output.

if consensus_level == 0.:
    suffix = '_full'
else:
    suffix = ''

def load_data():

    # Load the matched catalog from MongoDB

    client = MongoClient('localhost', 27017)
    db = client['radio'] 
    
    catalog = db['catalog']

    return catalog

def flat_version(catalog,full=False):

    # Write the MongoDB catalog to a CSV file, including only sources with a WISE and SDSS match.
    # Include all data fields.

    filename = '%s/csv/static_rgz_flat%s.csv' % (path,suffix)

    # Check the WISE and SDSS fields to see if they had a match; if not, return null values
    
    wise_default_dict = {'designation':-99, 'ra':-99, 'dec':-99, 'numberMatches':0, \
                          'w1mpro':-99, 'w1sigmpro':-99, 'w1snr':-99, 'w2mpro':-99, 'w2sigmpro':-99, 'w2snr':-99, \
                          'w3mpro':-99, 'w3sigmpro':-99, 'w3snr':-99, 'w4mpro':-99, 'w4sigmpro':-99, 'w4snr':-99}
    
    sdss_default_dict = {'objID':-99, 'ra':-99, 'dec':-99, 'numberMatches':0, 'redshift':-99, 'redshift_err':-99, 'redshift_type':-99, \
                         'u':-99, 'g':-99, 'r':-99, 'i':-99, 'z':-99, 'u_err':-99, 'g_err':-99, 'r_err':-99, 'i_err':-99, 'z_err':-99}

    with open(filename,'w') as f:

        # CSV file header

        header = 'catalog_id,rgz_name,zooniverse_id,first_id'
        consensus_keys = ['IR_ra', 'IR_dec', 'level', 'n_votes', 'n_total']
        sdss_keys  = [str(x) for x in sdss_default_dict.keys()]
        wise_keys  = [str(x) for x in wise_default_dict.keys()]
        radio_keys = ['ra', 'dec', 'totalFlux', 'totalFluxErr', 'outermostLevel', 'numberComponents', 'numberPeaks', 'maxAngularExtent', 'totalSolidAngle', \
                      'peakFluxErr', 'maxPhysicalExtent', 'totalCrossSection', 'totalLuminosity', 'totalLuminosityErr', 'peakLuminosityErr']
        component_keys = ['fluxes', 'fluxErrs', 'peakRas', 'peakDecs']
        peak_keys = ['fluxes', 'ras', 'decs']
        duplicate_keys = ['shareComponents', 'matchComponents', 'WISECATmismatch']

        arrs = (consensus_keys,radio_keys,component_keys,peak_keys,wise_keys,sdss_keys,duplicate_keys)
        labels = ('consensus','radio','components','peaks','wise','sdss','duplicateSources')

        for arr,label in zip(arrs,labels):
            arr.sort()
            extended_arr = ['']
            extended_arr.extend(arr)
            header += (',%s.'%label).join(extended_arr)

        print >> f,header
        good_entry,bad_entry = 0,0

        # Find all duplicate sources for deletion
        cids_for_removal = []
        for c in catalog.find({'duplicateSources.exactDuplicate':{'$exists':True}}):
            if c['catalog_id'] != min(c['duplicateSources']['exactDuplicate']):
                cids_for_removal.append(c['catalog_id'])
        
        # Select all matching galaxies (in this case, sources with optical and IR counterparts)
        args = {'catalog_id':{'$nin':cids_for_removal},'consensus.level':{'$gte':consensus_level}}#,'SDSS':{'$exists':True},'AllWISE':{'$exists':True}}
        
        # Loop over number of RGZ catalog entries that match the consensus requirements
        for c in catalog.find(args).sort([('catalog_id', 1)]):
            
            wiseval = c.setdefault('AllWISE',wise_default_dict)
            sdssval = c.setdefault('SDSS',sdss_default_dict)
            #sdssredshift = c['SDSS'].setdefault('redshift',-99.)
            #sdssredshifterr = c['SDSS'].setdefault('redshift_err',-99.)
            #sdssredshifttype = c['SDSS'].setdefault('redshift_type',-99)

            # Determine component strings
            component_strings = {'fluxes':'', 'fluxErrs':'', 'peakRas':'', 'peakDecs':''}
            for component in c['radio']['components']:
                component_strings['fluxes'] += '{0};'.format(component['flux'])
                component_strings['fluxErrs'] += '{0};'.format(component['fluxErr'])
                maxPeak = {'flux':-99, 'ra':-99, 'dec':-99}
                for peak in c['radio']['peaks']:
                    if component['raRange'][0] <= peak['ra'] <= component['raRange'][1] and \
                       component['decRange'][0] <= peak['dec'] <= component['decRange'][1] and \
                       peak['flux'] > maxPeak['flux']:
                        maxPeak = peak.copy()
                component_strings['peakRas'] += '{0};'.format(maxPeak['ra'])
                component_strings['peakDecs'] += '{0};'.format(maxPeak['dec'])
            for key in component_strings:
                component_strings[key] = component_strings[key][:-1]

            # Determine peak strings (only for single component sources)
            peak_strings = {'fluxes':'', 'ras':'', 'decs':''}
            if c['radio']['numberComponents'] == 1:
                for peak in c['radio']['peaks']:
                    peak_strings['fluxes'] += '{0};'.format(peak['flux'])
                    peak_strings['ras'] += '{0};'.format(peak['ra'])
                    peak_strings['decs'] += '{0};'.format(peak['dec'])
                for key in peak_strings:
                    peak_strings[key] = peak_strings[key][:-1]
            else:
                for key in peak_strings:
                    peak_strings[key] = '-99'

            # Determine overlap strings (when applicable)
            duplicate_strings = {'shareComponents':'', 'matchComponents':'', 'WISECATmismatch':''}
            if 'duplicateSources' in c:
                for key in duplicate_strings:
                    if key in c['duplicateSources']:
                        for cid in c['duplicateSources'][key]:
                            if cid not in cids_for_removal:
                                duplicate_strings[key] += '{0};'.format(cid)
                        duplicate_strings[key] = duplicate_strings[key][:-1]
            for key in duplicate_strings:
                if duplicate_strings[key] == '':
                    duplicate_strings[key] = -99

            # Combine sources (if duplicates exist)
            if 'duplicateSources' in c and 'exactDuplicate' in c['duplicateSources']:
                votes, total = 0, 0
                IR_ra, IR_dec = 0., 0.
                for d in catalog.find({'catalog_id': {'$in': c['duplicateSources']['exactDuplicate']}}):
                    votes += d['consensus']['n_votes']
                    total += d['consensus']['n_total']
                    if 'IR_ra' in d['consensus']:
                        IR_ra += d['consensus']['n_votes'] * d['consensus']['IR_ra']
                        IR_dec += d['consensus']['n_votes'] * d['consensus']['IR_dec']
                IR_ra /= votes
                IR_dec /= votes
                c['consensus']['n_votes'] = votes
                c['consensus']['n_total'] = total
                if IR_ra:
                    c['consensus']['IR_ra'] = IR_ra
                    c['consensus']['IR_dec'] = IR_dec
            
            try:
                # Print all values to new row in file. 
                row = [c['catalog_id'],c['rgz_name'],c['zooniverse_id'],c['first_id']]
                for entry in header.split(',')[4:]:
                    bothvar = entry.split('.')
                    if bothvar[0] == 'wise':
                        bothvar[0] = 'AllWISE'
                    if bothvar[0] == 'sdss':
                        bothvar[0] = 'SDSS'
                    try:
                        if bothvar[0] == 'components':
                            row.append(component_strings[bothvar[1]])
                        elif bothvar[0] == 'peaks':
                            row.append(peak_strings[bothvar[1]])
                        elif bothvar[0] == 'duplicateSources':
                            row.append(duplicate_strings[bothvar[1]])
                        else:
                            row.append(c[bothvar[0]][bothvar[1]])
                    except KeyError:
                        row.append(-99)

                prow = [str(x) for x in row]
                print >> f,','.join(prow)
                good_entry += 1
                
            except IndexError:
                # If couldn't find one or more of the fields selected
                bad_entry += 1
                print "Unable to print {0}".format(c['catalog_id'])

        # Print summary to screen

        print "{0:d} entries written to CSV file {1}".format(good_entry,filename)
        print "{0:d}/{1:d} had errors writing data to file".format(bad_entry,catalog.find(args).count())

    return None

#Hasn't been updated with new fields
def selected_fields(catalog,full=False):

    # Write the MongoDB catalog to a CSV file, but only include specific fields.

    filename = '%s/csv/static_rgz_selected%s.csv' % (path,suffix)

    # Check the WISE and SDSS fields to see if they had a match; if not, return null values

    wise_default_dict = catalog.find_one({'AllWISE':{'$exists':True}})['AllWISE']
    for k in wise_default_dict:
        wise_default_dict[k] = -99.
    wise_default_dict['designation'] = 'no_wise_match'
    wise_default_dict['numberMatches'] = 0

    sdss_default_dict = catalog.find_one({'SDSS':{'$exists':True},'SDSS.redshift':{'$exists':False}})['SDSS']
    for k in sdss_default_dict:
        sdss_default_dict[k] = -99.
    sdss_default_dict['objID'] = 'no_sdss_match'
    sdss_default_dict['numberMatches'] = 0

    with open(filename,'w') as f:

        # CSV file header
        print >> f,'source_id,zooniverse_id,n_radio,ra_min,ra_max,dec_min,dec_max,max_angular_extent,total_solid_angle,wise_id,wise_ra,wise_dec,wise_w1mag,redshift,redshift_err,redshift_type,sdss_id,sdss_ra,sdss_dec,consensus_level'
        good_entry,bad_entry = 0,0

        # Select all matching galaxies (in this case, double-lobed sources with optical and IR counterparts)
        args = {'consensus.level':{"$gte":consensus_level},'SDSS':{'$exists':True},'AllWISE':{'$exists':True}}

        # Loop over number of RGZ catalog entries that match the consensus requirements
        for c in catalog.find(args):
            wiseval = c.setdefault('AllWISE',wise_default_dict)
            sdssval = c.setdefault('SDSS',sdss_default_dict)
            sdssredshift = c['SDSS'].setdefault('redshift',-99.)
            sdssredshifterr = c['SDSS'].setdefault('redshift_err',-99.)
            sdssredshifttype = c['SDSS'].setdefault('redshift_type',-99)

            # Find the maximum RA, dec extent on the radio components

            components = c['radio']['components']
            ra_min,ra_max = components[0]['raRange']
            dec_min,dec_max = components[0]['decRange']
            for comp in components:
                ra_min = ra_min if comp['raRange'][0] > ra_min else comp['raRange'][0]
                ra_max = ra_max if comp['raRange'][1] < ra_max else comp['raRange'][1]
                dec_min = dec_min if comp['decRange'][0] > dec_min else comp['decRange'][0]
                dec_max = dec_max if comp['decRange'][1] < dec_max else comp['decRange'][1]

            try:
                # Print values to new row in file. -99 values for the first few roles can be removed if the fields are renumbered in this line.
                print >> f,'RGZ_{0},{14},{21:d},{22:.5f},{23:.5f},{24:.5f},{25:.5f},{17:.3f},{18:.3f},{10},{7:.5f},{8:.5f},{9:.2f},{11:.4f},{12:.4f},{13:d},{16},{19:.5f},{20:.5f},{15:.2f}'.format(\
                c['catalog_id'], 
                -99,
                -99,
                -99,
                -99,
                -99,
                -99,
                c['AllWISE']['ra'],
                c['AllWISE']['dec'],
                c['AllWISE']['w1mpro'],
                c['AllWISE']['designation'],
                c['SDSS']['redshift'],
                c['SDSS']['redshift_err'],
                c['SDSS']['redshift_type'],
                c['Zooniverse_id'],
                c['consensus']['level'],
                c['SDSS']['objID'],
                c['radio']['maxAngularExtent'],
                c['radio']['totalSolidAngle'],
                c['SDSS']['ra'],
                c['SDSS']['dec'],
                c['radio']['numberComponents'],
                ra_min,
                ra_max,
                dec_min,
                dec_max)
                good_entry += 1
            except IndexError:
                # If couldn't find one or more of the fields selected
                bad_entry += 1

        # Print summary to screen

        print "{0:d} entries written to CSV file {1}".format(good_entry,filename)
        print "{0:d}/{1:d} had errors writing data to file".format(bad_entry,catalog.find(args).count())

    return None

if __name__ == "__main__":

    # Make the static catalog from the command line

    catalog = load_data()
    flat_version(catalog)
