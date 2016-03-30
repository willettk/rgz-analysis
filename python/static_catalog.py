from pymongo import MongoClient

# Make a version of the Radio Galaxy Zoo catalog that's perusable as a flat FITS or CSV table

# Set the desired consensus level of the RGZ classifications (between 0 and 1)

path = '.'
consensus_level = 0.50

# Define a filestem that will be appended to the new catalog output

if consensus_level == 0.:
    suffix = '_full'
else:
    suffix = ''

def load_data():

    # Loads in the matched catalog from MongoDB

    client = MongoClient('localhost', 27017)
    db = client['radio'] 
    
    catalog = db['catalog']

    return catalog

def flat_version(catalog,full=False):

    # Write the MongoDB catalog to a CSV file including all fields with a fixed number of elements

    filename = '%s/csv/static_rgz_flat%s.csv' % (path,suffix)

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

        header = ''
        header += 'catalog_id,zooniverse_id,first_id'
        ex = catalog.find_one({'AllWISE':{'$exists':True},'SDSS':{'$exists':True}})
        consensus_keys = [str(x) for x in ex['consensus'].keys() if x not in ('label')]
        sdss_keys  = [str(x) for x in ex['SDSS'].keys()]
        wise_keys  = [str(x) for x in ex['AllWISE'].keys()]
        radio_keys = [str(x) for x in ex['radio'].keys() if x not in ('components','peaks')]

        arrs = (consensus_keys,radio_keys,wise_keys,sdss_keys)
        labels = ('consensus','radio','wise','sdss')

        for arr,label in zip(arrs,labels):
            arr.sort()
            extended_arr = ['']
            extended_arr.extend(arr)
            header += (',%s.'%label).join(extended_arr)

        print >> f,header
        good_entry,bad_entry = 0,0

        # Select all matching galaxies (in this case, sources with optical and IR counterparts)
        args = {'consensus.level':{"$gte":consensus_level},'SDSS':{'$exists':True},'AllWISE':{'$exists':True}}

        # Loop over number of RGZ catalog entries that match the consensus requirements
        for c in catalog.find(args):
            wiseval = c.setdefault('AllWISE',wise_default_dict)
            sdssval = c.setdefault('SDSS',sdss_default_dict)
            sdssredshift = c['SDSS'].setdefault('redshift',-99.)
            sdssredshifterr = c['SDSS'].setdefault('redshift_err',-99.)
            sdssredshifttype = c['SDSS'].setdefault('redshift_type',-99)

            try:
                # Print all values to new row in file. 

                row = ['RGZ_'+str(c['catalog_id']),c['zooniverse_id'],c['first_id']]
                for entry in header.split(',')[3:]:
                    bothvar = entry.split('.')
                    if bothvar[0] == 'wise':
                        bothvar[0] = 'AllWISE'
                    if bothvar[0] == 'sdss':
                        bothvar[0] = 'SDSS'
                    try:
                        row.append(c[bothvar[0]][bothvar[1]])
                    except KeyError:
                        row.append(-99)

                prow = [str(x) for x in row]
                print >> f,','.join(prow)
                good_entry += 1
            except IndexError:
                # If couldn't find one or more of the fields selected
                bad_entry += 1

        # Print summary to screen

        print "{0:d} entries written to CSV file {1}".format(good_entry,filename)
        print "{0:d}/{1:d} had errors writing data to file".format(bad_entry,catalog.find(args).count())

    return None

def selected_fields(catalog,full=False):

    # Write a version of the MongoDB catalog to a CSV file. Includes only fields that are manually selected here. 

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

    # Run catalog from command line
    catalog = load_data()
    flat_version(catalog)
    #selected_fields(catalog)
