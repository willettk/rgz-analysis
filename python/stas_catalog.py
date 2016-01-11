from pymongo import MongoClient
import pandas as pd

path = '/Users/willettk/Astronomy/Research/GalaxyZoo'

consensus_level = 0.50

if consensus_level == 0.:
    suffix = '_full'
else:
    suffix = ''

def load_data():

    client = MongoClient('localhost', 27017)
    db = client['radio'] 
    
    catalog = db['catalog']

    return catalog

def run_static(catalog,full=False):
    filename = '%s/rgz-analysis/csv/stas_rgz%s.csv' % (path,suffix)

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

        # Header
        print >> f,'source_id,zooniverse_id,n_radio,ra_min,ra_max,dec_min,dec_max,max_angular_extent,total_solid_angle,wise_id,wise_ra,wise_dec,wise_w1mag,redshift,redshift_err,redshift_type,sdss_id,sdss_ra,sdss_dec,consensus_level'
        good_entry,bad_entry = 0,0
        args = {'consensus.level':{"$gte":consensus_level},'SDSS':{'$exists':True},'radio.numberComponents':{'$gte':2}}
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
                print >> f,'RGZ_{0:},{14:},{21:d},{22:.5f},{23:.5f},{24:.5f},{25:.5f},{17:.3f},{18:.3f},{10:},{7:.5f},{8:.5f},{9:.2f},{11:.4f},{12:.4f},{13:d},{16:},{19:.5f},{20:.5f},{15:.2f}'.format(\
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
                bad_entry += 1

        print "{0:d} entries written to CSV file".format(good_entry)
        print "{0:d}/{1:d} had errors writing data to file".format(bad_entry,catalog.find(args).count())

    return None

if __name__ == "__main__":
    catalog = load_data()
    run_static(catalog)
