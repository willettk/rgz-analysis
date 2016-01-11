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
    filename = '%s/rgz-analysis/csv/static_catalog%s.csv' % (path,suffix)

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
        print >> f,'source_id zooniverse_id peak1_ra peak1_dec peak1_flux peak2_ra peak2_dec peak2_flux max_angular_extent total_solid_angle wise_designation wise_ra wise_dec wise_w1mag redshift redshift_err redshift_type sdss_id sdss_ra sdss_dec consensus_level'
        # Data requested by Larry for double-peaked sources
        #for c in catalog.find({'radio.numberComponents':2}):
        bad_entry = 0
        for c in catalog.find({'consensus.level':{"$gte":consensus_level}}):
            wiseval = c.setdefault('AllWISE',wise_default_dict)
            sdssval = c.setdefault('SDSS',sdss_default_dict)
            sdssredshift = c['SDSS'].setdefault('redshift',-99.)
            sdssredshifterr = c['SDSS'].setdefault('redshift_err',-99.)
            sdssredshifttype = c['SDSS'].setdefault('redshift_type',-99)

            try:
                print >> f,'RGZ_{0:} {14:} {1:.5f} {2:.5f} {3:.2f} {4:.5f} {5:.5f} {6:.2f} {17:.3f} {18:.3f} {10:} {7:.5f} {8:.5f} {9:.2f} {11:.4f} {12:.4f} {13:d} {16:} {19:.5f} {20:.5f} {15:.2f}'.format(\
                c['catalog_id'], 
                c['radio']['peaks'][0]['ra'],
                c['radio']['peaks'][0]['dec'],
                c['radio']['peaks'][0]['flux'],
                c['radio']['peaks'][1]['ra'],
                c['radio']['peaks'][1]['dec'],
                c['radio']['peaks'][1]['flux'],
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
                c['SDSS']['dec'])
            except IndexError:
                bad_entry += 1

        print "{0:d}/{1:d} had no results for radio, SDSS, or WISE".format(bad_entry,catalog.find().count())

    return None

def match_clusters():

    df1 = pd.read_csv('%s/radiogalaxyzoo/cluster_matching/MATCHED_PAIRS.tsv' % path,delim_whitespace=True)
    df2 = pd.read_csv('%s/rgz-analysis/csv/static_catalog%s.csv' % (path,suffix),delim_whitespace=True)

    # Keep only columns that Larry is interested in

    allcols = set(df1.columns)
    keep = set((u'cluster_id', u'cluster_ra', u'cluster_dec', u'rgz_id',u'cluster_best_z', u'cluster_rl', u'Separation'))
    df1.drop(list(allcols.difference(keep)),axis=1,inplace=True)
    df1.rename(columns={'Separation':'projected_sep'},inplace=True)

    # Rename the columns so I can match them against each other
    df2.rename(columns={'source_id':'rgz_id'},inplace=True)
    df2.rgz_id.replace("_","",regex=True,inplace=True)

    dfm = df2.merge(df1,on='rgz_id')

    # Saves file as hard copy, but not directly used
    dfm.to_csv('%s/rgz-analysis/csv/static_catalog2%s.csv' % (path,suffix),sep = ' ', index=False)

    return dfm

def load_angles(filename):

    df = pd.read_csv('%s/rgz-analysis/bending_angles/%s.csv' % (path,filename))
    df['angle_type'] = pd.Series((filename[7:],)*len(df), index=df.index)
    
    return df

def match_bending_angle(dfm):

    df1 = load_angles('angles_double_pixradio')
    #df2 = load_angles('angles_triple_pixradio')
    df3 = load_angles('angles_multipeaked_singles')
    #df4 = load_angles('angles_multipeaked_singles_no_optical')

    for x in (df1,df3):
        print len(x),x['angle_type'][0]

    dfall = pd.concat([df1,df3],ignore_index=True)

    dfba = dfm.merge(dfall,on='zooniverse_id')
    dfba.to_csv('%s/rgz-analysis/csv/static_catalog3%s.csv' % (path,suffix),sep = ' ', index=False)

    print "\n%i sources in static catalog" % len(dfba)

    return None
    
if __name__ == "__main__":
    catalog = load_data()
    run_static(catalog)
    dfm = match_clusters()
    match_bending_angle(dfm)
