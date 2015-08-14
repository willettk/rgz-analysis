from pymongo import MongoClient

path = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

def load_data():

    client = MongoClient('localhost', 27017)
    db = client['rgz'] 
    
    catalog = db['catalog']

    return catalog

def run_static(catalog):
    filename = '%s/csv/static_catalog.csv'
    with open(filename,'w') as f:

        # Header
        print >> f,'source_id peak1_ra peak1_dec peak1_flux peak2_ra peak2_dec peak2_flux wise_ra wise_dec wise_w1mag'
        # Data requested by Larry for double-peaked sources
        for c in catalog.find({'radio.numberComponents':2,'AllWISE':{'$exists':True}}):
            print >> f,'RGZ_%i' % c['catalog_id'],c['radio']['peaks'][0]['ra'],c['radio']['peaks'][0]['dec'],c['radio']['peaks'][0]['flux'],c['radio']['peaks'][1]['ra'],c['radio']['peaks'][1]['dec'],c['radio']['peaks'][1]['flux'],c['AllWISE']['ra'],c['AllWISE']['dec'],c['AllWISE']['w1mpro']

if __name__ == "__main__":
    catalog = load_data()
    run_static(catalog)
