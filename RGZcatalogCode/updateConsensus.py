'''
Updates the consensus database in Mongo so that new entries may be added to the catalog.
'''

import logging
from pymongo import MongoClient
import csv
from ast import literal_eval

csvPath = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/csv/consensus_rgz_first_75.csv'

def updateConsensus(csvPath):

    logging.basicConfig(filename='RGZcatalog.log', level=logging.DEBUG, format='%(asctime)s %(message)s')
    logging.captureWarnings(True)
    logging.info('New consensus collection added from %s', csvPath)

    with open(csvPath, 'r') as csvFile:
        
        #db = MongoClient()['radio']
        #db.drop_collection('consensus')
        #consensus = db['consensus']
        db = MongoClient()['platypus']
        db.drop_collection('test')
        test = db['test']

        consensusDict = csv.DictReader(csvFile)
        header = ['zooniverse_id', 'first_id', 'n_users', 'n_total', 'consensus_level', 'n_radio', 'label', 'bbox', 'ir_peak']
        for entry in consensusDict:
            row = {}
            for field in header:
                try:
                    entry_typed = literal_eval(entry[field])
                except (ValueError,SyntaxError) as e:
                    entry_typed = str(entry[field])
                row[field] = entry_typed
            test.insert(row)

        logging.info('%i entries added to test collection', test.count())

if __name__ == "__main__":
    updateConsensus(csvPath)
