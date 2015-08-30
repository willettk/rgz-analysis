import logging
from pymongo import MongoClient
import csv

def updateConsensus(csvPath):

    logging.basicConfig(filename='RGZcatalog.log', level=logging.DEBUG, format='%(asctime)s %(message)s')
    logging.captureWarnings(True)
    logging.info('New consensus collection added from %s', csvPath)

    with open(csvPath, 'r') as csvFile:
        
        db = MongoClient()['radio']
        db.drop_collection('consensus')
        consensus = db['consensus']

        consensusDict = csv.DictReader(csvFile)
        header = ['zooniverse_id', 'FIRST_id', 'n_users', 'n_total', 'consensus_level', 'n_radio', 'label', 'bbox', 'ir_peak']
        for entry in consensusDict:
            row = {}
            for field in header:
                row[field] = entry[field]
            consensus.insert(row)

        logging.info('%i entries added to consensus collection', consensus.count())
