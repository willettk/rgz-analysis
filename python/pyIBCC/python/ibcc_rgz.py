# Produce test data sets for use of RGZ data with Edwin Simpson's pyIBCC routine

from pymongo import MongoClient

def load_rgz_data():

    # Connect to Mongo database
    # Make sure to run mongorestore /path/to/database to restore the updated files
    # mongod client must be running locally
    
    client = MongoClient('localhost', 27017)
    db = client['ouroboros'] 
    
    subjects = db['radio_subjects'] 		# subjects = images
    classifications = db['radio_classifications']	# classifications = classifications of each subject per user
    users = db['radio_users']	# volunteers doing each classification (can be anonymous)

    return subjects,classifications,users

def make_rgz_sample(subjects):

    # Produce a sample dataset of 100 galaxies
    
    '''
    # Idea: make sure each classified at least once by our four top users?
    top_users = ('antikodon','planetari7','pamelaann','WizardHowl')
    '''

    n_subjects = 100

    zooniverse_ids = []
    for sub in list(subjects.find({'classification_count':{'$eq':20}}))[:n_subjects]:
        zooniverse_ids.append(sub['zooniverse_id'])

    return zooniverse_ids

# Call program from command line

if __name__ == '__main__':

    subjects,classifications,users = load_rgz_data()
    zooniverse_ids = make_rgz_sample(subjects)

