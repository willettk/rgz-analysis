# Reduce RGZ data with Edwin Simpson's pyIBCC routine

import rgz
import consensus
from astropy.io import ascii

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

def rgz_sample_input(subjects,classifications):

    # Try 10 galaxies from our top 10 users?
    
    n_subjects = 5
    n_users = 5

    topusers = ('antikodon', 'WizardHowl', 'kijkuit', 'planetari7', 'pamelaann', 'arianna99', 'orionsam1', 'mattson', 'bjornstromjennifet', 'snowysky')[:n_users]
    
    inputfile = open('%s/python/pyIBCC/data/rgz/input.csv' % rgz_dir,'wb')

    subject_dict = {}
    master_subcount = 0

    for userID,tu in enumerate(topusers):
        # Need userID, subjectID, score.
        
        # List of their classifications
        cdone = classifications.find({'user_name':tu})
        scount = 0

        # Find classifications they've done on completed subjects
        for c in list(cdone):
            s = subjects.find_one({'_id':c['subject_ids'][0]})
            if s['state'] == 'complete':
                scount += 1

                sid = s['_id']
                if sid not in subject_dict.keys():
                    subject_dict[sid] = master_subcount
                    master_subcount += 1

                subjectID = subject_dict[sid]
            
                oa = consensus.one_answer(s['zooniverse_id'],tu)
                score = len(oa['answer'])
 
                print '%i    -    http://radiotalk.galaxyzoo.org/#/subjects/%s' % (subjectID,s['zooniverse_id'])

                # For the gold standard data, see if any of these were in the expert 100 sample?

                inputfile.write('%i,%i,%i\n' % (userID,subjectID,score))

            if scount >= n_subjects:
                break

        print 'Wrote %i galaxies to input.csv for %s' % (n_subjects,tu) 

    inputfile.close()

    return None

def rgz_goldstandard_input(subjects,classifications):

    # Get data for all users who did galaxies in the gold standard sample. 
    # Compare this to the raw votes from consensus.py !!

    # Here, score is number of radio components in image, not number of galaxies

    gs = ascii.read('%s/goldstandard/gs_zids.txt' % rgz_dir,names=['name'],data_start=0)
    #gs = ascii.read('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir,names=['name'],data_start=0)
    
    inputfile = open('%s/python/pyIBCC/data/rgz/gs_input.csv' % rgz_dir,'wb')

    user_dict = {}
    master_usercount = 0

    experts=("42jkb", "ivywong", "stasmanian", "klmasters", "Kevin", "akapinska", "enno.middelberg", "xDocR", "DocR", "vrooje", "KWillett")

    for subjectID,gal in enumerate(gs['name']):

        # Find all non-expert users who classified the image

        sid = subjects.find_one({'zooniverse_id':gal})['_id']
        clist = classifications.find({'subject_ids':sid,'user_name':{'$nin':experts,'$exists':True}})

        for c in list(clist):

            try:
                uname = c['user_name']
            except KeyError:
                uname = 'Anonymous'

            if uname not in user_dict.keys():
                user_dict[uname] = master_usercount
                master_usercount += 1

            userID = user_dict[uname]

            oa = consensus.one_answer(gal,uname)
            score = 0
            for ans in oa['answer'].itervalues():
                score = max((score,len(ans['xmax'])))

            inputfile.write('%i,%i,%i\n' % (userID,subjectID,score))
        

        print 'Wrote %i classifications to input file for %s' % (clist.count(),gal) 

    inputfile.close()

    return None

# Call program from command line

if __name__ == '__main__':

    subjects,classifications = rgz.load_rgz_data()
    #rgz_sample_input(subjects,classifications)
    rgz_goldstandard_input(subjects,classifications)

