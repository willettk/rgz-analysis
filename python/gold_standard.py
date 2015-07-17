# How many galaxies from the gold standard sample are necessary to establish whether a user agrees with the experts?

# Hm - not many have done more than 1. Might have to expand the experiment to the full expert data set of 100 galaxies ...

# Plot (N/2 - M/2) vs. (M+N) and see where agreement lies!!

# Read in the list of gold standard images

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
import consensus
import rgz
import operator
import numpy as np
import cPickle as pickle
from matplotlib import pyplot as plt
import scipy.stats.distributions as dist
from collections import Counter

subjects,classifications,users = rgz.load_rgz_data()

experts=("42jkb", "ivywong", "stasmanian", "klmasters", "Kevin", "akapinska", "enno.middelberg", "xDocR", "DocR", "vrooje", "KWillett")
    
def get_galaxies():
    
    #with open('%s/goldstandard/gs_zids.txt' % rgz_dir) as f:
    with open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir) as f:
        gs = f.readlines()
    
    gals = [g.strip() for g in gs]

    return gals
    
def top_volunteers(gs):

    # How many unique users have done each galaxy in the expert sample?
    
    allusers = []
    
    for g in gs:
    
        s = subjects.find_one({'zooniverse_id':g})
        c = classifications.find({'subject_ids.0':s['_id']})
        
        ulist = []
        for cc in c:
            if cc.has_key('user_name'):
                ulist.append(cc['user_name'])
            else:
                ulist.append('Anonymous')
    
        setu = set(ulist)
    
        for u in setu:
            if u not in experts:
                allusers.append(u)
    
    cntr_users = Counter(allusers)
    
    # Limit to users who have done 10 or more of the expert sample
    
    ud = {}
    for k,v in cntr_users.iteritems():
        if v >= 10 and k != 'Anonymous':
            ud[k] = v
    
    return ud

def compare(zid,user_name,verbose=True):

    their_answer = consensus.one_answer(zid,user_name)
    expert_answer = consensus.checksum(zid)

    nt = len(their_answer)
    nx = len(expert_answer)

    if nt == nx:
        match_temp = True
        for tk in their_answer.iterkeys():
            if expert_answer.has_key(tk):
                match_temp = True and expert_answer.has_key(tk)

        if verbose:
            if match_temp:
                print "%s's answer for %s (%i sources) matches that of the experts." % (user_name,zid,nt)
            else:
                print "%s's answer for %s (%i sources) does not match that of the experts." % (user_name,zid,nt)

        match_expert = match_temp
    else:
        if verbose:
            print "%s's number of sources for %s does not match (%i vs %i)" % (user_name,zid,nt,nx)
        match_expert = False

    return match_expert

def all_compare(gs,user_name):

    # Run the comparison between a single volunteer and the expert RGZ science team

    matches = []
    for zid in gs:

        sub = subjects.find_one({'zooniverse_id':zid})
        imgid = sub['_id']

        c = classifications.find_one({"subject_ids": imgid, "updated_at": {"$gt": consensus.main_release_date},'user_name':user_name})
        
        if c is not None:
            match_expert = compare(zid,user_name,verbose=False)
            matches.append(match_expert)

    return matches

def ac_data():

    # Compute and save the levels of agreement for volunteers and experts for all overlapping galaxies

    gs = get_galaxies()
    ud = top_volunteers(gs)

    match_ratio = []
    for user_name in ud:
        matches = all_compare(gs,user_name)
        match_ratio.append(np.sum(matches,dtype=float)/len(matches))

    # List of the ratio matches for a single user and the RGZ science team
    with open('%s/goldstandard/gs_compare.pkl' % rgz_dir, 'wb') as output:
        pickle.dump(match_ratio, output)
    
    # Dictionary of the usernames and number of expert100 galaxies classified
    with open('%s/goldstandard/gs_topusers.pkl' % rgz_dir, 'wb') as output:
        pickle.dump(ud, output)
    
    return None

def barchart(expert=False):

    # Load saved data from pickled files

    #file_compare = 'gs_compare.pkl'
    #file_topusers = 'gs_topusers.pkl'
    file_compare = 'gs_compare_ex.pkl'
    file_topusers = 'gs_compare_vol.pkl'

    with open('%s/goldstandard/%s' % (rgz_dir,file_compare), 'rb') as pkl_file:
        data = pickle.load(pkl_file)

    with open('%s/goldstandard/%s' % (rgz_dir,file_topusers), 'rb') as output:
        ud = pickle.load(output)

    ind = np.arange(len(data))  # the x locations for the groups
    width = 0.75                # the width of the bars

    # Parameters for binomial uncertainty bars
    c = 0.683
    n = ud.values()
    k = [int(nn*mr) for nn,mr in zip(n,data)]
    p_lower = [dist.beta.ppf((1-c)/2., kk+1,nn-kk+1) for kk,nn in zip(k,n)]
    p_upper = [dist.beta.ppf(1-(1-c)/2., kk+1,nn-kk+1) for kk,nn in zip(k,n)]
    err_lower = [d-pl for d,pl in zip(data,p_lower)]
    err_upper = [pu-d for d,pu in zip(data,p_upper)]

    # Set up bar chart
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, data, width,yerr=[err_lower,err_upper])
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Percentage matched with experts')
    ax.set_title('Volunteer agreement on expert sample of 100 galaxies')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ud.keys(),rotation=45,ha='right' )
    ax.set_ylim(0,1)
    
    # Plot number of total galaxies above bars
    
    for rect,ngal in zip(rects1,ud.itervalues()):
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%i'%int(ngal),
                ha='center', va='bottom',fontsize=8)
    
    fig.subplots_adjust(bottom=0.2)
    plt.show() 

    return None

def update_gs_subjects(subjects): 

    gs_gals = get_galaxies()

    for gal in gs_gals:
        subjects.update({'zooniverse_id':gal},{'$set':{'goldstandard':True}})

    return None

def plot_gs_volunteers():

    # Plot the results of the gold standard classifications by the users, removing the expert classifications

    gs_gals = get_galaxies()
    for gal in gs_gals:
        print gal
        c0 = consensus.checksum(gal,experts_only=False,excluded=experts)
        consensus.plot_consensus(c0,save_fig = True)

    return None


