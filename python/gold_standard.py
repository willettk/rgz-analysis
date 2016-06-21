# How many galaxies from the gold standard sample are necessary to establish whether a user agrees with the experts?

# Hm - not many have done more than 1. Might have to expand the experiment to the full expert data set of 100 galaxies ...

# Plot (N/2 - M/2) vs. (M+N) and see where agreement lies!!

# Read in the list of gold standard images

import consensus
import rgz
import operator
import numpy as np
import cPickle as pickle
from matplotlib import pyplot as plt
import scipy.stats.distributions as dist
from collections import Counter

from mpl_defaults import red,blue,green

rgz_dir = consensus.rgz_path
subjects,classifications = rgz.load_rgz_data()
experts=("42jkb", "ivywong", "stasmanian", "klmasters", "Kevin", "akapinska", "enno.middelberg", "xDocR", "DocR", "vrooje", "KWillett")
    
def get_galaxies(expert=False):
    
    if expert:
        filename = "expert/expert_all_zooniverse_ids.txt"
    else:
        filename = "goldstandard/gs_zids.txt"
    
    with open('{0:}/{1:}'.format(rgz_dir,filename)) as f:
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
        consensus.plot_consensus(c0,savefig = True)

    return None

def plot_gs_classified(savefig=True):

    # How many users have completed all galaxies in the gold standard set of 20 subjects? (raised by Ivy during Sep 2015 telecon)

    counts,labels = [],[]
    for gal in subjects.find({'goldstandard':True}):
        counts.append(gal['classification_count'])
        labels.append(gal['zooniverse_id'])

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    N = len(counts)
    rects1 = ax.bar(np.arange(N),counts,color=blue)
    ax.set_xticks(np.arange(N))
    ax.set_xticklabels(labels,rotation=45,fontsize=8)
    ax.set_ylabel(r'$N_\mathrm{users}$',fontsize=20)
    ax.set_title("Number of RGZ users who classified each gold standard subject")

    if savefig:
        plt.savefig('{0}/plots/gs_classified.png'.format(consensus.rgz_path))
    else:
        plt.show()

    return None

def plot_agreement(savefig=True):

    import pandas as pd
    from matplotlib.colors import LogNorm
    
    data = pd.read_csv("{0}/csv/user_weights.csv".format(consensus.rgz_path))

    agreement = np.nan_to_num(data['agreed'] * 1./data['gs_seen'])

    fig,ax = plt.subplots(1,1)
    h = ax.hist2d(data['gs_seen'],agreement,bins=25,range=[[0,50],[0,1]],cmap=plt.cm.viridis,norm=LogNorm())
    fig.colorbar(h[3],ax=ax)

    ax.set_xlabel('Number of gold standard subjects seen')
    ax.set_ylabel(r'$N_\mathrm{agree}/N_\mathrm{seen}$',fontsize=20)

    if savefig:
        plt.savefig('{0}/plots/gs_agreement.png'.format(consensus.rgz_path))
    else:
        plt.show()

    return None

def weights():

    import pandas as pd
    data = pd.read_csv("{0}/csv/user_weights.csv".format(consensus.rgz_path))
    agreement = np.nan_to_num(data['agreed'] * 1./data['gs_seen'])
    weight_val = agreement * 10.

    subindex = classifications.create_index([('user_name',consensus.pymongo.ASCENDING)],name='subject_ids_idx')

    # What do I want to know?

    # Per Larry's suggestion:
    # As a function of N (minimum number of gold standard subjects seen), what's the distribution of user weights?

    fig,axarr = plt.subplots(2,2,sharex=True,sharey=True,figsize=(8,8))
    Narr = (3,5,10,20)

    for N,ax in zip(Narr,axarr.ravel()):
        wvlim = weight_val[np.array(data['gs_seen'] >= N)]
        ax.hist(wvlim,bins=10,histtype='step',lw=3,color=blue)

        ax.text(0,550,r'$N_\mathrm{wtd. users}=$'+'{0} ({1:.0f}%)'.format(len(wvlim),len(wvlim)*100./len(data)),ha='left')

        # What are the number of weighted classifications?

        wc = 0
        for u in data[data['gs_seen'] >= N]['user_name']:
            wc += classifications.find({'user_name':u}).count()
        tc = classifications.count()

        ax.text(0,500,r'$N_\mathrm{wtd. class}=$'+'{0:.0f}%'.format(wc * 100./tc),ha='left')
        ax.text(0,450,r'$\langle w\rangle=$'+'{0:.1f}'.format(np.mean(wvlim),ha='left',color='red'))

        ax.vlines(np.mean(wvlim),ax.get_ylim()[0],ax.get_ylim()[1],color='red',linestyles='--',lw=1)
        ax.set_xlim(0,10)
        ax.set_xlabel('Weight')
        ax.set_ylabel(r'$N_\mathrm{users}$',fontsize=20)
        ax.set_title(r'$N_\mathrm{seen}\geq$'+'{0}'.format(N))

    fig.tight_layout()

    plt.show()
