import consensus
import bending_angles as ba
import rgz

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

path = '/Users/willettk/Astronomy/Research/GalaxyZoo'

def load_rgz_data():

    client = rgz.MongoClient('localhost', 27017)
    db = client['radio'] 
    
    subjects = db['radio_subjects']
    classifications = db['radio_classifications']
    catalog = db['catalog']

    return subjects,classifications,catalog

# Analyze the results of the consensus algorithm and catalog aggregation

def consensus_hist(catalog,savefig=False):

    # Make histogram of the consensus level for sources

    df = pd.read_csv('%s/rgz-analysis/csv/static_catalog_full.csv' % path,delim_whitespace=True)

    fig = plt.figure(figsize=(15,8))

    all_entries = catalog.find()
    consensus_level = np.array([])
    n_components = []
    
    for a in all_entries:
        consensus_level = np.append(consensus_level,a['consensus']['level'])
        n_components = np.append(n_components,a['radio']['numberComponents'])

    ax1 = fig.add_subplot(121)
    ax1.hist(consensus_level,bins=20)

    ax1.set_yscale('log')
    ax1.set_xlabel('Consensus level',fontsize=16)
    ax1.set_ylabel('Count',fontsize=16)
    ax1.set_title('All RGZ completed subjects')

    ax2 = fig.add_subplot(122)
    for i in np.arange(4):
        ax2.hist(consensus_level[n_components == i+1], bins=20, alpha=0.5, label=r'$N=$%i' % (i+1))
    ax2.hist(consensus_level[n_components >= 5], bins=20, alpha=0.5, label=r'$N\geq$%i' % 5)

    ax2.set_yscale('log')

    ax2.set_xlabel('Consensus level',fontsize=16)
    ax2.set_ylabel('Count',fontsize=16)
    ax2.set_title('Grouped by number of radio components')
    ax2.legend(loc='upper left')

    '''
    ax1 = fig.add_subplot(121)
    df['consensus_level'].hist(kind='hist',ax=ax1,bins=20)

    ax2 = fig.add_subplot(122)
    df['consensus_level'].hist(kind='hist',by=df['ax=ax2,bins=20)
    '''

    if savefig:
        fig.savefig('%s/rgz-analysis/plots/analyze_consensus1.pdf' % path)
    else:
        plt.show()

    return None

def edge_down(catalog,savefig=False):

    # Make histogram of the consensus level for sources

    df = pd.read_csv('%s/rgz-analysis/csv/static_catalog_full.csv' % path,delim_whitespace=True)

    fig = plt.figure(figsize=(15,8))

    all_entries = catalog.find()
    consensus_level = np.array([])
    n_components = []
    
    for a in all_entries:
        consensus_level = np.append(consensus_level,a['consensus']['level'])
        n_components = np.append(n_components,a['radio']['numberComponents'])

    ax1 = fig.add_subplot(121)
    n,bins,patches = ax1.hist(consensus_level,bins=20,cumulative=True,histtype='step')
    ax1.set_yscale('log')

    for value in (0.50,0.75):
        ax1.axvline(value,ls='--',color='k')
        ax1.text(value,ax1.get_ylim()[1]*0.95,int(n[(np.abs(bins-value)).argmin()]),fontsize=10,va='top')

    ax1.set_xlabel('Consensus level',fontsize=16)
    ax1.set_ylabel('Count',fontsize=16)
    ax1.set_title('All RGZ completed subjects')

    ax2 = fig.add_subplot(122)
    for i in np.arange(4):
        n,bins,patches = ax2.hist(consensus_level[n_components == i+1], bins=20, alpha=0.5, label=r'$N=$%i' % (i+1),cumulative=True,histtype='step')
        for value in (0.50,0.75):
            print i+1,value,int(n[(np.abs(bins-value)).argmin()-1])

    n,bins,patches = ax2.hist(consensus_level[n_components >= 5], bins=20, alpha=0.5, label=r'$N\geq$%i' % 5,cumulative=True,histtype='step')
    for value in (0.50,0.75):
        ax2.axvline(value,ls='--',color='k')
        print "5+",value,int(n[(np.abs(bins-value)).argmin()-1])


    ax2.set_yscale('log')
    ax2.set_xlabel('Consensus level',fontsize=16)
    ax2.set_ylabel('Count',fontsize=16)
    ax2.set_title('Grouped by number of radio components')
    ax2.legend(loc='upper left')

    if savefig:
        fig.savefig('%s/rgz-analysis/plots/analyze_consensus2.pdf' % path)
    else:
        plt.show()

    return None

