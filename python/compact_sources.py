import rgz
from astropy.io import ascii,fits
from astropy import wcs
import requests

import numpy as np

from matplotlib import pyplot as plt
from astroML.plotting import hist as histML
    
subjects,classifications,users = rgz.load_rgz_data()

# Find the distribution of areas and/or max sizes 

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
rgz_75_file = '%s/csv/rgz_75.csv' % rgz_dir

def get_singles():

    # Load in all examples of single-component radio galaxies from RGZ

    rgz75 = ascii.read(rgz_75_file,format='csv')

    idx = (rgz75['n_lobes'] == 1)
    singles = rgz75[idx]

    return singles

def get_bbox(galaxy):

    sub = subjects.find_one({'zooniverse_id':galaxy['zid']})

    r = requests.get(sub['location']['contours'])
    try:
        contours = r.json()
    except ValueError:
        return []
    
    radio_components = contours['contours']

    try:
        assert len(radio_components) == 1, \
            'Radio data only %i components for %s' % (len(radio_components),galaxy['zid'])
    except AssertionError:
        return []

    bbox = radio_components[0][0]['bbox']

    return bbox

def area_comp(singles):

    area = []
    for s in singles:

        bbox = get_bbox(s)
        if len(bbox) == 4:
            xmax,ymax,xmin,ymin = bbox
            a = (xmax-xmin) * np.float(ymax-ymin)
            area.append(a)

    return area

def plot_area_hist(area):
    fig = plt.figure(2,(8,8))
    ax1 = fig.add_subplot(111)
    
    c1 = '#e41a1c'
    c2 = '#377eb8'
    c2 = '#a6cee3'
    c3 = '#386cb0'
    
    histML(area, bins=25, ax=ax1, histtype='stepfilled', alpha=1.0, color=c1, range=(0,180),label='RGZ singles')
    ax1.vlines(x=np.median(area),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color='k',linestyle='--')
    ax1.set_xlabel('bbox area [sq. pixels]',fontsize=24)
    ax1.set_ylabel('count',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    
    #ax1.legend()
    #plt.show()
    
    fig.tight_layout()
    fig.savefig('%s/plots/area_hist.pdf' % rgz_dir)

    return None

if __name__ == '__main__':

    singles = get_singles()
    area = area_comp(singles)
    plot_area_hist(area)
