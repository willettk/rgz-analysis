import rgz
import requests
import time
import numpy as np
import cPickle as pickle

zid = 'ARG0001n7u'
subjects,classifications = rgz.load_rgz_data()
sublist = list(subjects.find({'metadata.contour_count':1,'state':'complete'}))      # 26,950 as of the 2014-09-04 data dump

path = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/pkl'

# http://sundog.stsci.edu/first/catalogs/
FIRST_axmin = 5.4 # arcsec
FIRST_area = FIRST_axmin**2 * np.pi

def shoelace(points):
    # Find enclosed area of polygon in Cartesian plane using shoelace formula 
    # https://en.wikipedia.org/wiki/Shoelace_formula
    #
    #points = [[3,4],[5,11],[12,8],[9,5],[5,6]]
    rpoints = np.roll(points,-1,axis=0)
    
    mdet = 0.
    for p,r in zip(points,rpoints):
        mdet += np.linalg.det([p,r])
    
    area = 0.5 * np.abs(mdet)
    
    return area

def get_all_1component_areas():
    '''
    Need additional cut that there is no other radio component within, say, a few arcmin
    '''
    
    n_resolved = 0
    n_compact = 0
    
    resolved_name = []
    resolved_area = []
    compact_name = []
    compact_area = []
    
    print '\nArea of FIRST beam:    %.2f sq arcsec\n' % FIRST_area

    for idx,s in enumerate(sublist):
    
        r = requests.get(s['location']['contours'])
        contours = r.json()
        
        if abs(contours['width'] - contours['height']) > 1:
            print 'Width (%i pixels) and height (%i pixels) of image do not match for %s' % (contours['width'],contours['height'],s['zooniverse_id'])
    
        arcsec_per_pixel = 180./contours['width']
        
        # First entry (containing bounding box) should be outermost contour of only component
    
        outercontour = contours['contours'][0][0]
        assert outercontour['k'] == 0, \
            'Not outermost contour'
    
        points = [(foo['x'],foo['y']) for foo in outercontour['arr']]
        area_pix = shoelace(points)     # Area in sq. pixels
    
        area_sqarcsec = area_pix * arcsec_per_pixel**2
    
        if area_sqarcsec < np.sqrt(2.) * FIRST_area:
            #print '%s is compact  (%.2f sq arcsec)' % (s['zooniverse_id'],area_sqarcsec)
            compact_name.append(s['zooniverse_id'])
            compact_area.append(area_sqarcsec)
        else:
            #print '%s is resolved (%.2f sq arcsec)' % (s['zooniverse_id'],area_sqarcsec)
            resolved_name.append(s['zooniverse_id'])
            resolved_area.append(area_sqarcsec)
    
        if not idx % 1000:
            print idx,time.strftime('%c')
    
    # Save files
    
    pickle.dump(compact_name, open('%s/compact_name.pkl' % path,'wb'))
    pickle.dump(compact_area, open('%s/compact_area.pkl' % path,'wb'))
    pickle.dump(resolved_name,open('%s/resolved_name.pkl' % path,'wb'))
    pickle.dump(resolved_area,open('%s/resolved_area.pkl' % path,'wb'))

    return None

def plot_compact_hist():

    compact_area = pickle.load(open('%s/compact_area.pkl' % path,'rb'))
    resolved_area = pickle.load(open('%s/resolved_area.pkl' % path,'rb'))

    # Plot histogram
    
    import matplotlib.pyplot as plt
    
    nbins = 50
    maxarea = 1500
    fig = plt.figure(1,figsize=(8,8))
    ax = fig.add_subplot(111)
    n,bins,patches = ax.hist(resolved_area+compact_area,bins=nbins,range=(0,maxarea),cumulative=True,label='All 1-component images')
    ax.hist(compact_area,              bins=nbins,range=(0,maxarea),cumulative=True,label='Compact')
    
    ax.vlines(np.sqrt(2)*FIRST_area,ax.get_ylim()[0],ax.get_ylim()[1],linestyles='--',lw=2,color='k')

    # Note that the fraction on the y-axis is wrong unless ymax is reset
    
    ax.set_xlim(0,maxarea)
    ax.set_xlabel('Component area [sq arcsec]',fontsize=16)
    ax.set_ylabel('Count',fontsize=16)
    
    ax2 = ax.twinx()
    ntl = np.arange(6)/5.
    new_tick_locations = [x * max(n)/ax.get_ylim()[1] for x in ntl]
    new_tick_labels = ['%.2f' % x for x in ntl]
    
    ax2.set_yticks(new_tick_locations)
    ax2.set_yticklabels(new_tick_labels)
    ax2.set_ylabel('Fraction',rotation=270,labelpad=15)
    
    legend = ax.legend(loc='upper left', shadow=True)
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('medium')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
        plt.show()
    
    plt.show()
    
    fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/plots/compact_area.pdf')

if __name__ == '__main__':

    get_all_1component_areas()
    plot_compact_hist()
