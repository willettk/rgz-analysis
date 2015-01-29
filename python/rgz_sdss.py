from astropy.io import fits
from astropy.cosmology import WMAP9
from astropy import units as u
from matplotlib import pyplot as plt
import numpy as np

from scipy import stats

plt.ion()

def get_rgz_data():

    with fits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/rgz_wise_75_sdss.fits') as f:
        data = f[1].data

    return data

def get_na10_data():

    with fits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/fits/na10.fits') as f:
        data = f[1].data

    return data

def plot_rgz_sdss(rgzdata,nadata):

    # Plot comparison of the (g-r) CMD for both the RGZ 75% consensus sources (matched with WISE) and SDSS galaxies from NA10

    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)

    gr_color = rgzdata['absmagG'] - rgzdata['absmagR']
    absmagR = rgzdata['absmagR']

    na_gr_color = nadata['g-r']
    dm = WMAP9.distmod(nadata['z'])
    na_absmagR = (nadata['r'] * u.mag - dm).value

    h,xedges,yedges,im=plt.hist2d(absmagR,gr_color,bins=30,range=[[-24,-18],[0,1.5]],cmap=plt.cm.hot)
    #h,xedges,yedges = np.histogram2d(absmagR,gr_color,bins=20,range=[[-24,-18],[0,1.5]])
    xa,ya = np.array(xedges),np.array(yedges)
    xcen = xa[:-1] + (xa[1]-xa[0])/2.
    ycen = ya[:-1] + (ya[1]-ya[0])/2.

    #cs = plt.contour(xcen,ycen,h.T,10,vmin=0.,colors='k',linewidths=1)
    #ax.scatter(absmagR,gr_color,s=1,c='r')

    hn,xedgesn,yedgesn = np.histogram2d(na_absmagR,na_gr_color,bins=30,range=[[-24,-18],[0,1.5]])

    csn = plt.contour(xcen,ycen,hn.T,10,vmin=0.,colors='g',linewidths=2)
    #ax.scatter(na_absmagR,na_gr_color,s=1,c='b')

    #Create custom artists for legend
    c = csn.collections[0]
    art1 = plt.Line2D((0,1),(0,0), color='g')
    
    legend = ax.legend([art1],['Nair & Abraham (2010)'],loc='upper left', shadow=True)

    # Set final properties
    ax.set_xlim(-18,-24)
    ax.set_ylim(0.2,1.0)

    ax.set_xlabel(r'$M_r$',fontsize=20)
    ax.set_ylabel(r'$(g-r)$',fontsize=20)

    plt.show()

    # Try linear fits to the data, see if slope of the red sequence is same

    ind1 = (gr_color > 0.6) & (gr_color < 1.0)
    ind2 = (na_gr_color > 0.6) & (na_gr_color < 1.0)
    x1,y1 = absmagR[ind1],gr_color[ind1]
    x2,y2 = na_absmagR[ind2],na_gr_color[ind2]

    slope1, intercept1, r_value1, p_value1, slope_std_error1 = stats.linregress(x1, y1)
    slope2, intercept2, r_value2, p_value2, slope_std_error2 = stats.linregress(x2, y2)

    print 'RGZ:  slope = %.2f +- %.5f' % (slope1,slope_std_error1)
    print 'NA10: slope = %.2f +- %.5f' % (slope2,slope_std_error2)

    return None

