import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from scipy import stats

data = ascii.read('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/rgz_wise_75_sdss_magerr.csv')

goodg = (data['g'] < 22.2) & (data['err_g'] < 0.05)
goodr = (data['r'] < 22.2) & (data['err_r'] < 0.05)

print len(data)

data = data[goodg & goodr]

print len(data)

'''
# split by redshift
#
nbins = 20
zbins = np.arange(nbins+1)/float(nbins)
deltaz = zbins[1] - zbins[0]

fig = plt.figure(1,(15,10))

# Try fitting the red sequence at z=0.1
#

ind = (data['redshift'] >= 0.05) & (data['redshift'] < 0.10) & ((data['absmagG'] - data['absmagR']) > 0.7)
x1,y1 = data[ind]['absmagR'], data[ind]['absmagG'] - data[ind]['absmagR']
slope1, intercept1, r_value1, p_value1, slope_std_error1 = stats.linregress(x1, y1)
xarr = np.linspace(-24,-19,100)

for idx,zb in enumerate(zbins[:-1]):
    ax = fig.add_subplot(4,5,idx+1)
    ind = (data['redshift'] >= zb) & (data['redshift'] < zb+deltaz)
    ax.scatter(data[ind]['absmagR'], data[ind]['absmagG'] - data[ind]['absmagR'])
    ax.plot(xarr,intercept1 + xarr*slope1, linestyle='--',color='k')
    ax.set_xlim(-19,-24)
    ax.set_ylim(0.3,1.1)
    ax.set_title(r'$z=$%s' % zb,fontsize=12)

    if not idx % 5:
        ax.set_ylabel(r'$(g-r)$',fontsize=16)
    else:
        ax.set_yticklabels(('','','','','',''))
    if idx > 14:
        ax.set_xlabel(r'$M_R$',fontsize=16)
    else:
        ax.set_xticklabels(('','','','','','','','',''))

plt.show()

'''
fig2 = plt.figure(2)

ax = fig2.add_subplot(111)

ax.hist(data['redshift'],bins=20)
ax.set_ylabel('Count')
ax.set_xlabel(r'$z$')
ax.set_title('Redshift distribution of RGZ-SDSS optical matches')

plt.show()
