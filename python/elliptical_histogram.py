import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

plt.ion()

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

# Plot the normalized histogram of (W2-W3) colors for RGZ and WISE all-sky sources

filenames = ['%s/%s.fits' % (rgz_dir,x) for x in ('wise_allsky_2M','gurkan/gurkan_all','rgz_75_wise_16jan')]
labels = ('WISE all-sky sources','Gurkan+14 radio galaxies','RGZ 75% radio galaxies')
    
wise_snr = 5.
 
with fits.open(filenames[0]) as f:
    d = f[1].data
    maglim_w1 = (d['w1snr'] > wise_snr)
    maglim_w2 = (d['w2snr'] > wise_snr)
    maglim_w3 = (d['w3snr'] > wise_snr)
    wise = d[maglim_w1 & maglim_w2 & maglim_w3]
with fits.open(filenames[2]) as f:
    d = f[1].data
    rgz75 = d['ratio'] >= 0.75  
    snr_w1 = (d['snr1'] >= wise_snr)
    snr_w2 = (d['snr2'] >= wise_snr)
    snr_w3 = (d['snr3'] >= wise_snr)
    rgz = d[rgz75 & snr_w1 & snr_w2 & snr_w3]

# Cut both

wise_ell = wise[((wise['w2mpro'] - wise['w3mpro']) > -0.5) & ((wise['w2mpro'] - wise['w3mpro']) < 1.5)]
rgz_ell = rgz[((rgz['w2mpro'] - rgz['w3mpro']) > -0.5) & ((rgz['w2mpro'] - rgz['w3mpro']) < 1.5)]

fig = plt.figure(2,(8,8))
ax1 = fig.add_subplot(111)


from astroML.plotting import hist as histML

c1 = '#e41a1c'
c2 = '#377eb8'
c2 = '#a6cee3'
c3 = '#386cb0'

histML(wise_ell['w2mpro'] - wise_ell['w3mpro'], bins=25, ax=ax1, histtype='stepfilled', alpha=0.5, color=c2, weights=np.zeros(len(wise_ell)) + 1./len(wise_ell), range=(-0.5,1.5),label='WISE all-sky')
histML(rgz_ell['w2mpro'] - rgz_ell['w3mpro'], bins=25, ax=ax1, histtype='step', linewidth=3, color=c1, weights=np.zeros(len(rgz_ell)) + 1./len(rgz_ell), range=(-0.5,1.5),label='RGZ 75%')
ax1.set_ylim(0,0.30)
ax1.vlines(x=np.median(wise_ell['w2mpro'] - wise_ell['w3mpro']),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c3,linestyle='--')
ax1.vlines(x=np.median(rgz_ell['w2mpro'] - rgz_ell['w3mpro']),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c1,linestyle='-')
ax1.set_xlabel(r'$(W2-W3)$',fontsize=24)
ax1.set_ylabel('normalized count',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)

ax1.legend()

#plt.show()

fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/paper/figures/elliptical_histogram.eps')

