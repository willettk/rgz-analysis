# Analyze distribution of RGZ counterparts in WISE color-color space
#

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
from astropy.io import fits
import numpy as np

# WISE All-sky sample
#
filenames = ['%s/%s.fits' % (rgz_dir,x) for x in ('wise_allsky_2M','gurkan/gurkan_all','rgz_75_wise')]
labels = ('WISE all-sky sources','Gurkan+14 radio galaxies','RGZ 75% radio galaxies')

print ''
for fname,label in zip(filenames,labels):
    with fits.open(fname) as f:
        d = f[1].data
    
    if label == 'RGZ 75% radio galaxies':
        d = d[d['ratio']>=0.75]

    w1 = d['w1mpro']
    w2 = d['w2mpro']
    w3 = d['w3mpro']
    w4 = d['w4mpro']

    x = w2-w3
    y = w1-w2
    
    # AGN wedge is INCORRECTLY cited in Gurkan+14; check original Mateos+12 for numbers
    #
    wedge_lims = (y > -3.172*x + 7.624) & (y > (0.315*x - 0.222)) & (y < (0.315*x + 0.796))
    #
    # Very rough loci from Wright et al. (2010)
    stars_lims = (x > 0) & (x < 1) & (y > 0.1) & (y < 0.4)
    el_lims = (x > 0.5) & (x < 1.3) & (y > 0.) & (y < 0.2)
    sp_lims = (x > 1.5) & (x < 3.0) & (y > 0.1) & (y < 0.4)
    
    agn_frac = wedge_lims.sum()/float(len(d))
    stars_frac = stars_lims.sum()/float(len(d))
    el_frac = el_lims.sum()/float(len(d))
    sp_frac = sp_lims.sum()/float(len(d))

    print 'Fraction of %25s in AGN wedge: %4.1f percent' % (label,agn_frac*100)
    print 'Fraction of %25s in stars locus: %4.1f percent' % (label,stars_frac*100)
    print 'Fraction of %25s in elliptical locus: %4.1f percent' % (label,el_frac*100)
    print 'Fraction of %25s in spiral locus: %4.1f percent' % (label,sp_frac*100)
    print ''

print ''

'''
# Make empty arrays for TOPCAT
#
xb,yb = 2.250,0.487
xt,yt = 1.958,1.413

xab = np.linspace(xb,6,100)
xat = np.linspace(xt,6,100)
xal = np.linspace(xb,xt,100)

yab = 0.315*xab - 0.222
yat = 0.315*xat + 0.796
yal =-3.172*xal + 7.624

xall = np.append(xab,np.append(xat,xal))
yall = np.append(yab,np.append(yat,yal))

with open('%s/agn_wedge.txt' % rgz_dir,'w') as f:
    for x,y in zip(xall,yall):
        print >> f,x,y
'''

# Bin data and look at differences?
#

with fits.open(filenames[0]) as f:
    wise = f[1].data
with fits.open(filenames[2]) as f:
    rgz = f[1].data

bins_w2w3 = np.linspace(-1,7,25)
bins_w1w2 = np.linspace(-0.5,3,25)
hw,xedges,yedges = np.histogram2d(wise['w2mpro']-wise['w3mpro'],wise['w1mpro']-wise['w2mpro'],bins=(bins_w2w3,bins_w1w2))
hr,xedges,yedges = np.histogram2d(rgz['w2mpro']-rgz['w3mpro'],rgz['w1mpro']-rgz['w2mpro'],bins=(bins_w2w3,bins_w1w2))

from matplotlib import pyplot as plt
from matplotlib import cm
fig = plt.figure(1,(10,5))
fig.clf()

hw_norm = hw/float(np.max(hw))
hr_norm = hr/float(np.max(hr))

from numpy import ma

hw_norm_masked = ma.masked_array(hw_norm,mask=(hw <= 10))
hr_norm_masked = ma.masked_array(hr_norm,mask=(hr <= 10))

extent = [bins_w2w3[0],bins_w2w3[-1],bins_w1w2[0],bins_w1w2[-1]]

ax1 = fig.add_subplot(121)
cmap = cm.jet
cmap.set_bad('w')
im1 = ax1.imshow(hw_norm_masked.T, alpha=1.0, extent=extent, vmin = 0., vmax = 1., interpolation='nearest', origin='lower')
ax1.set_title('WISE All-Sky')
ax1.set_xlabel('(W2-W3)')
ax1.set_ylabel('(W1-W2)')
ax1.set_aspect('auto')

ax2 = fig.add_subplot(122)
cmap = cm.jet
im2 = ax2.imshow(hr_norm_masked.T, alpha=1.0, extent=extent,vmin = 0., vmax = 1., interpolation='nearest', origin='lower')
ax2.set_title('RGZ 75%')
ax2.set_xlabel('(W2-W3)')
ax2.set_aspect('auto')

position=fig.add_axes([0.92,0.1,0.02,0.80])
cb = plt.colorbar(im2,cax=position,orientation='vertical')
cb.set_label('Normalized ratio',fontsize=16)

'''
ax3 = fig.add_subplot(133)
cmap = cm.jet
im3 = ax3.imshow((np.log10(hr_norm/hw_norm)).T, alpha=1.0, extent=extent,interpolation='nearest', origin='lower')
ax3.set_title('RGZ/WISE ratio')
ax3.set_aspect('auto')

position=fig.add_axes([0.92,0.1,0.02,0.80])
cb = plt.colorbar(im3,cax=position,orientation='vertical')
cb.set_label('log(ratio)',fontsize=16)
'''

#plt.show()

fig.savefig('%s/wise_rgz_fractions.png' % rgz_dir)

