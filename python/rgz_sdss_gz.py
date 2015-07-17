import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from scipy import stats
from astropy.io.votable import parse_single_table
from astropy.cosmology import WMAP9

rgzdir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

def votable_read(fname):

    table = parse_single_table(fname)
    data = table.array

    return data

def na10_read():

    filename = '/Users/willettk/Astronomy/Research/GalaxyZoo/fits/na10.fits'

    with fits.open(filename) as f:
        data = f[1].data

    return data

# Plot the CMD as a function of GZ1 morphology for the RGZ 75% counterparts

data = ascii.read('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/rgz_wise_75_sdss_magerr.csv')

fname_spec = 'rgz_75_sdss_mags_spec.vot'
fname_nospec = 'rgz_75_sdss_mags_nospec.vot'

zoospec = votable_read('%s/%s' % (rgzdir,fname_spec))
zoonospec = votable_read('%s/%s' % (rgzdir,fname_nospec))

na10 = na10_read()

gr_na10 = na10['g-r']
distmod_na10 = WMAP9.distmod(na10['z'])
Mr_na10 = na10['r'] - distmod_na10.value
Mr_na10 = na10['M_g'] - gr_na10

# Generate contours of the NA10 sample

magbins = np.linspace(-24,-20,20)
colorbins = np.linspace(0.3,1.0,20)
hc,xc,yc = np.histogram2d(Mr_na10,gr_na10,bins=(magbins,colorbins))
levels=np.linspace(0,200,10)

# Make 2x2 figure

fig,axarr = plt.subplots(2,2,figsize=(10,10))

from scipy.ndimage.filters import gaussian_filter

spec_goodmags = (zoospec['modelmag_G'] < 22.2) & (zoospec['modelmag_R'] < 22.2) & (zoospec['modelmagerr_G'] < 0.10) & (zoospec['modelmagerr_R'] < 0.10)
sp_spec = (zoospec['spiral'] == 1) & spec_goodmags
el_spec = (zoospec['elliptical'] == 1) & spec_goodmags
un_spec = (zoospec['uncertain'] == 1) & (zoospec['p_mg'] < 0.4) & spec_goodmags
mg_spec = (zoospec['p_mg'] >= 0.4) & spec_goodmags

nospec_goodmags = (zoonospec['modelmag_G'] < 22.2) & (zoonospec['modelmag_R'] < 22.2) & (zoonospec['modelmagerr_G'] < 0.10) & (zoonospec['modelmagerr_R'] < 0.10)
sp_nospec = (zoonospec['p_cs'] >= 0.8) & nospec_goodmags
el_nospec = (zoonospec['p_el'] >= 0.8) & nospec_goodmags
un_nospec = (zoonospec['p_el'] < 0.8) & (zoonospec['p_cs'] < 0.8) & (zoonospec['p_mg'] < 0.4) & nospec_goodmags
mg_nospec = (zoonospec['p_mg'] >= 0.4) & nospec_goodmags

sp_all_mag = np.append(zoospec[sp_spec]['absmagR'].data,zoonospec[sp_nospec]['absmagR'].data)
el_all_mag = np.append(zoospec[el_spec]['absmagR'].data,zoonospec[el_nospec]['absmagR'].data)
un_all_mag = np.append(zoospec[un_spec]['absmagR'].data,zoonospec[un_nospec]['absmagR'].data)
mg_all_mag = np.append(zoospec[mg_spec]['absmagR'].data,zoonospec[mg_nospec]['absmagR'].data)

sp_all_color = np.append(zoospec[sp_spec]['absmagG'].data - zoospec[sp_spec]['absmagR'].data,zoonospec[sp_nospec]['absmagG'].data - zoonospec[sp_nospec]['absmagR'].data)
el_all_color = np.append(zoospec[el_spec]['absmagG'].data - zoospec[el_spec]['absmagR'].data,zoonospec[el_nospec]['absmagG'].data - zoonospec[el_nospec]['absmagR'].data)
un_all_color = np.append(zoospec[un_spec]['absmagG'].data - zoospec[un_spec]['absmagR'].data,zoonospec[un_nospec]['absmagG'].data - zoonospec[un_nospec]['absmagR'].data)
mg_all_color = np.append(zoospec[mg_spec]['absmagG'].data - zoospec[mg_spec]['absmagR'].data,zoonospec[mg_nospec]['absmagG'].data - zoonospec[mg_nospec]['absmagR'].data)

mags = [sp_all_mag,el_all_mag,un_all_mag,mg_all_mag]
colors = [sp_all_color,el_all_color,un_all_color,mg_all_color]
labels = ('Spiral','Elliptical','Uncertain','Merger')
pcolors = ('blue','red','purple','green')

for ax,mag,color,label,pcolor in zip(axarr.reshape(-1),mags,colors,labels,pcolors):

    fi = gaussian_filter(hc.T,0.5)
    CS = ax.contour(magbins[1:],colorbins[1:],fi,levels,colors='k',linewidths=1)
    ax.scatter(mag,color,color=pcolor)
    ax.set_xlim(-24,-20)
    ax.set_ylim(0.3,1.0)
    ax.set_xlabel(r'$M_R$',fontsize=14)
    ax.set_ylabel(r'$^{0.0}(g-r)$',fontsize=14)
    ax.text(-21.0,0.9,r'$N=$%i' % len(mag),ha='left')
    ax.set_xticks([-24,-23,-22,-21,-20])
    ax.set_title(label,fontsize=16)

#plt.show()

fig.savefig('%s/plots/rgz_sdss_gz.png' % rgzdir)

# Print spirals

with open('%s/rgz_spirals.csv' % rgzdir,'wb') as f:
    for w,x,y,z in zip(zoospec[sp_spec]['name'].data,zoospec[sp_spec]['ra'].data,zoospec[sp_spec]['dec'].data,zoospec[sp_spec]['optical_objId'].data):
        f.write('%s, %.3f, %.3f, %i\n' % (w,x,y,z))
        
    for w,x,y,z in zip(zoonospec[sp_nospec]['name'].data,zoonospec[sp_nospec]['ra'].data,zoonospec[sp_nospec]['dec'].data,zoonospec[sp_nospec]['optical_objId'].data):
        f.write('%s, %.3f, %.3f, %i\n' % (w,x,y,z))
    
    
