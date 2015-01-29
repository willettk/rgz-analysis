'''
In my script, this is accomplished by calculating the Median
Absolute Difference of an image and converting it to the standard
deviation by dividing by a certain factor. Using the MAD is much more
robust against the presence of strong sources than simply calculating
the standard deviation directly.
'''

from scipy.special import erfinv
import numpy as np
from astropy.io import fits

# Need image data - try downloading ELAIS/CDFS fields and see if they're there. If not, bug Ed or Ivy.

imgdir = '/Volumes/3TB/rgz/raw_images/ATLAS/2x2'
galname = 'CI0001C1'

radio_fits = '%s/%s_radio.fits' % (imgdir,galname)
with fits.open(radio_fits) as f:
    imgdata = f[0].data

# divide by this to convert mad into sigma in radio image noise
# calculation
mad2sigma=np.sqrt(2)*erfinv(2*0.75-1)

# the following is needed to derive some statistics when the input radio image is not a SNR image
med=np.median(imgdata)
mad=np.median(np.abs(imgdata-med))
sigma=mad/mad2sigma

import make_contours
from PIL import Image
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.path import Path
import matplotlib.patches as patches

cs = make_contours.contour(radio_fits, sigma)
cs['contours'] = map(points_to_dict, cs['contours'])
#print json.dumps(cs)

# Plot the infrared results

fig = plt.figure(1,(12,4))

# Download contour data

contours = cs

sf_x = 500./contours['width']
sf_y = 500./contours['height']

verts_all = []
codes_all = []
components = contours['contours']

for comp in components:

    for idx,level in enumerate(comp):
        verts = [((p['x'])*sf_x,(p['y']-1)*sf_y) for p in level['arr']]
        
        codes = np.ones(len(verts),int) * Path.LINETO
        codes[0] = Path.MOVETO
    
        verts_all.extend(verts)
        codes_all.extend(codes)

path = Path(verts_all, codes_all)
patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)

# Scaling factor for FITS to radio files

radio_ir_scaling_factor = 500./132

# Display IR and radio images

im_standard = Image.open('%s/%s_heatmap.png' % (imgdir,galname))
ax1 = fig.add_subplot(131)
ax1.imshow(im_standard,origin='upper')
ax1.add_patch(patch_black)
ax1.set_title('IR image')

im_radio = Image.open('%s/%s_heatmap+contours.png' % (imgdir,galname))
ax2 = fig.add_subplot(132)
ax2.imshow(im_radio,origin='upper')
ax2.set_title('Original levels')

ax4.set_xlim([0, 500])
ax4.set_ylim([500, 0])
ax4.set_title('My contour levels')

ax4.set_aspect('equal')

ax4.add_patch(patch_black)
ax4.yaxis.tick_right()

ax1.get_xaxis().set_ticks([0,100,200,300,400])
ax2.get_xaxis().set_ticks([0,100,200,300,400])
ax3.get_xaxis().set_ticks([0,100,200,300,400])

plt.subplots_adjust(wspace=0.02)

plt.show()
  

# Save hard copy of the figure
#fig.savefig('%s/%s' % (plot_dir,writefile))


