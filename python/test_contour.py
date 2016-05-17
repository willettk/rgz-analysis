import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import json
import numpy as np
from pymongo import MongoClient
import requests
from StringIO import StringIO
from PIL import Image

zid = 'ARG0001n7u'

client = MongoClient('localhost', 27017)
db = client['radio'] 
subjects = db['radio_subjects'] 		# subjects = images

s = subjects.find_one({'zooniverse_id':zid})
r = requests.get(s['location']['contours'])
contours = r.json()

sf_x = 500./contours['width']
sf_y = 500./contours['height']

verts_all = []
codes_all = []

components = contours['contours']

for comp in components:

    for idx,level in enumerate(comp):
        verts = [((p['x'])*sf_x,(p['y']-1)*sf_y) for p in level['arr']]
        
        #verts.append(verts[0])  # Place beginning contour point at end again to close path?

        codes = np.ones(len(verts),int) * Path.LINETO
        codes[0] = Path.MOVETO
        #codes[-1] = Path.CLOSEPOLY

        '''
        for v,c in zip(verts,codes):
            print idx,v,c
        '''
    
        verts_all.extend(verts)
        codes_all.extend(codes)

path = Path(verts_all, codes_all)
patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)
patch_orange = patches.PathPatch(path, facecolor = 'none', edgecolor='orange', lw=1)

fig = plt.figure()

url_standard = s['location']['standard']
r_standard = requests.get(url_standard)
im_standard = Image.open(StringIO(r_standard.content))
ax1 = fig.add_subplot(121,aspect='equal')
ax1.imshow(im_standard,origin='upper')
ax1.add_patch(patch_black)
ax1.set_title('WISE')

url_radio = s['location']['radio']
r_radio = requests.get(url_radio)
im_radio = Image.open(StringIO(r_radio.content))
ax2 = fig.add_subplot(122,aspect='equal')
ax2.imshow(im_radio,origin='upper')
ax2.add_patch(patch_orange)
ax2.set_title('FIRST')

'''
ax2.set_xlim(220,280)
ax2.set_ylim(280,220)
'''
plt.show()
