# Make a gallery of images showing the RGZ consensus double sources, sorted by bending angle.

from astropy.io import ascii

path = '/Users/willettk/Astronomy/Research/GalaxyZoo'
data = ascii.read('{:}/rgz-analysis/csv/static_catalog3.csv'.format(path),delimiter=' ')

import bending_angles as ba
import numpy as np

pathdict = ba.make_pathdict()

def bending_examples():
    for a in np.linspace(0,80,9):
        bdata = data[(data['bending_angle'] >= a) & (data['bending_angle'] < a+10.)]
        count,errcount = 0,0
        if len(bdata) > 0:
            for b in bdata:
                zid = b['zooniverse_id']
                try:
                    if b['angle_type'] == 'multipeaked_singles':
                        angle_type = 'mps'
                    else:
                        angle_type = 'radio'
                    ba.plot_one_double(zid,pathdict,save_fig=True,anglepath='{0:.0f}_{1:.0f}/'.format(a,a+10),dbltype=angle_type)
                    count += 1
                except ValueError as inst:
                    print "ValueError,",inst.args,zid
                    errcount += 1
        print '{:d} galaxies with bending angle, {:d} with errors for angles between {:.0f} and {:.0f}'.format(count,errcount,a,a+10)
