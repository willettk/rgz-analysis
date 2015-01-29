import numpy as np

randomsize = 2e6
index = 0

ra_lim = [x * 15. for x in (10,15)]
dec_lim = (0,60)

with open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/random_points.csv','w') as f:
    print >> f,'ra,dec'

    while index < randomsize:
        thisindex = index 
        while index == thisindex:
            ra = np.random.random_sample() * 360.
            dec = np.arcsin(2.*np.random.random_sample() - 1.) * (180./ np.pi)
            if (ra >= ra_lim[0]) and (ra <= ra_lim[1]) and (dec >= dec_lim[0]) and (dec <= dec_lim[1]):
                print >> f,'%.4f,%.4f' % (ra,dec)
                index += 1

    
