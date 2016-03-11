# Does size matter?

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
import rgz

catalog = rgz.load_catalog()

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()

sizearr = []
consarr = []
comparr = []

for nc in (1,2,3,4,5):
    for c in catalog.find({'radio.numberComponents':nc}):
        components = c['radio']['components']
        minang = components[0]['angularExtent']
        for comp in components:
            if comp['angularExtent'] < minang:
                minang = comp['angularExtent']
        sizearr.append(minang)
        consarr.append(float(c['consensus']['level']))
        comparr.append(nc)

# Convert to pandas dataframe for use with seaborn

import pandas as pd
import numpy as np

d = {'min_angular_size':sizearr,'consensus_level':consarr,'numberComponents':comparr,'logsize':np.log10(sizearr)}
df = pd.DataFrame(d)

g = sns.lmplot(x='logsize',y='consensus_level',data=df,y_jitter=0.025,hue='numberComponents',col='numberComponents',col_wrap=3,size=4,truncate=True,markers='.')

g.set_axis_labels("log(Minimum angular size)", "RGZ consensus level")

for ax,nc in zip(g.axes.ravel(),(1,2,3,4,5)):
    ax.set_xlim(-2.0,0.5)
    ax.set_ylim(-0.1,1.1)
    a = df[df['numberComponents']==nc]
    rho = a['logsize'].corr(a['consensus_level'], method='spearman')
    tau = a['logsize'].corr(a['consensus_level'], method='kendall')
    ax.text(-1.9,0.90,r'$\rho=$%.3f' % rho)
    ax.text(-1.9,0.82,r'$\tau=$%.3f' % tau)

g.savefig('%s/plots/size_consensus.png' % rgz_dir)
