import rgz
import datetime
import numpy as np
import scipy.optimize

from matplotlib import pyplot as plt
import matplotlib.dates as mdates

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

# When will we finish with the RGZ classifications?

subjects,classifications,users = rgz.load_rgz_data()
jan1 = datetime.datetime(2015, 01, 01, 0, 0, 0, 0)

# How many images are there?

N = subjects.count()

# How many are done?

N_done = subjects.find({'state':'complete'}).count()

# how many are left to do?

N_active = subjects.find({'state':'active'}).count()
N_inactive = subjects.find({'state':'inactive'}).count()

# How many classifications per day over the last three months? Is it constant?

nclass_3months = classifications.find({"updated_at": {"$gt": jan1}}).count()
today = datetime.datetime.today()
ndays = (today - jan1).days

classrate = nclass_3months / float(ndays)

# Plot classification rate for last three months

crarr = []
datearr = []
for d in range(ndays):
    dstart = jan1 + datetime.timedelta(d)
    dend = dstart + datetime.timedelta(1)
    nclass = classifications.find({"updated_at": {"$gte": dstart,"$lt": dend}}).count()
    datearr.append(dstart)
    crarr.append(nclass)

# Plot it


fig = plt.figure(1,(11,6))
ax = fig.add_subplot(111)

ax.plot_date(datearr,crarr)
ax.set_ylabel('Number of classifications')
ax.set_xlabel('Date')

fig.savefig('%s/plots/rate_of_classifications_full.png' % rgz_dir)

ax.set_ylim(0,5000)
fig.savefig('%s/plots/rate_of_classifications_ylim.png' % rgz_dir)

# Plot classifications per day for entire lifetime of RGZ

main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)
all_class = []
all_dates = []
alldays = (today - main_release_date).days
for d in range(alldays):
    dstart = main_release_date + datetime.timedelta(d)
    dend = dstart + datetime.timedelta(1)
    nclass = classifications.find({"updated_at": {"$gte": dstart,"$lt": dend}}).count()
    all_dates.append(dstart)
    all_class.append(nclass)


def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y

def exponential(p,x):
    N0,lamb,x0,b=p
    y = N0 * np.exp(-1*lamb*(x)) + b
    return y

def residuals(p,x,y):
    return y - exponential(p,x)

x = mdates.date2num(all_dates)
xnorm = x - x[0]

p_guess=(10000.,1.,0,100.)
p, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,p_guess,args=(x,all_class),full_output=1)

x0,y0,c,k=p
xp = np.linspace(x[0], x[-1], 1500)
pxp=exponential(p,xp)

fig = plt.figure(2,(15,6))
ax = fig.add_subplot(111)

ax.plot_date(all_dates,all_class)

# Plot fit
ax.plot(xp, pxp, '-')

ax.set_ylabel('Number of classifications')
ax.set_xlabel('Date')

fig.savefig('%s/plots/rate_of_all_classifications_full.png' % rgz_dir)

# Assume a flat rate

rate = np.median(crarr)

# How many classifications are left?

'''
1. active with 1 component
2. active with 2+ components
3. inactive with 1 component
4. inactive with 2+ components
'''

nas = []
active_single = subjects.find({'state':'active','metadata.contour_count':1})
for x in active_single:
    nas.append(5 - x['classification_count'])
nas_total = np.sum(nas)

nam = []
active_multiple = subjects.find({'state':'active','metadata.contour_count':{"$gt":1}})
for x in active_multiple:
    nam.append(20 - x['classification_count'])
nam_total = np.sum(nam)

n_inactive_single = subjects.find({'state':'inactive','metadata.contour_count':1}).count()
nis_total = n_inactive_single * 5

n_inactive_multiple = subjects.find({'state':'inactive','metadata.contour_count':1}).count()
nim_total = n_inactive_multiple * 20

remaining_classifications = nas_total + nam_total + nis_total + nim_total

remaining_time = remaining_classifications / rate

print '%i days left' % remaining_time
print 'Will finish on %s' % (today + datetime.timedelta(int(remaining_time))).strftime('%b %d, %Y')


