# When will we finish with the RGZ classifications?
#
# Kyle Willett, 3 Mar 2015

import rgz
import datetime
import numpy as np
import scipy.optimize

from matplotlib import pyplot as plt
import matplotlib.dates as mdates

# Global variables and paths

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)
today = datetime.datetime.today()

# Mathematical functions for fitting time-series data

def exponential(p,x):
    N0,lamb,b=p
    y = N0 * np.exp(-1*lamb*x) + b
    return y

def residuals(p,x,y):
    return y - exponential(p,x)

def count_classifications(subjects,classifications):

    # Takes several minutes to run on the 15 month RGZ data

    # Calculate the rate at which classifications are done
    
    all_class = []
    all_dates = []
    alldays = (today - main_release_date).days
    for d in range(alldays):
        dstart = main_release_date + datetime.timedelta(d)
        dend = dstart + datetime.timedelta(1)
        nclass = classifications.find({"updated_at": {"$gte": dstart,"$lt": dend}}).count()
        all_dates.append(dstart)
        all_class.append(nclass)

    return all_dates,all_class
    
def how_many_left(subjects,classifications,verbose=False):

    rate = find_rate(classifications)

    # How many classifications are left?
    
    '''
    Four types of subjects in the system:

    - active with 1 component
    - active with 2+ components
    - inactive with 1 component
    - inactive with 2+ components
    '''

    nc_single = 5
    nc_multiple = 20
    
    nas = []
    active_single = subjects.find({'state':'active','metadata.contour_count':1})
    for x in active_single:
        nas.append(nc_single - x['classification_count'])
    nas_total = np.sum(nas)
    
    nam = []
    active_multiple = subjects.find({'state':'active','metadata.contour_count':{"$gt":1}})
    for x in active_multiple:
        nam.append(nc_multiple - x['classification_count'])
    nam_total = np.sum(nam)
    
    n_inactive_single = subjects.find({'state':'inactive','metadata.contour_count':1}).count()
    nis_total = n_inactive_single * nc_single
    
    n_inactive_multiple = subjects.find({'state':'inactive','metadata.contour_count':1}).count()
    nim_total = n_inactive_multiple * nc_multiple
    
    remaining_classifications = nas_total + nam_total + nis_total + nim_total
    finished_classifications = classifications.count()
    
    remaining_time = remaining_classifications / rate
    
    if verbose:
        print '%i classifications to do' % remaining_classifications
        print '%i classifications finished' % finished_classifications
        print '%i days left' % remaining_time
        print 'Will finish on %s' % (today + datetime.timedelta(int(remaining_time))).strftime('%b %d, %Y')

    return remaining_classifications,finished_classifications

def find_rate(classifications,rdays=90):

    # Find the rate of classifications per day, assuming a constant rate over a recent span of time. Default is last three months (90 days).

    c3 = classifications.find({"updated_at": {"$gte": today - datetime.timedelta(rdays),"$lt": today}}).count()
    rate = c3 / rdays

    return rate

def plot_time_left(all_dates,all_class,remaining_classifications,finished_classifications,rate,savefig=False):

    # Plot results
    
    p_guess=(10000.,1.,100.)
    xnorm = mdates.date2num(all_dates) - mdates.date2num(all_dates[0]) + 1
    p, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,p_guess,args=(xnorm,all_class),full_output=1)
    
    xp = np.linspace(xnorm[0], xnorm[-1], 1000)
    plot_dates = mdates.num2date(xp+mdates.date2num(all_dates[0]))
    
    fig = plt.figure(2,(14,12))

    # Panel 1: plot number of classifications per day over lifetime of project

    ax1 = fig.add_subplot(211)
    
    ax1.plot_date(all_dates,all_class)
    ax1.plot(plot_dates, exponential(p,xp), '-')
    
    ax1.set_ylabel('Number of RGZ classifications',fontsize=16)
    ax1.set_xlabel('Date',fontsize=16)
    ax1.set_ylim(1,4.5e4)
    #ax1.set_yscale('log')
    
    # Panel 2: plot cumulative number of classifications

    ax2 = fig.add_subplot(212)

    fc = remaining_classifications + finished_classifications
    cumclass = [np.sum(all_class[:i]) for i in range(len(all_class))]
    ax2.plot_date(all_dates,cumclass,'-',color='orange',linewidth=2)
    ax2.hlines(fc,all_dates[0],datetime.datetime(2019,1,1),linestyle='--',color='k')
    
    ndays = 2000
    y_remaining = finished_classifications + np.arange(ndays) * rate
    today = datetime.datetime.today()
    dates_remaining = [today + datetime.timedelta(i) for i in range(ndays)]
    ax2.plot_date(dates_remaining,y_remaining,'-.',color='gray',linewidth=1)
    
    hiplotlim = (np.round(float(fc)/(10**int(np.log10(fc))),1)+0.2) * 10**int(np.log10(fc))
    ax2.set_xlim(all_dates[0],datetime.datetime(2019,1,1))
    ax2.set_ylim(1,hiplotlim)
    ax2.set_ylabel('Cumulative RGZ classifications',fontsize=16)
    ax2.set_xlabel('Date',fontsize=16)

    if savefig:
        fig.savefig('%s/plots/enddate.pdf' % rgz_dir)
    else:
        plt.show()

    return None

if __name__ == '__main__':

    # Load the data
    subjects,classifications = rgz.load_rgz_data()

    # Find out how many classifications per day the project averages
    rate = find_rate(classifications)

    # Compute how many classifications are left to do, given the size of the active dataset
    remaining_classifications,finished_classifications = how_many_left(subjects,classifications)

    # Find the number of classifications per day over the lifetime of the project
    all_dates,all_class = count_classifications(subjects,classifications)

    # Load the data
    plot_time_left(all_dates,all_class,remaining_classifications,finished_classifications,rate,savefig=True)
