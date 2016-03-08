from matplotlib import pyplot as plt
import numpy as np
import datetime

'''

Data made with following commands:

    within Mongo shell:
    
        > use radio
        > var result = db.radio_classifications.aggregate({$match:{user_id:{$exists:false}}},{$group:{_id:{year:{$year:"$created_at"},month:{$month:"$created_at"},day:{$dayOfMonth:"$created_at"}},count:{$sum:1}}})
        > db.nouserid.insert(result.result)
        > var result = db.radio_classifications.aggregate({$match:{user_id:{$exists:true}}},{$group:{_id:{year:{$year:"$created_at"},month:{$month:"$created_at"},day:{$dayOfMonth:"$created_at"}},count:{$sum:1}}})
        > db.userid.insert(result.result)
    
    command line:
    
        > mongoexport -d radio -c nouserid -f "_id.year,_id.month,_id.day,count" --csv -o ~/Astronomy/Research/GalaxyZoo/rgz-analysis/csv/nouserid.csv
        > mongoexport -d radio -c userid -f "_id.year,_id.month,_id.day,count" --csv -o ~/Astronomy/Research/GalaxyZoo/rgz-analysis/csv/userid.csv

'''

path = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

with open('%s/csv/nouserid.csv' % path,'rb') as f:
    header = f.next()
    dates_nouserid = []
    count_nouserid = []
    for line in f:
        l = line.strip().split(',')
        dates_nouserid.append(datetime.datetime(int(l[0]),int(l[1]),int(l[2])))
        count_nouserid.append(float(l[-1]))

with open('%s/csv/userid.csv' % path,'rb') as f:
    header = f.next()
    dates_userid = []
    count_userid = []
    for line in f:
        l = line.strip().split(',')
        dates_userid.append(datetime.datetime(int(l[0]),int(l[1]),int(l[2])))
        count_userid.append(float(l[-1]))

fig = plt.figure(figsize=(6,6))
'''
ax1 = fig.add_subplot(121)

ax1.plot_date(dates_nouserid,count_nouserid,'-',alpha=1.0)
ax1.plot_date(dates_userid,count_userid,'-',alpha=1.0)
'''

ax2 = fig.add_subplot(111)
ax2.plot_date(dates_nouserid,np.cumsum(count_nouserid),'-',alpha=1.0,label='user_id')
ax2.plot_date(dates_userid,np.cumsum(count_userid),'-',alpha=1.0,label='No user_id')

ax2.legend()

plt.show()
