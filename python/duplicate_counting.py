import rgz
import matplotlib.pyplot as plt
import numpy as np

'''
Some quick code intended to explore whether users are seeing the same subject
more than once in RGZ.

- KWW, 29 Apr 2014

'''

def duplicate_class_count(c):

    a = classifications.find({'subject_ids':c})
    dfa = pd.DataFrame(list(a))
    dup_ip = dfa.duplicated(cols='user_ip')
    dups = dup_ip.sum()
    if dups > 0:
        timestamps = dfa[dup_ip]['created_at']
    else:
        timestamps=''

    return dups,timestamps

subjects,classifications = rgz.load_rgz_data()
complete_subjects = subjects.find({'state':'complete'})
cs = pd.DataFrame(list(complete_subjects))

# There are 19232 complete subjects; this inefficient loop
# takes an hour or two to run. 

duplist,timestamplist = [],[]
cs_small = cs[:1000]

for c in cs_small['_id']:
    d,t = duplicate_class_count(c)
    duplist.append(d)
    timestamplist.append(t)

# Plot histogram of total number of duplicate users
dupgood,timegood = [d,t for d,t in zip(duplist,timestamplist) if d <= 20]
plt.hist(dupgood,bins=len(dupgood))
plt.xlabel('Duplicate user classifications per completed RGZ subject')
plt.ylabel('Count')
plt.xlim(0,20)

plt.show()


# Checking the timestamps of duplicate classifications

tgood = [t for t in timestamplist if type(t) != str]
tseries = pd.concat(tgood)

thist = pd.Series()
for t in tgood:
    thist = thist.append(pd.Series(np.ones(len(t))*len(t),index=t.tolist()))

thist_sorted = thist.sort_index()
plt.plot(thist_sorted.index,thist_sorted.values)
plt.xlabel('Timestamp of duplicate classification',fontsize=20)
plt.ylabel('Number of duplicate users for this classification',fontsize=20)

# Results from looking at specific list of subjects sent by JeanTate,
# forwarded in email from Ivy Wong. 

'''
ARG0001wxp - JeanTate has classified it twice.
ARG0002ycm - No logged-in user has classified it more than once.
ARG0001oo1 - No logged-in user has classified it more than once.
ARG00019n4 - No logged-in user has classified it more than once.
ARG0000r4u - mpetruck has classified it twice.
ARG0001wxp - JeanTate has classified it twice.
ARG00029to - No logged-in user has classified it more than once.
ARG0000p1g - No logged-in user has classified it more than once.
ARG0003nns - No logged-in user has classified it more than once.
'''
