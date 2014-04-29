import rgz
import matplotlib.pyplot as plt
import numpy as np

def duplicate_class_count(c):

    a = classifications.find({'subject_ids':c})
    dfa = pd.DataFrame(list(a))
    dups = dfa.duplicated(cols='user_ip').sum()

    return dups

subjects,classifications,users = rgz.load_rgz_data()
complete_subjects = subjects.find({'state':'complete'})
cs = pd.DataFrame(list(complete_subjects))

for c in cs['_id']:
    duplist.append(duplicate_class_count(c))

dupgood = [x for x in duplist if x <= 20]
plt.hist(dupgood,bins=len(dupgood))
plt.xlabel('Duplicate user classifications per completed RGZ subject')
plt.ylabel('Count')
plt.xlim(0,20)

plt.show()

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
