import shutil
filename = '/Users/willettk/Desktop/rgz_subjects.csv'

#import bending_angles
#pathdict = bending_angles.make_pathdict()

with open(filename,'r') as f:
    for line in f:
        _id,zooniverse_id,first_id = line.rstrip().split(',')
        try:
            #newpath = '/'.join(pathdict[first_id].split('/')[:-2])
            #shutil.move('/Users/willettk/Downloads/contours/{0:}.json'.format(_id),'{0:}/CONTOURS/{1:}.json'.format(newpath,first_id))
            shutil.move('/Users/willettk/Downloads/contours/{0:}.json'.format(_id),'/Volumes/REISEPASS/rgz/raw_images/ATLAS/CONTOURS/{0:}.json'.format(first_id))
        except IOError as inst:
            print "Error",_id,zooniverse_id,first_id,inst.args
