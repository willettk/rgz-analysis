path = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
m = open('%s/rgz1000_first_zooniverse.txt' % path,'wb')

with open('%s/rgz1000_firstIDs.txt' % path,'rb') as f:
    for firstid in f.readlines():
        s = subjects.find({'metadata.source':firstid.rsplit()[0]})
        ss = s.next()
        print >> m,'%20s, %10s' % (firstid.rsplit()[0], ss['zooniverse_id'])

m.close()
