import json
import datetime

# Running a query on Ourobors - search Talk for particular tags/keywords
# within a specific date range. 

# Actual Talk query is done by a Ruby script (talk_search.rb) - this Python program
# just limits the data range and extracts the Zooniverse ID to construct a URL list.

# Parse the entries by date

start_date = datetime.datetime(2014,8,10)
end_date = datetime.datetime(2014,12,16)

datapath = 'talk_searches'

# Example keywords to search on
rgz_types = ('giant','large','huge','kpc','mpc','overedge')

for r in rgz_types:
    d = json.load(open('%s/radio_%s.json' % (datapath,r),'r'))

    keep = []
    w = open('%s/url_lists/url_%s.txt' % (datapath,r),'w')
    for gal in d:
        updated = datetime.datetime.strptime(gal['updated_at'],'%Y-%m-%dT%H:%M:%SZ')

        if (updated >= start_date) & (updated <= end_date):
            keep.append(gal)
            print >> w,'http://radiotalk.galaxyzoo.org/#/subjects/%s' % gal['name']

    w.close()

    with open('%s/date_limited/dated_%s.json' % (datapath,r),'w') as f:
        json.dump(keep,f)
