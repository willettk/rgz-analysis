from consensus import db

consensus = db['consensus_multi_ir']
catalog = db['catalog_dr1']

count = 0
with open('ir_clicks.csv', 'w') as f:
    print >> f, 'catalog_id,zooniverse_id,ir_ra,ir_dec,ir_level,radio_level,n_total,n_radio,n_ir,clicks'
    for con in consensus.find({'first_id':{'$exists':True}}):
        count += 1
        if not count%1000:
            print count
        cat = catalog.find_one({'zooniverse_id':con['zooniverse_id'], 'consensus.label':con['label']})
        output = '{},{},{},{},{},{},{},{},{}'.format(cat['catalog_id'], cat['zooniverse_id'], cat['consensus']['ir_ra'] if 'ir_ra' in cat['consensus'] else -99, \
                                                     cat['consensus']['ir_dec'] if 'ir_dec' in cat['consensus'] else -99, cat['consensus']['ir_level'], \
                                                     cat['consensus']['radio_level'], cat['consensus']['n_total'], cat['consensus']['n_radio'], \
                                                     cat['consensus']['n_ir'])
        for i in con['ir_count']:
            print >> f, '{},{}'.format(output, i)
