from astropy import coordinates as coord, units as u
import numpy as np
import itertools, os, csv
import matplotlib.pyplot as plt
from consensus import db, rgz_path

catalog_0 = db['catalog']
catalog_5 = db['catalog_5']
catalog_10 = db['catalog_10']
catalog_20 = db['catalog_20']
catalogs = {'catalog_0':catalog_0, 'catalog_5':catalog_5, 'catalog_10':catalog_10, 'catalog_20':catalog_20}

def cids_for_removal():
    
    # Find duplicate sources (so ones not in the static catalog)
    
    cids = {}
    
    for weight in [0, 5, 10, 20]:
        cids['cids_{}'.format(weight)] = []
        for c in catalogs['catalog_{}'.format(weight)].find({'duplicate_sources.exact_duplicate':{'$exists':True}}):
            if c['catalog_id'] != min(c['duplicate_sources']['exact_duplicate']):
                cids['cids_{}'.format(weight)].append(c['catalog_id'])
    
    return cids

######### Begin IR block #########

def get_ir_mismatch(cids_unwt, cids_wt, weight):
    
    radio_pos_unwt, radio_pos_wt = set(), set()
    
    # Get set of all radio positions in each catalog
    print '    Finding radio positions'
    for c in catalog_0.find({'catalog_id':{'$nin':cids_unwt}}):
        radio_pos_unwt.add(c['rgz_name'])
    for c in catalogs['catalog_{}'.format(weight)].find({'catalog_id':{'$nin':cids_wt}}):
        radio_pos_wt.add(c['rgz_name'])
    
    # Find radio positions that are only in weighted catalogs
    unique_to_wt = radio_pos_wt.difference(radio_pos_unwt)
    
    # Iterate over every entry in the unweighted catalog
    no_counterpart_unwt_wt, ir_mismatch = compare_ir(names=radio_pos_unwt, cat_1=catalog_0, cat_2=catalogs['catalog_{}'.format(weight)], \
                                                     cids_1=cids_unwt, cids_2=cids_wt)
    
    # Iterate over every entry that's only in the weighted catalogs
    no_counterpart_wt_unwt = no_counterparts(names=unique_to_wt, cat_1=catalogs['catalog_{}'.format(weight)], cids_1=cids_wt)
    
    return ir_mismatch, no_counterpart_unwt_wt, no_counterpart_wt_unwt

def get_ir(c):
    if 'ir_ra' in c['consensus']:
        return coord.SkyCoord(c['consensus']['ir_ra'], c['consensus']['ir_dec'], unit=(u.deg,u.deg))
    else:
        return coord.SkyCoord(0, 0, unit=(u.deg,u.deg))

def compare_ir(names, cat_1, cat_2, cids_1, cids_2):
    
    # Iterate over every entry in catalog 1 for matches in catalog 2
    # This returns the lists no_counterpart and ir_mismatch
    
    print '    Comparing {} and {} IR positions'.format(cat_1.name, cat_2.name)
    count = 0
    
    no_counterpart, ir_mismatch = [], []
    
    for rgz_name in names:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(names))
        
        # Get all entries with this position from each catalog
        c1s, c2s = [], []
        for c1 in cat_1.find({'rgz_name':rgz_name, 'catalog_id':{'$nin':cids_1}}):
            c1s.append(c1)
        for c2 in cat_2.find({'rgz_name':rgz_name, 'catalog_id':{'$nin':cids_2}}):
            c2s.append(c2)
        
        # If the weighted catalog has no match, log in 'no_counterpart'
        if not len(c2s):
            for c1 in c1s:
                no_counterpart.append((c1['catalog_id'], c1['zooniverse_id']))
            
        # Otherwise check each combination for an IR match
        else:
            for c1, c2 in itertools.product(c1s, c2s):
                ir_1 = get_ir(c1)
                ir_2 = get_ir(c2)
                
                # When a pair is found with same IR position, remove them from consideration
                if ir_1.separation(ir_2).arcsecond <= 3.:
                    try:
                        c1s.remove(c1)
                    except ValueError:
                        pass
                    try:
                        c2s.remove(c2)
                    except ValueError:
                        pass
            
            # For each combination that didn't match, log in 'ir_mismatch'
            for c1, c2 in itertools.product(c1s, c2s):
                ir_mismatch.append((c1['catalog_id'], c2['catalog_id']))
    
    return no_counterpart, ir_mismatch

def no_counterparts(names, cat_1, cids_1):
    
    # Iterate over every entry that's only in the weighted catalogs
    # These can just be logged in 'no_counterpart' without checking (by definition)
    
    print '    No counterpart between catalog and {}'.format(cat_1.name)
    count = 0
    
    no_counterpart = []
    
    for rgz_name in names:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(names))
        
        for c1 in cat_1.find({'rgz_name':rgz_name, 'catalog_id':{'$nin':cids_1}}):
            no_counterpart.append((c1['catalog_id'], c1['zooniverse_id']))
    
    return no_counterpart

def print_ir(filename, ir_mismatch):
    
    # Print the list of IR mismatches to a CSV
    
    with open(filename, 'w') as f:
        print >> f, 'unweighted_catalog_id,weighted_catalog_id'
        for c1, c2 in ir_mismatch:
            print >> f, '{},{}'.format(c1, c2)
    
    return None
                
######### End IR block; begin radio block #########

def get_radio_mismatch(no_counterpart_1, no_counterpart_2, cat_1, cat_2, cids_1, cids_2):
    
    print '    Find radio mismatch between {} and {}'.format(cat_1.name, cat_2.name)
    count = 0
    
    # Sort the sources by subject
    czids_1 = {}
    for cid, zid in no_counterpart_1:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(no_counterpart_1)+len(no_counterpart_2))
        
        if zid in czids_1:
            czids_1[zid].append(cid)
        else:
            czids_1[zid] = [cid]
    
    czids_2 = {}
    for cid, zid in no_counterpart_2:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(no_counterpart_1)+len(no_counterpart_2))
        
        if zid in czids_2:
            czids_2[zid].append(cid)
        else:
            czids_2[zid] = [cid]
    
    return czids_1, czids_2

def print_radio(filename, czid_unweighted, czid_weighted, weight):
    
    # Print the subject fields that contain radio sources without matches
    
    zids = set().union(czid_unweighted.keys(), czid_weighted.keys())
    
    with open(filename, 'w') as f:
        print >> f, 'zooniverse_id,catalog_id,catalog_weight'
        for zid in zids:
            if zid in czid_unweighted:
                for cid in czid_unweighted[zid]:
                    print >> f, '{},{},0'.format(zid, cid)
            if zid in czid_weighted:
                for cid in czid_weighted[zid]:
                    print >> f, '{},{},{}'.format(zid, cid, weight)
    
    return None

######### End radio block; begin post-processing block #########

def make_plots(n=None):

    cids = cids_for_removal()
    
    output_dir = '{}/csv/weighting_checks'.format(rgz_path)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plot_dir = '{}/plots'.format(output_dir)
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    if n is None or n == 1:
        consensus_level_prob(cids, plot_dir)

    if n is None or n == 2:
        multi_component_change(cids, plot_dir)

    if n is None or n == 3:
        angular_extent_vs_ir_change(cids, plot_dir)

    return None

def consensus_level_prob(cids, plot_dir):

    # What is the relationship between unweighted consensus level and the probability of change? (separately for IR and radio changes)

    print 'Comparing unweighted consensus level vs probability of consensus change (IR or radio)'

    for weight in [5, 10, 20]:

        print '    Processing catalog with weight = {}'.format(weight)

        # Load IR mismatches from file
        ir_dict = {}
        ir_mismatch_file = '/data/tabernacle/larry/RGZdata/rgz-analysis/csv/weighting_checks/ir_mismatch_{}.csv'.format(weight)
        with open(ir_mismatch_file, 'r') as f:
            raw_dict = csv.DictReader(f)
            count = 0
            for row in raw_dict:
                for i in row:
                    row[i] = int(row[i])
                ir_dict[count] = row
                count += 1

        # Get distribution of consensus level in general
        levels_orig = []
        for c in catalog_0.find({'catalog_id':{'$nin':cids['cids_0']}}):
            levels_orig.append(c['consensus']['radio_level'])

        # Get consensus levels for unweighted counterparts with changed IR
        levels_ir = []
        for ix, row in ir_dict.iteritems():
            c = catalog_0.find_one({'catalog_id':row['unweighted_catalog_id']})
            levels_ir.append(c['consensus']['radio_level'])

        # Group the sources into 5% bins for statistics (multiplied by 20 for ints)
        bins = np.array(range(21))

        binned_counts_orig = {}
        for b in bins:
            binned_counts_orig[b] = 0

        for level in levels_orig:
            binned_counts_orig[int(level*20)] += 1

        total = 0
        for b in binned_counts_orig:
            total += binned_counts_orig[b]

        assert total == len(levels_orig), '    For unweighted catalog: {} sources total, but {} were binned'.format(len(levels_orig), total)

        binned_counts_ir = {}
        for b in bins:
            binned_counts_ir[b] = 0

        for level in levels_ir:
            binned_counts_ir[int(level*20)] += 1

        total = 0
        for b in binned_counts_ir:
            total += binned_counts_ir[b]

        assert total == len(levels_ir), '    For catalog with weight {}: {} sources changed IR, but {} were binned'.format(weight, len(levels_ir), total)

        # Plot the distribution of consensus levels
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.hist(levels_orig, bins=bins/20., color='b', label='All sources')
        ax1.hist(levels_ir, bins=bins/20., color='g', label='Changed sources\nweight = {}'.format(weight))
        ax2.hist(levels_orig, bins=bins/20., color='b')
        ax2.hist(levels_ir, bins=bins/20., color='g')
        # zoom-in / limit the view to different portions of the data
        ax1.set_ylim(46000, 56000)  # outliers only
        ax2.set_ylim(0, 10000)  # most of the data
        # hide the spines between ax1 and ax2
        ax1.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax1.xaxis.tick_top()
        ax1.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        ax2.set_xlabel('Consensus level (radio)')
        fig.text(0.03, 0.5, '# sources', va='center', rotation='vertical')
        ax1.set_title('Consensus level of sources with changed IR')
        ax1.legend(loc='upper left')
        outfile = '{}/consensus_level_ir_{}.png'.format(plot_dir, weight)
        print '        Saving {}'.format(outfile)
        plt.savefig(outfile)

        # Find the fractional change as a function of consensus level
        frac_change = {}
        frac_change_err = {}
        for b in bins:
            if b in binned_counts_ir and binned_counts_ir[b]>0:
                frac_change[b] = 1.0*binned_counts_ir[b]/binned_counts_orig[b]
                frac_change_err[b] = frac_change[b] * np.sqrt( 1./binned_counts_ir[b] + 1./binned_counts_orig[b] )
            else:
                frac_change[b] = 0.
                frac_change_err[b] = 0.

        # Plot fractional change
        plt.clf()
        plt.gca().set_ylim([0, 0.45])
        plt.plot(bins/20., [frac_change[b] for b in bins], 's-', label='Fraction changed\nweight = {}'.format(weight))
        plt.errorbar(bins/20., [frac_change[b] for b in bins], yerr=[frac_change_err[b] for b in bins], linestyle='None', color='b')
        plt.xlabel('Consensus level (radio)')
        plt.ylabel('Fraction changed')
        plt.title('Fraction changed of sources with changed IR')
        plt.legend(loc='upper right')
        outfile = '{}/consensus_level_ir_frac_{}.png'.format(plot_dir, weight)
        print '        Saving {}'.format(outfile)
        plt.savefig(outfile)

        # Load radio mismatches from file
        radio_dict = {}
        radio_mismatch_file = '/data/tabernacle/larry/RGZdata/rgz-analysis/csv/weighting_checks/no_counterpart_{}.csv'.format(weight)
        with open(radio_mismatch_file, 'r') as f:
            raw_dict = csv.DictReader(f)
            count = 0
            for row in raw_dict:
                for i in ['catalog_id', 'catalog_weight']:
                    row[i] = int(row[i])
                radio_dict[count] = row
                count += 1

        # Get consensus levels for unweighted counterparts with changed radio
        levels_radio = []
        for ix, row in radio_dict.iteritems():
            if row['catalog_weight'] == 0:
                c = catalog_0.find_one({'catalog_id':row['catalog_id']})
                levels_radio.append(c['consensus']['radio_level'])

        # Group the sources into 5% bins for statistics (multiplied by 20 for ints)
        binned_counts_radio = {}
        for b in bins:
            binned_counts_radio[b] = 0

        for level in levels_radio:
            binned_counts_radio[int(level*20)] += 1

        total = 0
        for b in binned_counts_radio:
            total += binned_counts_radio[b]

        assert total == len(levels_radio), '    For catalog with weight {}: {} sources changed radio, but {} were binned'.format(weight, len(levels_radio), total)

        # Plot the distribution of consensus levels, splitting outlying 100% bin
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.hist(levels_orig, bins=bins/20., color='b', label='All sources')
        ax1.hist(levels_radio, bins=bins/20., color='g', label='Changed sources\nweight = {}'.format(weight))
        ax2.hist(levels_orig, bins=bins/20., color='b')
        ax2.hist(levels_radio, bins=bins/20., color='g')
        # zoom-in / limit the view to different portions of the data
        ax1.set_ylim(46000, 56000)  # outliers only
        ax2.set_ylim(0, 10000)  # most of the data
        # hide the spines between ax1 and ax2
        ax1.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax1.xaxis.tick_top()
        ax1.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        ax2.set_xlabel('Consensus level (radio)')
        fig.text(0.03, 0.5, '# sources', va='center', rotation='vertical')
        ax1.set_title('Consensus level of sources with changed radio')
        ax1.legend(loc='upper left')
        outfile = '{}/consensus_level_radio_{}.png'.format(plot_dir, weight)
        print '        Saving {}'.format(outfile)
        plt.savefig(outfile)

        # Find the fractional change as a function of consensus level
        frac_change = {}
        frac_change_err = {}
        for b in bins:
            if b in binned_counts_radio and binned_counts_radio[b]>0:
                frac_change[b] = 1.0*binned_counts_radio[b]/binned_counts_orig[b]
                frac_change_err[b] = frac_change[b] * np.sqrt( 1./binned_counts_radio[b] + 1./binned_counts_orig[b] )
            else:
                frac_change[b] = 0.
                frac_change_err[b] = 0.

        # Plot fractional change
        plt.clf()
        plt.gca().set_ylim([0, 0.8])
        plt.plot(bins/20., [frac_change[b] for b in bins], 's-', label='Fraction changed\nweight = {}'.format(weight))
        plt.errorbar(bins/20., [frac_change[b] for b in bins], yerr=[frac_change_err[b] for b in bins], linestyle='None', color='b')
        plt.xlabel('Consensus level (radio)')
        plt.ylabel('Fraction changed')
        plt.title('Fraction changed of sources with changed radio')
        plt.legend(loc='upper right')
        outfile = '{}/consensus_level_radio_frac_{}.png'.format(plot_dir, weight)
        print '        Saving {}'.format(outfile)
        plt.savefig(outfile)

    return None

def multi_component_change(cids, plot_dir):

    # Looking *only* at sources with more than 1 radio component, (how many are there?), how many of them change in radio?
    # For the ones with matched radio, how many of them change in IR?

    print 'Counting how many multi-component sources change in radio and IR'

    cids_multi = []
    for c in catalog_0.find({'radio.number_components':{'$gt':1}, 'catalog_id':{'$nin':cids['cids_0']}}):
        cids_multi.append(c['catalog_id'])

    counts_total = [len(cids_multi), len(cids_multi), len(cids_multi)]
    counts_ir, counts_radio = [], []
    weights = [5, 10, 20]

    for weight in weights:

        print '    Processing catalog with weight = {}'.format(weight)

        # Load IR mismatches from file
        ir_dict = {}
        ir_mismatch_file = '/data/tabernacle/larry/RGZdata/rgz-analysis/csv/weighting_checks/ir_mismatch_{}.csv'.format(weight)
        with open(ir_mismatch_file, 'r') as f:
            raw_dict = csv.DictReader(f)
            count = 0
            for row in raw_dict:
                for i in row:
                    row[i] = int(row[i])
                ir_dict[count] = row
                count += 1

        # Load radio mismatches from file
        radio_dict = {}
        radio_mismatch_file = '/data/tabernacle/larry/RGZdata/rgz-analysis/csv/weighting_checks/no_counterpart_{}.csv'.format(weight)
        with open(radio_mismatch_file, 'r') as f:
            raw_dict = csv.DictReader(f)
            count = 0
            for row in raw_dict:
                for i in ['catalog_id', 'catalog_weight']:
                    row[i] = int(row[i])
                radio_dict[count] = row
                count += 1

        # Count the number of multi-component sources in each category
        count_ir = 0
        for ix, row in ir_dict.iteritems():
            if row['unweighted_catalog_id'] in cids_multi:
                count_ir += 1
        
        counts_ir.append(count_ir)

        count_radio = 0
        for ix, row in radio_dict.iteritems():
            if row['catalog_weight'] == 0 and row['catalog_id'] in cids_multi:
                count_radio += 1
        
        counts_radio.append(count_radio)

    # Plot number of sources
    plt.clf()
    plt.plot(weights, counts_total, 'bs-', label='All sources (unweighted)')
    plt.plot(weights, counts_ir, 'go-', label='IR change')
    plt.plot(weights, counts_radio, 'r^-', label='Radio change')
    plt.xticks(np.arange(0, 21, 5))
    plt.xlabel('Weight')
    plt.ylabel('# sources')
    plt.title('Number of changed multi-component sources')
    plt.legend(loc='best')
    outfile = '{}/multi_comp_counts.png'.format(plot_dir)
    print '        Saving {}'.format(outfile)
    plt.savefig(outfile)

    # Plot fractional change
    frac_ir = 1.*np.array(counts_ir)/counts_total[0]
    frac_err_ir = [f * np.sqrt( 1./c + 1./counts_total[0] ) for f, c in zip(frac_ir, counts_ir)]
    frac_radio = 1.*np.array(counts_radio)/counts_total[0]
    frac_err_radio = [f * np.sqrt( 1./c + 1./counts_total[0] ) for f, c in zip(frac_radio, counts_radio)]
    
    plt.clf()
    plt.errorbar(weights, frac_ir, yerr=frac_err_ir, fmt='go-', label='IR change')
    plt.errorbar(weights, frac_radio, yerr=frac_err_radio, fmt='r^-', label='Radio change')
    plt.xticks(np.arange(0, 21, 5))
    plt.xlabel('Weight')
    plt.ylabel('Fraction changed')
    plt.title('Fraction changed of multi-component sources')
    plt.legend(loc='best')
    outfile = '{}/multi_comp_frac.png'.format(plot_dir)
    print '        Saving {}'.format(outfile)
    plt.savefig(outfile)
    
    return None

def angular_extent_vs_ir_change(cids, plot_dir):

    # Looking *only* at sources with 1 radio component, what is the relationship between source angular extent and the probability of an IR change?

    print 'Comparing angular extent vs probabilty of IR change'

    for weight in [5, 10, 20]:

        print '    Processing catalog with weight = {}'.format(weight)

        # Load IR mismatches from file
        ir_dict = {}
        ir_mismatch_file = '/data/tabernacle/larry/RGZdata/rgz-analysis/csv/weighting_checks/ir_mismatch_{}.csv'.format(weight)
        with open(ir_mismatch_file, 'r') as f:
            raw_dict = csv.DictReader(f)
            count = 0
            for row in raw_dict:
                for i in row:
                    row[i] = int(row[i])
                ir_dict[count] = row
                count += 1

        # Get angular extents for all 1-component sources
        extents_changed = []
        for ix, row in ir_dict.iteritems():
            c = catalog_0.find_one({'catalog_id':row['unweighted_catalog_id']})
            if c['radio']['number_components'] == 1:
                extents_changed.append(c['radio']['max_angular_extent'])

        extents_orig = []
        for c in catalog_0.find({'radio.number_components':1, 'catalog_id':{'$nin':cids['cids_0']}}):
            extents_orig.append(c['radio']['max_angular_extent'])

        # Group the sources into 10" bins for statistics
        bins = 10*np.array(range(int(max(extents_orig)/10.)+2))

        binned_counts_changed = {}
        for b in bins:
            binned_counts_changed[b] = 0

        for extent in extents_changed:
            binned_counts_changed[10*int(extent/10.)] += 1

        total = 0
        for b in binned_counts_changed:
            total += binned_counts_changed[b]

        assert total == len(extents_changed), '    For catalog with weight {}: {} sources changed, but {} were binned'.format(weight, len(extents_changed), total)

        binned_counts_orig = {}
        for b in bins:
            binned_counts_orig[b] = 0

        for extent in extents_orig:
            binned_counts_orig[10*int(extent/10.)] += 1

        total = 0
        for b in binned_counts_orig:
            total += binned_counts_orig[b]

        assert total == len(extents_orig), '    For catalog with weight {}: {} sources originally, but {} were binned'.format(weight, len(extents_orig), total)

        # Plot the distribution of angular extents
        plt.clf()
        plt.hist(extents_orig, bins=bins, label='All sources')
        plt.hist(extents_changed, bins=bins, label='Changed sources\nweight = {}'.format(weight))
        plt.xlabel('Angular extent (")')
        plt.ylabel('# sources')
        plt.title('Angular extent of 1-component sources')
        plt.legend()
        outfile = '{}/angular_extent_{}.png'.format(plot_dir, weight)
        print '        Saving {}'.format(outfile)
        plt.savefig(outfile)

        # Find the fractional change as a function of angular extent
        frac_change = {}
        frac_change_err = {}
        for b in bins:
            if b in binned_counts_changed and binned_counts_changed[b]>0:
                frac_change[b] = 1.0*binned_counts_changed[b]/binned_counts_orig[b]
                frac_change_err[b] = frac_change[b] * np.sqrt( 1./binned_counts_changed[b] + 1./binned_counts_orig[b] )
            else:
                frac_change[b] = 0.
                frac_change_err[b] = 0.

        # Plot fractional change
        plt.clf()
        plt.gca().set_ylim([0, 1])
        plt.plot(bins, [frac_change[b] for b in bins], 's-', label='Fraction changed\nweight = {}'.format(weight))
        plt.errorbar(bins, [frac_change[b] for b in bins], yerr=[frac_change_err[b] for b in bins], linestyle='None', color='b')
        plt.xlabel('Angular extent (")')
        plt.ylabel('Fraction changed')
        plt.title('Fraction changed of 1-component sources')
        plt.legend(loc='upper left')
        outfile = '{}/angular_extent_frac_{}.png'.format(plot_dir, weight)
        print '        Saving {}'.format(outfile)
        plt.savefig(outfile)

    return None

######### End post-processing block #########

if __name__ == "__main__":

    output_dir = '{}/csv/weighting_checks'.format(rgz_path)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    print 'Finding duplicate entries'
    cids = cids_for_removal()
    print '    Removed {}, {}, {}, and {} duplicate entries'.format(len(cids['cids_0']), len(cids['cids_5']), len(cids['cids_10']), len(cids['cids_20']))

    for weight in [5, 10, 20]:

        print 'Testing catalog with weight {}'.format(weight)

        # IR
        print 'Finding IR mismatches'
        ir_mismatch, no_counterpart_unwt_wt, no_counterpart_wt_unwt = get_ir_mismatch(cids_unwt=cids['cids_0'], cids_wt=cids['cids_{}'.format(weight)], weight=weight)
        ir_filename = '{}/ir_mismatch_{}.csv'.format(output_dir, weight)
        print 'IR checks complete; printing {} IR mismatches from catalog_{} to {}'.format(len(ir_mismatch), weight, ir_filename)
        print_ir(ir_filename, ir_mismatch)

        # Radio
        print 'Finding radio mismatches'
        czid_weighted, czid_unweighted = get_radio_mismatch(no_counterpart_1=no_counterpart_wt_unwt, no_counterpart_2=no_counterpart_unwt_wt, \
                                                            cat_1=catalogs['catalog_{}'.format(weight)], cat_2=catalogs['catalog_0'], \
                                                            cids_1=cids['cids_{}'.format(weight)], cids_2=cids['cids_0'])
        radio_filename = '{}/no_counterpart_{}.csv'.format(output_dir, weight)
        print 'Radio checks complete; printing {} radio mismatches from catalog_{} to {}'.format(len(czid_weighted)+len(czid_unweighted), weight, radio_filename)
        print_radio(radio_filename, czid_unweighted, czid_weighted, weight)
