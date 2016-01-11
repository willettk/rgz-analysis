import consensus

zids=("hymor01_ARG0003o4o", "hymor01_ARG0003o4w", "hymor01_ARG0003o4t", "hymor02_ARG0003hkm", "hymor02_ARG0003hlo", "hymor02_ARG0003hk3", "hymor03_ARG0000gng", "hymor04_ARG0003czn", "hymor04_ARG0003cyx", "hymor04_ARG0003cyr", "hymor04_ARG0003cyc", "hymor05_ARG0001arg", "hymor05_ARG0001ark", "hymor05_ARG0001arv", "hymor05_ARG0001aqn", "hymor06_ARG00021rf", "hymor07_ARG00024bn", "hymor07_ARG00024d2", "hymor07_ARG0002495", "hymor07_ARG000248y", "hymor08_ARG0003j0n", "hymor08_ARG0003j0a", "hymor08_ARG0003j0z", "hymor08_ARG0003j1d", "hymor09_ARG00027o7", "hymor09_ARG00027p0", "hymor09_ARG00027mr", "hymor10_ARG0003ph9", "hymor10_ARG0003pgv", "hymor10_ARG0003pgw", "hymor10_ARG0003pht", "hymor11_ARG0002p17", "hymor11_ARG0002p1h", "hymor11_ARG0002p1x", "hymor11_ARG0002ozg", "hymor11_ARG0002oza" "hymor12_ARG00007vb", "hymor12_ARG00007us", "hymor12_ARG00007uw", "hymor12_ARG00007up", "hymor12_ARG00007vp", "hymor13_ARG0001hz4", "hymor13_ARG0001hyx", "hymor13_ARG0001hzl", "hymor13_ARG0001hys")

for zid in zids:
    c = consensus.checksum(zid.split('_')[1])
    if c != None:
        consensus.plot_consensus(c,save_fig=zid)
    else:
        print "Couldn't do %s" % zid
    
