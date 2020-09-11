import pyBigWig
from collections import defaultdict
from matplotlib import pyplot as plt
import numpy as np

site = defaultdict(list)
with open('../metadata/reference/tss_sites.txt') as fin:
    for line in fin:
        line = line.strip().split('-')
        site[line[0]].append([line[1],line[2],int(line[3])])

#ChIP-seq H3K27ac        
bws = {}
# 1-H3K27AC_combined: ChIP-seq_for_HeLa-AB1_H3K27ac_AM_25uw_low
# 2-H3K27AC_combined: ChIP-seq_for_HeLa-AB1_H3K27ac_AM_25uw_high
# 3-H3K27AC_combined: ChIP-seq_for_HeLa-AB1_H3K27ac_Dark_control
prefix = ['1-H3K27AC_combined','2-H3K27AC_combined','3-H3K27AC_combined']
for p in prefix:
    bws[p] = pyBigWig.open('../metadata/chip-seq/{}.SeqDepthNorm.bw'.format(p))

fig, ax = plt.subplots()
length = 1000
data = {}
for p in prefix:
    depth = bws[p]
    total = np.array(depth.values(site['mRN'][0][1],
                                  site['mRN'][0][2]-length,
                                  site['mRN'][0][2]+length))
    for s in site['mRN'][1:]:
        values = depth.values(s[1],s[2]-length,s[2]+length)
        if s[0] == 'forward':
            total += np.array(values)
        else:
            values.reverse()
            total += np.array(values)
    data[p] = total
    plt.plot(np.array(range(0,2*length))-length,total)
plt.legend(prefix)

#ATAC-seq
bws = {}
# high: ATAC-seq_for_HeLa-AB1_AM_25uw_high
# low: ATAC-seq_for_HeLa-AB1_AM_25uw_low
# A1: ATAC-seq_for_HeLa-AB1_Dark_control
prefix = ['high','low','A1']
for p in prefix:
    bws[p] = pyBigWig.open('../metadata/atac-seq/{}.SeqDepthNorm.bw'.format(p))
    
fig, ax = plt.subplots()
length = 1000
for p in prefix:
    depth = bws[p]
    total = np.array(depth.values(site['mRN'][0][1],
                                  site['mRN'][0][2]-length,
                                  site['mRN'][0][2]+length))
    for s in site['mRN'][1:]:
        values = depth.values(s[1],s[2]-length,s[2]+length)
        if s[0] == 'forward':
            total += np.array(values)
        else:
            values.reverse()
            total += np.array(values)
    data[p] = total
    plt.plot(np.array(range(0,2*length))-length,total)
plt.legend(prefix)

#scATAC-seq

prefix = []
with open('../rawdata/scatac-seq/HGC20191230001-0003_lane7/L7_md5sum.check.out') as fin:
    pattern = '_R1_001.fastq'
    for line in fin:
        if pattern in line and 'HW' in line:
            line = line.split(' ')[0]
            line = line.replace(pattern,' ')
            line = line.split(' ')[0]
            prefix.append(line)
with open('../rawdata/scatac-seq/HGC20191230001-0003_lane8/L8_md5sum.check.out') as fin:
    pattern = '_R1_001.fastq'
    for line in fin:
        if pattern in line and 'HW' in line:
            line = line.split(' ')[0]
            line = line.replace(pattern,' ')
            line = line.split(' ')[0]
            prefix.append(line)

bws = {}
for p in prefix:
    bws[p] = pyBigWig.open('../metadata/scatac-seq/{}.SeqDepthNorm.bw'.format(p))

# HW1: scATAC-seq_for_HeLa-AB1_150min_light_on_80uw_384
# HW2: scATAC-seq_for_HeLa-AB1_dark_control_384
# HW3: scATAC-seq_for_HeLa-AB1_600min_light_on_80uw_384
# HW4: scATAC-seq_for_HeLa-AB1_750min_light_on_80uw_384
# HW5: scATAC-seq_for_HeLa-AB1_1200min_light_on_80uw_384
# HW6: scATAC-seq_for_HeLa-AB1_1350min_light_on_80uw_384
exp = ['HW1','HW2','HW3','HW4','HW5','HW6']
length = 1000
data = {}
for e in exp:
    total_exp = []
    for p in prefix:
        if e in p:
            depth = bws[p]
            try:
                total = np.array(depth.values(site['mRN'][0][1],
                                          site['mRN'][0][2]-length,
                                          site['mRN'][0][2]+length))
            except RuntimeError:
                print('mRN',e,p,s)
            for s in site['mRN'][1:]:
                try:
                    values = depth.values(s[1],s[2]-length,s[2]+length)
                    if s[0] == 'forward':
                        total += np.array(values)
                    else:
                        values.reverse()
                        total += np.array(values)
                except RuntimeError:
                    print('mRN',e,p,s)
            total_exp.append(total)
    total_exp = np.array(total_exp)
    data[e] = total_exp