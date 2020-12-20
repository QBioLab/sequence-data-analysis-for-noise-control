# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:47:47 2019

@author: chena
"""

from pyfaidx import Fasta
import pandas as pd

insert = Fasta('../rawdata/reference/inserts1.fa')
hg19 = Fasta('../../../Reference/hg19.fa')
insert_sites = pd.read_csv('../rawdata/reference/insert_sites.csv')
insert_sites = insert_sites.sort_values(by=['chromosome','start'])
output = '../metadata/reference/hg19_insert.fa'
tss_output= '../metadata/reference/tss_sites.txt'
new_site_output = '../metadata/reference/new_sites.txt'
tss_site = []
new_site = []
with open(output,'w') as fout:
    for i in hg19:
        name = i.name
        seq = str(i).upper()
        potential_insert = insert_sites[insert_sites['chromosome'] == name]
        potential_insert_len = len(potential_insert)
        if potential_insert_len > 0:
            bias = 0
            for j in range(potential_insert_len):
                d = potential_insert.iloc[j,1]
                insert_name = potential_insert.iloc[j,0]
                if d == 'forward':
                    insert_seq = insert[potential_insert.iloc[j,0]][:].seq
                    if insert_name == 'CEG':
                        tss = 2472
                    elif insert_name == 'mRN':
                        tss = 767
                elif d == 'reverse':
                    insert_seq = insert[potential_insert.iloc[j,0]][:].complement.reverse.seq
                    if insert_name == 'CEG':
                        tss = 4928
                    elif insert_name == 'mRN':
                        tss = 2597
                start = potential_insert.iloc[j,3] + bias - 1
                end = potential_insert.iloc[j,4] + bias - 1
                tss_site.append('{}-{}-{}-{}'.format(insert_name, d, name, 
                                start+tss))
                new_site.append('{},{},{},{},{}'.format(insert_name, d, name, 
                                start+1,end+1+len(insert_seq)))
                seq = '{}{}{}'.format(seq[:start+1],insert_seq,seq[end:])
                bias += len(insert_seq) - 4
        fout.write('>{}\n'.format(name))
        fout.write('{}\n'.format(seq))
with open(tss_output,'w') as fout:
    for s in tss_site:
        fout.write('{}\n'.format(s))
with open(new_site_output,'w') as fout:
    for s in new_site:
        fout.write('{}\n'.format(s))
        
        