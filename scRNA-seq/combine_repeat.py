# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 21:49:01 2020

@author: chena
"""
import pandas as pd

def combine(prefix):
    table1 = pd.read_csv("../result/scRNA_seq_2020_11_1-7/{}_1.csv".format(prefix))
    table2 = pd.read_csv("../result/scRNA_seq_2020_11_1-7/{}_2.csv".format(prefix))
    barcode = table1.columns[0]
    table = pd.concat([table1,table2],axis=0)
    del table1
    del table2
    result = table.groupby(barcode).sum()
    del table
    result.to_csv('../result/scRNA_seq_2020_11_1-7/{}_withInde.csv'.format(prefix))
    
prefixs = ['A1A2','B1A3','C1C2']
for p in prefixs:
    combine(p)
    
    