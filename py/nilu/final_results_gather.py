# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:23:56 2023

@author: mmar
"""

import pathlib
import os
import sys


def write_result_file(directory, fname, header, data):
    with open(directory / fname, 'w') as fil:
        fil.write(','.join(header))
        fil.write('\n')
        for dat in data:
            fil.write(','.join(dat))
            fil.write('\n')
            
            
            
def read_result_file(directory, fname):
    fields = ['NUTS_ID', 'CNTR_CODE', 'NAME_LATN', 'NUTS_NAME', 'POLLUTANT',
              'YEAR', 'GNF_SECTOR', 'EMIS(kTons)']
    header = []
    data = []
    with open(directory / fname, 'r') as fil:
        for (i, line) in enumerate(fil):
            line = line.strip().split(',')
            if i == 0:
                header = line[:]
            else:
                data.append(line)
    
    return header, data
    
results_path = pathlib.Path(
    'N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/results/') # directory with single pollutant/sector combination!!!
final_path = pathlib.Path(
    'N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/results/final') # must be different directory than results for single pollutant/sector combination!!!


project = 'EMISYS'
org = 'NILU'

tag = 'NUTS2016'
result_fname = f'{project}_{org}_{tag}.csv'
########################

header = []
datas = []
for fname in os.listdir(results_path):
    if fname.startswith(project) and tag in fname:
        print(fname)
        h, d = read_result_file(results_path, fname)
        
        if not header:
            header = h[:]
        
        if h != header:
            print('Headers are not the same!!', fname, header, h)
        
        datas += d
    

write_result_file(final_path, result_fname, header, datas)
        


tag = 'FUA2020'
result_fname = f'{project}_{org}_{tag}.csv'
########################

header = []
datas = []
for fname in os.listdir(results_path):
    if fname.startswith(project) and tag in fname:
        print(fname)
        h, d = read_result_file(results_path, fname)
        
        if not header:
            header = h[:]
        
        if h != header:
            print('Headers are not the same!!', fname, header, h)
        
        datas += d
    

write_result_file(final_path, result_fname, header, datas)
        
        



