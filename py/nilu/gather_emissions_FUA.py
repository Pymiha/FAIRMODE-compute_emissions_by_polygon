# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:02:46 2022

@author: mmar
"""


from shapely.geometry import LineString, mapping, Polygon, Point, MultiPolygon, box
from shapely.validation import make_valid
import fiona
import pathlib
import os
import sys
from pyproj import CRS
from pyproj import Transformer
import numpy.polynomial.polynomial as pol
import rasterio
import datetime
import numpy as np
import pickle

import gather_emissions_routines as ger




shape_id_field = 'URAU_CODE'

shapes_path = pathlib.Path('N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/FAIRMODE-compute_emissions_by_polygon-master/data/polygons/URAU_RG_100K_2020_4326_FUA')
shape_fname = 'URAU_RG_100K_2020_4326_FUA.shp' 

# path to raster files
grid_path = pathlib.Path(
    'N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/split/small')

results_path = pathlib.Path(
    'N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/results')

year = 2019

unit_factor = 10**-3 # if input data in tons - 10**-3; output is always in kton.

# project and tag get added to filename of single sector/pollutant result file. 
# in code final_results_gather, they are used to collect correct files.
project = 'EMISYS'
tag = 'FUA2020'
########################

shapes, attributes, meta = ger.read_shapefile_fiona_disjoint_polygons(shapes_path, shape_fname, fix_geometry=True,
                                                                      crs_in=4326, crs_out=3035, transf=True)

shapes, attributes = ger.filter_shapes(shapes, attributes, 
                                   filter_fields=['CNTR_CODE'], 
                                   filter_values=[['NO']], # NO is CNTR_CODE for Norway,
                                   min_area=2*10**6) # min area - smallest shape to use

# plot shapefile to make sure correct polygons are used.
ger.write_shapefile_fiona(shapes, attributes, meta, grid_path, 'Norway_FUA.shp')

shape_id_dictionary = ger.make_shape_id_dictionary(attributes, shape_id_field='URAU_CODE', 
                             other_fields=['CNTR_CODE', 'URAU_NAME', 'NUTS3_2016'])

c = 0
sizes = []
for shape in shapes:
    for sh in shape:
        c += 1
        sizes.append(sh.area)
print('Number of shapes', c)

########
raster_fnames = ['EMISYS_NILU_NO_NOX_GNFRF_epsg3035_2019_Traffic.tif'] # example of raster file name

raster_fnames = []
for fname in os.listdir(grid_path):
    if fname.startswith(project):
        raster_fnames.append(fname)
        
#raster_fnames = ['EMISYS_NILU_NO_NMVOC_GNFRF_epsg3035_2019_Traffic.tif']
########

for raster_fname in raster_fnames:
    print('\n', raster_fname)
    org = raster_fname.split('_')[1]
    pollutant = raster_fname.split('_')[3]
    sector = raster_fname.split('_')[4]
    r_data = ger.read_raster(grid_path, raster_fname)
    
    res_fname = f'{project}_{org}_{sector}_{pollutant}_{tag}.csv'
    try:
        with open(grid_path / 'cell_connections_FUA.dct', 'rb') as fil:
            cell_conn = pickle.load(fil)
            print('Number of already connected grid cells', len(list(cell_conn.keys())))
    except:
        cell_conn = {}
        
    try:
        with open(grid_path / 'never_connect_FUA.set', 'rb') as fil:
            never_connect = pickle.load(fil)
            print('Number of grid cells already determined to not be overlapping with any polygon', len(never_connect))
    except:
        never_connect = set([])
    
    # FUA covers only small part of Norway, so allowed_unallocated=0.9. For NUTS it is allowed_unallocated=0.05 (also default value in function)
    result, cell_conn, never_connect = ger.intersections(shapes, attributes, r_data, cell_conn,
                                      field=shape_id_field, raster_fname=raster_fname, unit_factor=unit_factor,
                                      allowed_unallocated=0.9, never_connect=never_connect)[:] 
    
    for shp_id in shape_id_dictionary:
        if shp_id not in result:
            print(shp_id, 'not in result. adding it with a 0.')
            result[shp_id] = 0
    #print(result)
    
    with open(grid_path / 'cell_connections_FUA.dct', 'wb') as fil:
        pickle.dump(cell_conn, fil)
    
    with open(grid_path / 'never_connect_FUA.set', 'wb') as fil:
        pickle.dump(never_connect, fil)
        
    ger.write_result_file(results_path, res_fname, result, shape_id_dictionary, year, pollutant, sector,
                          fields=['URAU_CODE', 'CNTR_CODE', 'URAU_NAME', 'NUTS3_2016', 'POLLUTANT',
                                    'YEAR', 'GNF_SECTOR', 'EMIS(kTons)'])

