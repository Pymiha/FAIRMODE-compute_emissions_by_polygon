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



shape_id_field = 'NUTS_ID'

shapes_path = pathlib.Path('N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/FAIRMODE-compute_emissions_by_polygon-master/data/polygons/NUTS_RG_01M_2016_4326')
shape_fname = 'NUTS_RG_01M_2016_4326.shp' 

# path to raster files
grid_path = pathlib.Path(
    'N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/split/small')

results_path = pathlib.Path(
    'N:/Inby/Aktive-prosjekter/b120011_EMISYS/Emis_2_FAIRMODE/results')

year = 2019

unit_factor = 10**-3 # if input data in tons - 10**-3; output is always in kton.

project = 'EMISYS'
tag = 'NUTS2016'
########################

shapes, attributes, meta = ger.read_shapefile_fiona_disjoint_polygons(shapes_path, shape_fname, fix_geometry=True,
                                                                      crs_in=4326, crs_out=3035, transf=True)

shapes, attributes = ger.filter_shapes(shapes, attributes, 
                                   filter_fields=['CNTR_CODE', 'LEVL_CODE'], 
                                   filter_values=[['NO'], [3]], # NO is CNTR_CODE for Norway, 3 is LEVL_CODE of NUTS 
                                   min_area=2*10**6) # min area - smallest shape to use

# plot shapefile to make sure correct polygons are used.
ger.write_shapefile_fiona(shapes, attributes, meta, grid_path, 'Norway_NUTS.shp')

shape_id_dictionary = ger.make_shape_id_dictionary(attributes, shape_id_field='NUTS_ID', 
                             other_fields=['CNTR_CODE', 'NAME_LATN', 'NUTS_NAME'])

c = 0
sizes = []
for shape in shapes:
    for sh in shape:
        c += 1
        sizes.append(sh.area)
print('Number of shapes', c)

########
raster_fnames = ['EMISYS_NILU_NO_NOX_GNFRF_epsg3035_2019_Traffic.tif']

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
        with open(grid_path / 'cell_connections_NUTS.dct', 'rb') as fil:
            cell_conn = pickle.load(fil)
            print('Number of already connected grid cells', len(list(cell_conn.keys())))
    except:
        cell_conn = {}
    
    # NUTS does not need never_connect, since NUTS covers whole Norway. If there were many offshore cells it 
    # could make sense to use it for NUTS as well. In case of FUA it always make sense to use it.
    result, cell_conn, _ = ger.intersections(shapes, attributes, r_data, cell_conn,
                                      field=shape_id_field, raster_fname=raster_fname, unit_factor=unit_factor,
                                      allowed_unallocated=0.05)[:]
    
    for shp_id in shape_id_dictionary:
        if shp_id not in result:
            print(shp_id, 'not in result. adding it with a 0.')
            result[shp_id] = 0
    #print(result)
    
    
    with open(grid_path / 'cell_connections_NUTS.dct', 'wb') as fil:
        pickle.dump(cell_conn, fil)
        
    ger.write_result_file(results_path, res_fname, result, shape_id_dictionary, year, pollutant, sector,
                          fields=['NUTS_ID', 'CNTR_CODE', 'NAME_LATN', 'NUTS_NAME', 'POLLUTANT',
                                    'YEAR', 'GNF_SECTOR', 'EMIS(kTons)'])

