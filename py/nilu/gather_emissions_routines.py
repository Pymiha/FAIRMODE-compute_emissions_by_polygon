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



def write_shapefile_fiona(shapes, attributes, meta, path_out: pathlib.Path,
                    name_out: str):
    """
    Write shapefile.

    Parameters
    ----------
    shapes : list of Shapes
        shapes to write

    attributes: list of attributes - one dictionary per shape
        attributes of shapes
    path_out : pathlib.Path
        path to resulting file
    meta : dict
        produced when Fiona reads shapefile
    name_out : str
        file name - .shp included 

    Returns
    -------
    None.
    Writes shapefile to provided path. 
    """
    if len(shapes) != len(attributes):
        raise Exception(
            f'Lenght of shapes is not the same as of attributes! {len(shapes)}, {len(attributes)}')
    meta.update({'encoding': 'utf-8'})
    succeeded = True
    with fiona.open(path_out / name_out, 'w', **meta) as dst:
        for (shp, attr) in zip(shapes, attributes):
            #print(shp)
            try:
                for elem in shp:
                    dst.write({'geometry': mapping(elem),
                              'properties': attr
                               }
                              )
            except Exception as e:
                #print(e)
                try:
                    for item in attr:
                        
                        try:
                            attr[item] = attr[item].decode('utf-8')
                        except:
                            pass
                        
                    dst.write({'geometry': mapping(shp),
                              'properties': attr
                               }
                              )
                    
                except Exception as e:
                    #print(2, e)
                    succeeded = False
                    print('attr where FAILED:', attr, e)
    
    if not succeeded:
        print(f'Saving file {name_out} FAILED. Deleting. Have to fix the issue! Most likely encoding of values')
        for fname in os.listdir(path_out):
            if name_out.replace('.shp', '') in fname:
                os.remove(path_out / fname)
        print('Exiting')
        sys.exit()

    '''
    with open(path_out / name_out.replace('.shp', '_field_names.txt'), 'w') as fil:
        for attribute in meta['schema']['properties']:
            fil.write(f'{attribute} {attribute[:10]}\n')
    '''
    
    
def prep_coords(coords, proj, transf: bool = True, multi: bool = True):
    '''
    Pay attention to these two lines:
        extr = [proj.transform(p[1], p[0])[::-1] for p in extr]
        holes = [proj.transform(p[1], p[0])[::-1] for p in holes]
    after transforming coordinate get switched - maybe 3035 has them in y, x order,
    maybe other crs as well. Writting out filtered shapefile and seeing it on map 
    is best way of knowing if it is done correctly. 
    
    
    MultiPolygon([(((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)),
                   [((0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1))])]
                 )

    MultiPolygon([(((0, 0), (3, 3), (0, 3)), 
                   [((0.5, 0.5), (1, 1), (0.5, 1)), 
                    ((1.5, 1.5), (2, 2), (1.5, 2))])])
    
    MultiPolygon([((coordinates of polygon), 
                   [(coordinates of 1st hole), 
                    (coordinates of 2nd hole)])])
    
    MultiPolygon([((coordinates of polygon), 
                   [(coordinates of 1st hole), 
                    (coordinates of 2nd hole)]),
                  ((2. coordinates of polygon), 
                  [(2. coordinates of 1st hole), 
                  (2. coordinates of 2nd hole)])])
    
    '''
    ret = []
    extr = coords[0]
    holes = coords[1:]
    #print(extr, holes)
    
    # forget about holes!!
    holes = []
    
    #print(holes)
    
    if transf:
        '''
        for p in extr:
            print(p, proj.transform(p[1], p[0]))
        '''
        # [::-1] - 3035 needs this. other projectioons maybe not.
        extr = [proj.transform(p[1], p[0])[::-1] for p in extr]
        holes = [proj.transform(p[1], p[0])[::-1] for p in holes]
        
    if multi:
        ret = [((extr), [tuple(h) for h in holes])]
    
    else:
        ret = extr
    
    return ret

def read_shapefile_fiona_disjoint_polygons(file_path: pathlib.Path = None, file_name: str = None, fix_geometry=True,
                                           crs_in: int = 4326, crs_out: int = 3035, transf: bool = True):
    '''
    Read shapefile using fiona and correctly handle multipolygons. Multipolygons are made from 2 or more unconnected polygons.
    Have to be done correctly, so that all individual polygons in shapefile are used in all steps of process - for instance 
    when connecting shapes with grid cells.

    Parameters
    ----------
    file_path : pathlib.Path, optional
        path to shapefile
    file_name : str, optional
        file name of shapefile
    fix_geometry : bool, optional
        Flag to fix geometry - advisable to keep it True. The default is True.
    crs_in : int, optional
        crs of shapefile. The default is 4326.
    crs_out : int, optional
        Which crs to reproject shapefile to. The default is 3035.
    transf : bool, optional
        Flag to actually reproject. If it shapefile does not need to be reprojected,
        can keep crs in arguments, but put this to False.. The default is True.

    Returns
    -------
    ret_polys : TYPE
        polygons and multipolygons
    ret_attributes : TYPE
        attributes - one entry per polygon or multipolygon.
    ret_meta : TYPE
        Meta of shapefile, already has crs data coresponding to crs_out (if reproject).

    '''
    # https://gis.stackexchange.com/questions/70591/creating-shapely-multipolygons-from-shapefile-multipolygons
    # https://shapely.readthedocs.io/en/stable/manual.html#MultiPolygon
    ret_polys = []
    ret_attributes = []
    ret_meta = None
    used_at_zero = 0
    
    proj = None
    if transf:
        crs1 = CRS.from_epsg(crs_in)
        crs2 = CRS.from_epsg(crs_out)
        proj = Transformer.from_crs(crs1, crs2)
    
    print(f'Reading file {file_name} from drive')
    with fiona.open(file_path / file_name) as src:
        ret_meta = src.meta
        for line in src:
            preped = None
        
            pom_polys = []
            polis = line['geometry']['coordinates']
            
            if len(polis) == 1:
                preped = prep_coords([polis[0]], proj, transf, multi=False)
                poly = Polygon(preped)
                pom_polys.append(poly)
            
            else:
                if isinstance(polis[0][0], tuple):
          
                    preped = prep_coords(polis, proj, transf)
                    
                    mpoly = MultiPolygon(preped)
                    pom_polys.append(mpoly)
                    
                else:
                    for elem in polis:
                        preped = prep_coords(elem, proj, transf)
                        
                        '''
                        if len(preped[0][1]) < 4:
                            print('\n')
                        print('pr 01', preped[0][1], type(preped[0][1]))
                        '''
                        
                        mpoly = MultiPolygon(preped)
                        
                        pom_polys.append(mpoly)
                        
                       
            
            if fix_geometry:                
                for (pi, pl) in enumerate(pom_polys):
                    if not pl.is_valid:
                        #print('Invalid polygon')
                        pom_polys[pi] = make_valid(pl)
                        #print(pom_polys[pi].is_valid)
                        
            
            
            #print(len(pom_polys))
            ret_polys.append(pom_polys)
            ret_attributes.append(line['properties'])
          
    if transf:
        crs = CRS(crs_out)
        ret_meta['crs_wkt'] = crs.to_wkt()
            
    return ret_polys, ret_attributes, ret_meta


def filter_shapes(shapes, attributes, filter_fields=['CNTR_CODE', 'LEVL_CODE'], filter_values=[['NO'], [3]], min_area=10**6):
    '''
    Take shapefile and attributes and keep only relevant ones.

    Parameters
    ----------
    shapes : list
        list of polygons and multipolygons, coming from read_shapefile_fiona_disjoint_polygons(...)
    attributes : list
        list of attributes, one per polygon and multipolygon
    filter_fields : list, optional
        List of fields in attributes of shapefile to use for filtering. The default is ['CNTR_CODE', 'LEVL_CODE'].
    filter_values : list[list], optional
        List of lists. Each filter_field in filter_fields must have corresponding list of allowed 
        values in filter_values. For instance if Sweden, Norway and Finland were to be used in NUTS, filter_values would be:
            filter_values = [['SE', 'NO', 'FI'], [3]]. The default is [['NO'], [3]].
    min_area : int, optional
        Smallest area of polygon for it to still be used. In Norway it removes many tiny islands. The default is 10**6.

    Returns
    -------
    ret_sh : list
        selected polygons and multipolygons
    ret_attr : list
        selected attributes.
    
    len(ret_sh) = len(ret_attr)
    '''
    
    ret_sh = []
    ret_attr = []
    
    assert len(filter_fields) == len(filter_values), 'function filter_shapes(): Filter fields and filter values have different lenghts.'
    
    for sh, attr in zip(shapes, attributes):
        sw = True
        # check if fields have desired attributes
        for (ff, fv) in zip(filter_fields, filter_values):
            if attr[ff] not in fv:
                sw = False
        if sw:
            # if they do, check if polygon(s) large enough
            sh_pom = []
            
            if isinstance(sh, Polygon):
                if sh.area > min_area:      
                    sh_pom.append(sh)
                    
            else:
                for s in sh:
                    if s.area > min_area:      
                        sh_pom.append(s)
        
            ret_sh.append(sh_pom)
            ret_attr.append(attr)
            
            '''
            if attr['NUTS_ID'] == 'NO021':
                print(attr['NUTS_ID'], sh_pom)
            '''
            
    return ret_sh, ret_attr



def intersections(shapes, attributes, raster_data, connections, field='NUTS_ID',
                  raster_fname=None, unit_factor=10**-3, allowed_unallocated=0.05,
                  never_connect=None, largest_distance_proximity=2 * 10**3):
    '''
    For each rid cell find overlaping polygons and split emission according to share of overlaping.
    
    Some grid cells may be above water, outside of shapefile, but still within administrative domain. 
    These are detected by the fact most or all of their emissions do not get allocated, and 
    if such a cell is closer than largest_distance_proximity from nearest polygon, remaining emissions are allocated to that polygon.
    
    Function uses dictionary connections, which keeps track of which grid cells have allready been connected to polygons. 
    In case of FUA that only covers smaller part of area, never_connect set keeps track of all cells that 
    are not inside any of the polygons.
    
    

    Parameters
    ----------
    shapes : list
        polygons and multipolygons
    attributes : list
        attributes of polygons and multipolygons
    raster_data : list[list, numpy.ndarray, tuple]
        list of cell data - [coordinates of corners, emissions amount, (j, i) tuple]
        (j, i) tuple are indices of cell in array and are used as uniq ID of cells.
    connections : dict
        Dictionary keeping track of how much already visited (in previous runs or rasters) cells
        are covered by which polygons. 
    field : str, optional
        What is attribute name that uniquly identifies shapefiule polygons. For NUTS it is NUTS_ID,
        for FUA URAU_CODE. The default is 'NUTS_ID'.
    raster_fname : str, optional
        Name of raster file. It is only used to print it out in case too many emissions are left unallocated.
        The default is None.
    unit_factor : float, optional
        Factor to go from source units to kton. If source is in tons, unit_factor=10**-3. The default is 10**-3.
    allowed_unallocated : float, optional
        What share of emissions can stay unallocated (unless it is Industry, due to offshore emissions) 
        before execution is terminated and user prompted to handle the issie. The default is 0.05.
    never_connect : set, optional
        Primarily for FUA, to keep track of cell ids that are not within any polygon. The default is None.
    largest_distance_proximity : int, optional
        What is largest allowed distance from cell center to nearest polygon, to allocate emissions to that polygon. 
        It is only used if more than 0.5% of cell emissions did not get allocated. Hapens in fjords, for instance.
        The default is 2 * 10**3.

    Returns
    -------
    ret : dict
        For each administrative unit, sum of emissions.
    connections : dict
        For each cell id dictionary of intersecting polygons with shares - {(0, 0): {'NO053': 0.3, 'NO052': 0.7}}
    never_connect : set
        Set of cell ids that are not within any polygon. 

    '''
    ret = {}
    all_amount = 0
    allocated = 0
    for rd in raster_data:
        con_sw = False
        bc = rd[0]
        value = rd[1] * unit_factor
        cell_id = rd[2]
       
        if never_connect:
            if cell_id in never_connect:
                continue # must be continue - to skip remaining code in loop, break would jump out of current loop.
        
        all_amount += value
        loct = 0
        if cell_id in connections:
            con_sw = True
            conn = connections[cell_id]
            for (key, ratio) in conn.items():
                if key not in ret:
                    ret[key] = 0
                ret[key] += value * ratio
                allocated += value * ratio
                loct += value * ratio
        else:
            center = [0.5 * (bc[0][0] + bc[2][0]), 0.5 * (bc[0][1] + bc[2][1])]
            point = Point(center[0], center[1])
            poly_box = box(bc[0][0], bc[0][1], bc[2][0], bc[2][1])
            s0 = poly_box.area
            
            dist_min = 10**6
            sid = None
            for (sh, attr) in zip(shapes, attributes):
                
                for i, elem in enumerate(sh):
                    # distnace 0 if cell within polygon; find closest polygon - used for those at sea close to land
                    if dist_min > 0:
                        dist = point.distance(elem)
                        if dist < dist_min:
                            dist_min = dist
                            sid = attr[field]
                        
                    if poly_box.intersects(elem):
                        con_sw = True
                        intrsct = poly_box.intersection(elem)
                        inter_s = intrsct.area
                        ratio = inter_s / s0
                        
                        shp_id = attr[field]
                        if cell_id not in connections:
                            connections[cell_id] = {}
                        
                        if shp_id not in connections[cell_id]:
                            connections[cell_id][shp_id] = 0
                            
                        connections[cell_id][shp_id] += ratio
                        #print(intrsct, inter_s, ratio, value, attr)
                        
                        if attr[field] not in ret:
                            ret[shp_id] = 0
                        
                        ret[shp_id] += value * ratio
                        
                        allocated += value * ratio
                        
                        loct += value * ratio
                        
                    if loct >=value * 0.998:
                        #print('Found 99.8 percent, going to next grid')
                        break
                    
            if loct < value * 0.995:
                #print('Unallocated, using closeness', loct, value, dist_min, sid)
                if dist_min < largest_distance_proximity:
                    if sid not in ret:
                        ret[sid] = 0
                    ret[sid] += value - loct
                    allocated += value - loct
                    
                    if cell_id not in connections:
                        connections[cell_id] = {}
                    
                    if sid not in connections[cell_id]:
                        connections[cell_id][sid] = 0
                        
                    connections[cell_id][sid] += (value - loct) / value # by proximity
                    
  
        # never_connect is a set and current cell didnt connec to any polygon
        if never_connect is not None and not con_sw:
            never_connect.add(cell_id)
            
    
    if all_amount:
        print('All amount in grid (within polygons)', round(all_amount, 6), 'Allocated', round(allocated, 6), 'Amount left', round(all_amount - allocated, 6), 'Relative left [%]', round((all_amount - allocated) / all_amount * 100, 3))
        
        if (all_amount - allocated) / all_amount > allowed_unallocated:
            if 'Industry' not in raster_fname:
                print('Difference in all emissions in grid and allocated to polygons exceeds 5%!!!!! Almost certainly something wrong. Exiting')
                print('\n EXITING; Filename:', raster_fname)
                sys.exit()
            else:
                print('File is for industrial sources, allowing due to offshore emissions', raster_fname)
    else:
        print('No emissions')
    return ret, connections, never_connect
        

def make_shape_id_dictionary(attributes, shape_id_field='NUTS_ID', 
                             other_fields=['CNTR_CODE', 'NAME_LATN', 'NUTS_NAME']):
    '''
    Make dictionary gathering other fields at some shape_id_field. Used when writting result files.
    '''
    ret = {}
    for attr in attributes:
        ret[attr[shape_id_field]] = {of: attr[of] for of in other_fields}
    return ret

def write_result_file(directory, fname, data_dict, attributes_dict, year, pollutant, sector,
                      fields):
    '''
    Write csv result file.

    Parameters
    ----------
    directory : pathlib.Path
        Folder where to store file
    fname : name of file
        
    data_dict : dct
        result of function intersections(...)
    attributes_dict : dct
        Dictionary connecting the other attributes to each administrative unit, 
        made by make_shape_id_dictionary(...)
    year : int
        DESCRIPTION.
    pollutant : str
        DESCRIPTION.
    sector : str
        GNF sector
    fields : list
        Fields to write to csv file. Keep them as they are in codes using this,
        attributes_dict must have them for each administrative unit.
        Different for NUTS and FUA.

    Returns
    -------
    None.

    '''
    
    with open(directory / fname, 'w') as fil:
        fil.write(','.join(fields))
        fil.write('\n')
        for sid, amount in data_dict.items():
            fil.write(','.join([sid, attributes_dict[sid][fields[1]], attributes_dict[sid][fields[2]],
                                attributes_dict[sid][fields[3]], pollutant, 
                                str(year), sector, str(round(amount, 6))]))
            fil.write('\n')
    

def read_raster(file_path, file_name):
    '''
    Read raster file and prepeare values in functional way.
    
    Only uses grids cells with non zero values.

    Parameters
    ----------
    file_path : pathlib.Path
        To raster file
    file_name : str
        File name of raster file.

    Returns
    -------
    procd_data : list
        Each element in list has for all cells with non zero emissions:
            list of all four corners, emission value, tuple of indices in array. 

    '''

    raster = rasterio.open(file_path / file_name)

    transf = raster.transform

    x0 = transf.c
    y0 = transf.f

    dx = transf.a
    dy = transf.e
    
    '''
    bounds = raster.bounds
    print('bounds', bounds)
    
    point_min = (bounds.left, bounds.top)
    point_max = (bounds.right, bounds.bottom)
    '''
   
    
    data = raster.read(1)
    procd_data = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i, j]:
                ll = [x0 + j * dx, y0 + i * dy]
                ul = [x0 + j * dx, y0 + (i + 1) * dy]
                lr = [x0 + (j + 1) * dx, y0 + i * dy]
                ur = [x0 + (j + 1) * dx, y0 + (i + 1) * dy]
                
                #print([ll, lr, ur, ul])
                procd_data.append([[ll, lr, ur, ul], data[i, j], (i, j)])
    
    print('Number of cells', len(procd_data), 'Sum of data', np.sum(data), 'in source units')
    return procd_data




