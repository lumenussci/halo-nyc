"""
 ## emissions_resampler.py

    ### purpose: regrids emissions
    ### author: Sean Crowell
    ### input:
            - input file to be regridded 
            - output file with lat/lon arrays
    ### output: 
            - output file augmented with emissions
"""

import os,sys,time
import glob,pdb
import importlib.util
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, MultiPolygon, Feature
from h5py import File
import numpy as np
import datetime as dt
import netCDF4 as nc
from scipy.interpolate import RectBivariateSpline,RegularGridInterpolator

sys.path.append('../wrf-stilt')
wrf_stilt_utils = __import__("wrf-stilt-utils")


def NN_regrid_emissions(ID=None,emis_in={'emis':None,'lat':None,'lon':None},emis_out={'emis':None,'lat':None,'lon_out':None,'loc_inds':None}):
    """
    Takes input and output emissions dictionaries and regrids the input emissions to a point or a set of points
    """
    fint = RegularGridInterpolator((emis_in['lat'],emis_in['lon']),emis_in['emis'],method='nearest')
    
    emis_out['emis'] = np.zeros((len(emis_out['lat']),len(emis_out['lon'])))
    loc_inds = emis_out['loc_inds']
    pts = np.zeros((loc_inds.shape[0],2))
    for i in range(loc_inds.shape[0]):
        pts[i] = [emis_out['lat'][loc_inds[i,0]],emis_out['lon'][loc_inds[i,1]]]
    eint = fint(pts)
    emis_int = np.zeros((fp_lat.shape[0],fp_lon.shape[0]))
    for i in range(len(loc_inds)):
        emis_int[*loc_inds[i,:]] = eint[i]

    emis_out['emis'] = emis_int
    return emis_out

def compute_xgas_enhancements(ID=None,emis={},footprint=None):
    
    emis_names = list(emis.keys())
    enh = dict.fromkeys(emis_names)
    for ky in emis_names:
        if len(emis[ky].shape) == 2:
            enh[ky] = (emis[ky]*footprint).sum()
        else:
            enh[ky] = (emis[ky]*footprint[None]).sum((1,2))
    return enh

if __name__ == "__main__":
    """
    What is returned is the set of boundary value file paths that were created by the sampler 
    """
    #try:
    #    emis_in_file = sys.argv[1]     # WRF domain: d01 or some higher resolution domain
    #    lat_str = sys.argv[2].split(',')
    #    lat_out = np.arange(np.float(lat_str[0]),np.float(lat_str[1]),np.float(lat_str[3]))
    #    lon_str = sys.argv[2].split(',')
    #    lon_out = np.arange(np.float(lon_str[0]),np.float(lon_str[1]),np.float(lon_str[3]))

    #except IndexError:
    #    print("Usage: python boundary_condition_pipeline.py domain altitude flight1 flight2 ...")
    #    sys.exit()

    emis_file = File('nyc_ch4_emissions.h5')
    emis = {}
    emis['lat'] = emis_file['lat'][:]
    emis['lon'] = emis_file['lon'][:]
    for v in ['edgar','epa','pitt']:
        emis[v] = emis_file[v][:]
        emis[f'{v}_cats'] = emis_file.attrs[f'{v}_categories'].split(';')
    flt='20230726_F1'
    fp_file = File(glob.glob(f'/scratch/07351/tg866507/halo/out_300m/col_footprints/{flt}/*foot.nc')[0],'r')
    fp_lat = fp_file['lat'][:]
    fp_lon = fp_file['lon'][:]
    loc_inds = np.zeros((len(fp_lat)*len(fp_lon),2)).astype(int)
    for ilat in range(len(fp_lat)):
        loc_inds[ilat*len(fp_lon):(ilat+1)*len(fp_lon),0] = ilat
        loc_inds[ilat*len(fp_lon):(ilat+1)*len(fp_lon),1] = range(len(fp_lon))

    args_list = []
    for iv,v in enumerate(['edgar','epa','pitt']):
        for ic,c in enumerate(emis[f'{v}_cats']):
            args_list.append(
                    {'ID':f'{v}: {c}','emis_in':{'emis':emis[v][ic,:,:],'lat':emis['lat'],'lon':emis['lon']},
                        'emis_out':{'emis':None,'lat':fp_lat,'lon':fp_lon,'loc_inds':loc_inds}}    
                )

    start_time = time.time()
    #results_dict = wrf_stilt_utils.run_function_in_parallel(NN_regrid_emissions,args_list[:])
    results_dict = NN_regrid_emissions(ID=args_list[1]['ID'],emis_in=args_list[1]['emis_in'],emis_out=args_list[1]['emis_out'])
    
    emis_regrid = {} 
    emis_source = np.array([args_list[i]['ID'].split(':')[0] for i in range(len(args_list))])
    emis_cats = np.array([args_list[i]['ID'].split(':')[1] for i in range(len(args_list))])
    emis_source_cats_str = {}
    for iv,v in enumerate(['edgar','epa','pitt']):
        source_inds = np.where(emis_source == v)
        emis_source_cats = emis_cats[source_inds]
        emis_regrid[v] = np.zeros((len(emis_source_cats),fp_lat.shape[0],fp_lon.shape[0]))
        for ic,c in enumerate(emis[f'{v}_cats']):
            source_cat_inds = np.where(emis_source_cats == c)
            emis_regrid[v][ic] = results_dict[f'{v}: {c}']['result']['emis'][:]
        emis_source_cats_str[v] = ';'.join(emis[f'{v}_cats'])

    with File('./nyc_emissions_regrid.h5','w') as f:
        for iv,v in enumerate(['edgar','epa','pitt']):
            f.create_dataset(v,data=emis_regrid[v][:])
            f.attrs[f'{v}_categories'] = emis_source_cats_str[v]
        f.create_dataset('lat',data=fp_lat[:])
        f.create_dataset('lon',data=fp_lon[:])
 
