"""
 ## halo_column_xgas_bc.py

    ### purpose: Creates column xgas boundary conditions for the HALO receptors
    ### author: Sean Crowell
    ### input:
            - wrf_domain: string like d01 or d02
            - flts: flights to read in like 20230726_F1 - folder names
    ### output: 
            - netCDF4 files with column concentrations 
"""

import os,sys,time
import glob,pdb
import importlib.util
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, MultiPolygon, Feature
from h5py import File
import numpy as np
import datetime as dt

sys.path.append('../wrf-stilt')
wrf_stilt_utils = __import__("wrf-stilt-utils")

def create_halo_bnd_list(wrf_domain='d01',flts=[],bnd_dir='./'):
    """
    This function returns a list of dictionaries with the information needed to compute the column xgas values.
    """

    args_list = []
    for flt in flts:
        bnd_folders = sorted(glob.glob(f'{bnd_dir}/{flt}/*'))
        for ob_num in bnd_folders[:]:
            levs = sorted([int(fi.split('/')[-1].split('_')[3]) for fi in glob.glob(f'{ob_num}/*wrf{wrf_domain}.h5')])
            files = {}
            for lev in levs:
                files[lev] = glob.glob(f'{ob_num}/*_{lev}_*wrf{wrf_domain}.h5')[0]
            args_list.append({'ID':f'{flt}_{ob_num.split('/')[-1]}','files':files})
    return args_list

def compute_pressure_wgt_avg(ID=None,files=[],var_dict={}):
    """
    Takes in a list of files for a single latitude and longitude and computes the column average values
    for any variable in the keys of var_dict
    """
    
    levs = sorted(list(files.keys()))
    var_names = list(var_dict.keys())
    tdict = {}
    for v in var_names:
        tdict[v] = np.zeros(len(levs))
    tdict['pressure'] = np.zeros(len(levs))

    for iz,z in enumerate(levs):
        f = File(files[z])
        for v in var_names:
            gp = var_dict[v]['bnd_group']
            if f[gp][v].shape[0] > 1:
                tdict[v][iz] = np.median(f[gp][v][:].flatten())
            else:
                tdict[v][iz] = f[gp][v][:].flatten()[0]
        tdict['pressure'][iz] = f['receptor']['pressure'][:].flatten()[0]
    tdict['psfc'] = f['receptor']['psfc'][:].flatten()[0]
    tdict['lat'] = f['receptor']['lat'][:].flatten()[0]
    tdict['lon'] = f['receptor']['lon'][:].flatten()[0]
    tdict['time'] = f['receptor']['time'][:].flatten()[0]

    #Compute the column values
    p_inds = np.argsort(tdict['pressure'])
    p_ext = [tdict['psfc']]
    p_ext.extend(tdict['pressure'][p_inds[::-1]])
    p_sort = np.array(p_ext)
    dp = []
    dp.extend(-np.diff(p_sort))
    dp = np.array(dp)
    for v in var_names:
        tdict[f'x{v}'] = np.dot(tdict[v],dp)/dp.sum()
    
    return_dict = {}
    for v in var_names:
        return_dict[f'x{v}'] = tdict[f'x{v}']
    for v in ['lat','lon','time','psfc']:
        return_dict[v] = tdict[v]
    
    return return_dict


if __name__ == "__main__":
    """
    What is returned is the set of boundary value file paths that were created by the sampler 
    """
    try:
        wrf_domain = sys.argv[1]     # WRF domain: d01 or some higher resolution domain
        flts = sys.argv[2:]
    except IndexError:
        print("Usage: python boundary_condition_pipeline.py domain altitude flight1 flight2 ...")
        sys.exit()

    file_list = create_halo_bnd_list(flts=flts,wrf_domain=wrf_domain,bnd_dir=f'./bnd/wrf_{wrf_domain}/')
    args_list = [
            {'ID':fl['ID'],'files':fl['files'],'var_dict':{'wrf_bc_co2':{'bnd_group':'boundary'},'wrf_bc_ch4':{'bnd_group':'boundary'}}} for fl in file_list]

    start_time = time.time()
    results_dict = wrf_stilt_utils.run_function_in_parallel(compute_pressure_wgt_avg,args_list[:])
    kys = np.array(sorted(list(results_dict.keys())))
    flt_kys = np.array([ri[:11] for ri in kys])
    halo_out = {}
    for flt in flts:
        halo_out[flt] = {}
        inds = np.where(flt_kys == flt)[0]
        for v in ['xwrf_bc_co2','xwrf_bc_ch4','lat','lon','time','psfc']:
            halo_out[flt][v] = np.array([results_dict[ky]['result'][v] for ky in kys[inds].tolist()])
        halo_out[flt]['obs_id'] = np.array([int(ky.split('_')[-1]) for ky in kys[inds].tolist()])
        
        with File(f'bnd/wrf_{wrf_domain}/{flt}_xbnd.h5','w') as f:
            for v in halo_out[flt].keys():
                f.create_dataset(v,data=halo_out[flt][v][:])
        
    print(f'Boundary Sampler Execution Time: {time.time()-start_time} seconds')

    
