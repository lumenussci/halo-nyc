"""
 ## halo_column_fp.py

    ### purpose: Creates column average footprints for the HALO receptors
    ### author: Sean Crowell
    ### input:
            - wrf_domain: string like d01 or d02
            - flts: flights to read in like 20230726_F1 - folder names
    ### output: 
            - netCDF4 files with column concentrations z
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

sys.path.append('../wrf-stilt')
wrf_stilt_utils = __import__("wrf-stilt-utils")

def create_halo_file_list(wrf_domain='d01',flts=[],bnd_dir='./bnd/',fpt_dir='./fpt/'):
    """
    This function returns a list of dictionaries with the information needed to compute the column average footprints.
    """

    file_list = []
    for flt in flts:
        bnd_folders = sorted(glob.glob(f'{bnd_dir}/wrf_{wrf_domain}/{flt}/*'))
        for ob_num_fol in bnd_folders[:]:
            ob_num = ob_num_fol.split('/')[-1]
            levs = sorted([int(fi.split('/')[-1].split('_')[3]) for fi in glob.glob(f'{ob_num_fol}/*wrf{wrf_domain}.h5')])
            bnd_files = {}
            fp_files = {}
            for lev in levs:
                bnd_files[lev] = glob.glob(f'{ob_num_fol}/*_{lev}_*wrf{wrf_domain}.h5')[0] 
                fp_file_list = glob.glob(f'/scratch/07351/tg866507/halo/out_300m/{flt}*_{lev}.0/footprints/*_{ob_num}*.nc')
                fp_files[lev] = None
                if len(fp_file_list) > 0:
                    fp_files[lev] = fp_file_list[0]
            file_list.append({'ID':f'{flt}_{ob_num}','bnd_files':bnd_files,'fp_files':fp_files})
    return file_list

def compute_pressure_wgt_avg_footprint(ID=None,bnd_files={},fp_files={},col_fp_save_dir='./'):
    """
    Takes in a list of files for a single latitude and longitude and computes the column average footprints
    """

    levs = sorted(list(bnd_files.keys()))
    nhrs = 0
    for lev in levs:
        if fp_files[lev] == None: continue
        nhrs = np.max((nhrs,nc.Dataset(fp_files[lev]).dimensions['time'].size))
        lat = nc.Dataset(fp_files[lev]).dimensions['lat'].size
        lon = nc.Dataset(fp_files[lev]).dimensions['lon'].size
        fp_filename = fp_files[lev]
    if nhrs == 0:
        return 'No valid footprints'

    fpts = np.zeros((len(levs),nhrs,lat,lon))
    for ilev,lev in enumerate(levs):
        if fp_files[lev] == None: continue
        lev_nhrs = nc.Dataset(fp_files[lev]).dimensions['time'].size
        fpts[ilev,-lev_nhrs:,:,:] = nc.Dataset(fp_files[lev])['foot'][:]
        latv = nc.Dataset(fp_files[lev])['lat'][:]
        lonv = nc.Dataset(fp_files[lev])['lon'][:]

    tdict={}
    tdict['pressure'] = np.zeros(len(levs))
    for iz,z in enumerate(levs):
        f = File(bnd_files[z])
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

    pwcol_fpt = (dp[:,None,None,None]*fpts).sum(0)/dp.sum()
    
    fn_comps = fp_filename.split('/')[-1].split('_')
    flt = '_'.join(fn_comps[4:6])
    col_fn_1 = '_'.join(fn_comps[:3])
    col_fn_2 = '_'.join(fn_comps[4:])
    col_fn = col_fp_save_dir+flt+'/'+col_fn_1+'_col_'+col_fn_2
    with File(col_fn,'w') as f:
        f.create_dataset('pwcol_fp',data=pwcol_fpt[:])
        f.create_dataset('lat',data=latv[:])
        f.create_dataset('lon',data=lonv[:])
        f.close()

    return_dict = {'dp':dp,'pwcol_fpt':pwcol_fpt,'fp_sum':fpts.sum((1,2,3))}
    for v in ['lat','lon','time']:
        return_dict[v] = tdict[v]
    
    return col_fn #return_dict


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

    file_list = create_halo_file_list(flts=flts,wrf_domain=wrf_domain,bnd_dir='./bnd/')
    args_list = [
            {'ID':fl['ID'],'bnd_files':fl['bnd_files'],'fp_files':fl['fp_files'],'col_fp_save_dir':f'/scratch/07351/tg866507/halo/out_300m/col_footprints/'} for fl in file_list]

    pdb.set_trace()
    start_time = time.time()
    results_dict = wrf_stilt_utils.run_function_in_parallel(compute_pressure_wgt_avg_footprint,args_list[:])
    #results_dict = compute_pressure_wgt_avg_footprint(ID=args_list[0]['ID'],bnd_files=args_list[0]['bnd_files'],fp_files=args_list[0]['fp_files'],col_fp_save_dir=args_list[0]['col_fp_save_dir'])

    print(f'Column Avg Footprint Computation Complete after {time.time()-start_time} seconds!')
    
