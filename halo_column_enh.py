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

def create_halo_file_list(wrf_domain='d01',flts=[],fpt_dir='./fpt/'):
    """
    This function returns a list of dictionaries with the information needed to compute the column average footprints.
    """

    file_list = []
    for flt in flts:
        fpt_files = sorted(glob.glob(f'{fpt_dir}/{flt}/*'))
        ob_num = np.array([int(fi.split('/')[-1].split('_')[-2]) for fi in fpt_files])
        file_list.append({'ID':f'{flt}','fp_files':fpt_files,'ob_num':ob_num})
    return file_list

def compute_pressure_wgt_avg_enh(ID=None,nhours=0,emis={},fp_files=[],ob_num=[],col_enh_save_dir='./'):
    """
    Takes in a list of files for a single latitude and longitude and computes the column average footprints
    """

    # Not all footprint files have the same number of time steps
    with File(fp_files[0],'r') as fp_f:
        fp_f = File(fp_files[0],'r')
        fp_lat = fp_f['lat'][:]
        fp_lon = fp_f['lon'][:]
    
    inds = np.argsort(ob_num)
    fpt_files_sort = np.array(fp_files)[inds]
    ob_num_sort = np.array(ob_num)[inds]
    dxch4 = {}
    for ky in list(emis.keys()):
        dxch4[ky] = np.zeros((len(fpt_files_sort),emis[ky].shape[0]))
    for ifi,fi in enumerate(fpt_files_sort[:]):
        with File(fi,'r') as fp_f:
            fnh = fp_f['pwcol_fp'][:].shape[0]
            ih = np.min((fnh,nhours))
            #fpt[ifi,-ih:] = fp_f['pwcol_fp'][:]
            ft = fp_f['pwcol_fp'][:]
            for ky in list(emis.keys()):
                dxch4[ky][ifi,:] = np.array([(emis[ky][i_cat]*ft).sum() for i_cat in range(emis[ky].shape[0])])
    # emis is a dictionary of dictionaries with the individual emission sources in it
    # These could have time and category axes - so we have to standarize
    # emis[emis_name] = emis_array, with emis_array.shape = [n_cats, n_times, n_lat, n_lon]
    dxch4['ID'] = ID 
    dxch4['ob_num'] = ob_num
    
    with File(f'{col_enh_save_dir}/{ID}_dxch4.h5','w') as f:
        for ky in emis.keys():
            f.create_dataset(ky,data=dxch4[ky][:])
        f.create_dataset('ob_num',data=ob_num[:])
        f.close()
    return dxch4 


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

    fp_file = File(glob.glob(f'/scratch/07351/tg866507/halo/out_300m/col_footprints/20230726_F1/*foot.nc')[0],'r')
    fp_lat = fp_file['lat'][:]
    fp_lon = fp_file['lon'][:]

    emis = {}
    emis_names = ['pitt','epa','edgar']
    with File('/scratch/07351/tg866507/halo-staaqs/nyc_ch4_emissions.h5','r') as f:
        for nm in emis_names:
            emis[nm] = {}
            e_lat = f['lat'][:]
            e_lon = f['lon'][:]
            fch4 = f[nm][:]
            #emis[nm]['fch4'] = fch4[:]#np.array([RegularGridInterpolator((e_lat,e_lon),fch4[i],method='nearest')(fp_lat,fp_lon) for i in range(fch4.shape[0])])
            emis[nm] = fch4[:]
            emis[nm+'_categories'] = f.attrs[f'{nm}_categories'].split(';')

    emis_regrid = {}
    with File('/scratch/07351/tg866507/halo-staaqs/nyc_emissions_regrid.h5','r') as f:
        for nm in emis_names:
            emis_regrid[nm] = {}
            fp_lat = f['lat'][:]
            fp_lon = f['lon'][:]
            emis_regrid[nm] = f[nm][:]
            emis_regrid[nm+'_categories'] = f.attrs[f'{nm}_categories'].split(';')

    file_list = create_halo_file_list(flts=flts,wrf_domain=wrf_domain,fpt_dir=f'../halo/out_300m/col_footprints/')
    args_list = [
            {'ID':fl['ID'],'fp_files':fl['fp_files'],'ob_num':fl['ob_num'],'nhours':24,'emis':{'edgar':emis_regrid['edgar'][:,None,:,:],
                'epa':emis_regrid['epa'][:,None,:,:],'pitt':emis_regrid['pitt'][:,None,:,:]}} for fl in file_list]

    start_time = time.time()
    results_dict = wrf_stilt_utils.run_function_in_parallel(compute_pressure_wgt_avg_enh,args_list[:])
    #results_dict = compute_pressure_wgt_avg_enh(**args_list[0])

    print(f'Column Avg Footprint Computation Complete after {time.time()-start_time} seconds!')
    
