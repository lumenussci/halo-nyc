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

def compute_xgas_enhancements(ID=None,emis={},footprint=None):
    
    emis_names = list(emis.keys())
    enh = dict.fromkeys(emis_names)
    for ky in emis_names:
        if len(emis[ky].shape) == 2:
            enh[ky] = (emis[ky]*footprint).sum()
        else:
            enh[ky] = (emis[ky]*footprint[None]).sum((1,2))
    return enh

def read_footprint_files(ID=None,flt=None,fp_save_dir='./'):

    fp_files = glob.glob(f'{fp_save_dir}/{flt}/*col*')
    obs_nums = np.array([int(fi.split('/')[-1].split('_')[-2]) for fi in fp_files])
    inds = np.argsort(obs_nums)
    fpt = File(fp_files[0],'r')['pwcol_fp'][:].sum(0)
    flt_fpt = np.zeros((len(obs_nums),fpt.shape[0],fpt.shape[1]))
    for iob,ob in enumerate(inds):
        flt_fpt[iob] = File(fp_files[ob],'r')['pwcol_fp'][:].sum(0)
    return flt_fpt

if __name__ == "__main__":

    emis = {}
    emis_names = ['edgar','epa','pitt']
    with File('./nyc_emissions_regrid.h5','r') as f:
        for iv,v in enumerate(emis_names):
            emis[v] = f[v][:]
            emis[f'{v}_categories'] = f.attrs[f'{v}_categories']  
    flts = ['20230726_F1','20230726_F2','20230728_F1','20230728_F2','20230805_F1','20230809_F1']
    fpt = {}
    enh = {}
    args_list = [{'ID':flt,'flt':flt,'fp_save_dir':'/scratch/07351/tg866507/halo/out_300m/col_footprints/'} for flt in flts]
    results_dict = wrf_stilt_utils.run_function_in_parallel(read_footprint_files,args_list[:])
        fpt[flt] = read_footprint_files(flt=flt,fp_save_dir='/scratch/07351/tg866507/halo/out_300m/col_footprints/')
        
