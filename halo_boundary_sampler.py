"""
 ## halo_boundary_sampler.py

    ### purpose: samples wrfout files to pull out boundary conditions for STILT trajectories for the HALO data
    ### author: Sean Crowell
    ### input:
            - wrf_domain: string like d01 or d02
            - receptor_altitude: receptor altitude to distinguish file names
            - flts: flights to read in like 20230726_F1 - folder names
            - domain_box: either the string "wrf" to use the full WRF domain boundary or a lat/lon box
            of the form [lon lower bound, lon upper bound, lat lower bound, lat upper bound]
    ### output: 
            - trajectories inside the domain box as HDF5 file and boundary point HDF5 file with wrfout samples
"""

import os,sys,time
import glob,pdb
import importlib.util
from turfpy.measurement import boolean_point_in_polygon
from turfpy.transformation import intersect
from geojson import Point, Polygon, Feature
from h5py import File
import numpy as np
import datetime as dt

sys.path.append('../wrf-stilt')
wrf_stilt_utils = __import__("wrf-stilt-utils")

def create_halo_rds_arg_list(args_list=[],flts=[],receptor_altitude=50,bbox=None):
    """
    This function returns a list of dictionaries with the information needed to run the sampler routine.
    - Find all the RDS files in the STILT output directory
    - Screen those files for lats/lons that fall within the bounding box bbox (defined in main)
    - Create dictionary with ID, RDS file path, bbox, save directory for trajectory files
    - Return list of dictionaries with arguments
    """

    for flt in flts:
        flt_files = glob.glob(f'/scratch/07351/tg866507/halo/out_300m/{flt}*_{receptor_altitude}.0/particles/*.rds')
        lats = [fi.split('/')[-1].split('_')[2] for fi in flt_files]
        lons = [fi.split('/')[-1].split('_')[1] for fi in flt_files]
        points = [Feature(geometry=Point([float(lat),float(lon)])) for (lat,lon) in zip(lats,lons)]
        in_domain = np.where(list(map(boolean_point_in_polygon,points,[bbox for i in range(len(flt_files))])))[0]
        for fi in [flt_files[i] for i in in_domain]:
            args_list.append({'ID':fi.split('/')[-1].split('_traj.rds')[0],'trajectory_rds_filename':fi,'bbox':bbox,'save_dir':f'trj/{flt}/{fi.split("/")[-1].split('_')[-2]}/','write_files':True})
    return args_list

if __name__ == "__main__":
    """
    In order to run the sampler functions in wrf-stilt-utils.py:
        1. create bounding box for domain
        2. create list of RDS files with domain screening
        3. define a list of variables desired from the receptor location/time
        4. create a list of dictionaries with all of the variables needed by the sampler
    What is returned is the set of boundary value file paths that were created by the sampler 
    """
    try:
        wrf_domain = sys.argv[1]     # WRF domain: d01 or some higher resolution domain
        receptor_altitude = sys.argv[2]   # Altitude of receptor used to subset filenames
        flts = sys.argv[3:]
        domain_box = "wrf"#[-75.2,-72,40,42] #Either a lat/lon box or the boundaries of the WRF domain
    except IndexError:
        print("Usage: python boundary_condition_pipeline.py domain altitude flight1 flight2 ...")
        sys.exit()

    #WRF Domain
    files = glob.glob(f'/work2/07655/tg869546/stampede3/nyc-chem/2023/wrfout/*{wrf_domain}*')
    wrf_f = File(files[0])
    wrf_lat,wrf_lon = wrf_f['XLAT'][:][0],wrf_f['XLONG'][:][0]
    poly_verts_wrf = [(float(wrf_lat[0,0]),float(wrf_lon[0,0])),(float(wrf_lat[-1,0]),float(wrf_lon[-1,0])),(float(wrf_lat[-1,-1]),float(wrf_lon[-1,-1])),(float(wrf_lat[0,-1]),float(wrf_lon[0,-1]))]
    
    #STILT Domain
    fp_file = File(glob.glob(f'/scratch/07351/tg866507/halo/out_300m/col_footprints/20230726_F1/*foot.nc')[0],'r')
    fp_lat = fp_file['lat'][:]
    fp_lon = fp_file['lon'][:]
    lon_lb,lon_ub,lat_lb,lat_ub = fp_lon.min(),fp_lon.max(),fp_lat.min(),fp_lat.max()
    poly_verts_stilt = [(float(lat_lb),float(lon_lb)),(float(lat_ub),float(lon_lb)),(float(lat_ub),float(lon_ub)),(float(lat_lb),float(lon_ub))]
    
    bbox_wrf = Feature(geometry=Polygon((poly_verts_wrf,)))
    bbox_stilt = Feature(geometry=Polygon((poly_verts_stilt,)))
    bbox = intersect([bbox_wrf,bbox_stilt])

    for flt in flts:
        rds_args_list = []
        rds_args_list = create_halo_rds_arg_list(rds_args_list,flts=[flt],bbox=bbox,receptor_altitude=receptor_altitude)
        start_time = time.time()
        receptor_wrf_sample_vars = ['pressure','psfc','CO2_BCK','CO2_ANT','CO2_BIO','CH4_BCK','CH4_ANT']
        args_list = [
            {
            'trajectory_rds_filename':arg['trajectory_rds_filename'], \
            'ID':arg['ID'], \
            'bbox':bbox, \
            'rec_sample_vars':receptor_wrf_sample_vars, \
            'wrf_path':'/work2/07655/tg869546/stampede3/nyc-chem/2023/wrfout/', \
            'wrf_domain':wrf_domain,
            'trj_save_dir':arg['save_dir'], \
            'bnd_save_dir':f'./bnd/{arg["ID"][-17:-6]}/wrf_{wrf_domain}/{arg["ID"].split("_")[-1]}/', \
            'overwrite_bnd':True, \
            'overwrite_trj':True \
             } \
             for arg in rds_args_list[:]]

        results_dict = wrf_stilt_utils.run_function_in_parallel(wrf_stilt_utils.stilt_boundary_wrfout_sampler,args_list[:])
#    results_dict = wrf_stilt_utils.stilt_boundary_wrfout_sampler(**args_list[0])
    print(f'Boundary Sampler Execution Time: {time.time()-start_time} seconds')


