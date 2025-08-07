"""
 ## sample_wrfout_ghg.py

    ### purpose: samples wrfout files to pull out boundary conditions for STILT trajectories for the HALO data
    ### author: Sean Crowell
    ### input:
            - domain: string like d01 or d02
            - altitude: receptor altitude to distinguish file names
            - flts: flights to read in like 20230726_F1 - folder names
    ### output: 
            - boundary point HDF5 file with wrfout samples
"""

import os,sys,time
import glob
import importlib.util
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, MultiPolygon, Feature
from h5py import File

sys.path.append('../wrf-stilt')
wrf_stilt_utils = __import__("wrf-stilt-utils")

def create_halo_rds_arg_list(args_list=[],flts=[],receptor_altitude=50,bbox=[-75.2,-72,40,42]):
    for flt in flts:
        flt_files = glob.glob(f'/scratch/07351/tg866507/halo/out_300m/{flt}*_{receptor_altitude}.0/particles/*.rds')
        for fi in flt_files:
            args_list.append({'ID':fi.split('/')[-1].split('_traj.rds')[0],'trajectory_rds_filename':fi,'bbox':bbox,'save_dir':f'trj/{flt}/','write_files':True})
    return args_list

if __name__ == "__main__":
    try:
        wrf_domain = sys.argv[1]     # WRF domain: d01 or some higher resolution domain
        receptor_altitude = sys.argv[2]   # Altitude of receptor used to subset filenames
        flts = sys.argv[3:]
        domain_box = "wrf"#[-75.2,-72,40,42] #Either a lat/lon box or the boundaries of the WRF domain
    except IndexError:
        print("Usage: python sample_wrfout.py domain altitude flight1 flight2 ...")
        sys.exit()

    files = glob.glob(f'/work2/07655/tg869546/stampede3/nyc-chem/2023/wrfout/*{wrf_domain}*')
    wrf_f = File(files[0])
    wrf_lat,wrf_lon = wrf_f['XLAT'][:][0],wrf_f['XLONG'][:][0]
    poly_verts = [(float(wrf_lat[0,0]),float(wrf_lon[0,0])),(float(wrf_lat[-1,0]),float(wrf_lon[-1,0])),(float(wrf_lat[-1,-1]),float(wrf_lon[-1,-1])),(float(wrf_lat[0,-1]),float(wrf_lon[0,-1]))]
    bbox = Feature(geometry=MultiPolygon([(poly_verts,)]))


    args_list = []
    for flt in flts:
        args_list = create_halo_rds_arg_list(args_list,flts=[flt],bbox=bbox,receptor_altitude=receptor_altitude)
    start_time = time.time()
    results_dict = wrf_stilt_utils.run_function_in_parallel(wrf_stilt_utils.locate_trajectory_boundary_points,args_list)
    print(f'Boundary Location Execution Time: {time.time()-start_time} seconds')

    start_time = time.time()
    # results_dict looks like results_dict['*.rds'] = {result,warnings,error}
    # result is the output of the function, in this case bnd = result
    # bnd is the boundary point information for each trajectory that is needed for the WRF sampler routines:
    #   bnd_lat: vector of boundary point latitudes (n_particles x 1)
    #   bnd_lon:                          longitudes
    #   bnd_alt:                          altitudes
    #   bnd_t:                            times (seconds after 1970-1-1)
    #   t:                                times (seconds since release time)
    #   obs_lat: latitude of receptor
    #   obs_lon: longitude  "" 
    #   obs_alt: zagl       ""
    #   obs_t:   time       "" (seconds after 1970-1-1)

    args_list_2 = []
    for ky in list(results_dict.keys()):
        ky_comp = ky.split('_')
        flt = '_'.join(ky_comp[-3:-1])
        dct = results_dict[ky]['result']
        dct['ID'] = ky.split('_traj.rds')[0]
        dct['wrf_path'] = '/work2/07655/tg869546/stampede3/nyc-chem/2023/wrfout/'
        dct['wrf_domain'] = wrf_domain
        for v in ['pressure','psfc','CO2_BCK','CO2_ANT','CO2_BIO','CH4_BCK','CH4_ANT']:
            dct['receptor_loc_vars'][v] = None
        dct['save_dir'] = f'./bnd/{flt}/{dct["ID"]}/'
        dct['overwrite'] = True
        args_list_2.append(dct)

    results_dict_2 = wrf_stilt_utils.run_function_in_parallel(wrf_stilt_utils.stilt_boundary_wrfout_sampler,args_list_2)
    print(f'Boundary Sampler Execution Time: {time.time()-start_time} seconds')
