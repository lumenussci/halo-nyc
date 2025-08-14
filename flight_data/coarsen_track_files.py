import netCDF4 as nc
from h5py import File
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import glob,sys
from geopy.distance import geodesic

try:
    halo_files = sorted(sys.argv[1:])
except:
    print('Usage: python coarse_track_files.py file_1 file_2 ...')
    sys.exit()
    
t_lims = [[13.55,16.11],[18.2,20.75],[13.55,16.07],[18.2,20.75],[13.7,16.85],[16.98,20.15]]
for ifi,fi in zip(range(6),halo_files[:]):
    lon_avg,lat_avg,t_avg,xch4_avg = [],[],[],[]
    fid = nc.Dataset(fi,'r')
    xch4 = fid['CH4DataProducts']['XCH4_clear'][:]
    good_inds = range(len(xch4))#np.where(1-np.isnan(xch4))
    t = fid['Nav_Data']['gps_time'][:][good_inds]
    lat = fid['Nav_Data']['gps_lat'][:][good_inds]
    lon = fid['Nav_Data']['gps_lon'][:][good_inds]
    xch4 = xch4[good_inds]
    
    inds = np.where((t > t_lims[ifi][0])*(t < t_lims[ifi][1]))[0]
    tlat = lat[inds]
    tlon = lon[inds]
    tt = t[inds]
    tch4 = xch4[inds]
    ii = 0
    while ii < len(inds)-1:
        ub = min((200,len(inds)-ii-1))
        dd = np.array([geodesic((tlat[ii],tlon[ii]),(tlat[ii+jj],tlon[ii+jj])).m for jj in range(1,ub)])
        jnds = [ii]
        jnds.extend(ii+np.where(dd < dx)[0])
        lon_avg.append(tlon[jnds].mean())
        lat_avg.append(tlat[jnds].mean())
        t_avg.append(tt[jnds].mean())
        xch4_avg.append(np.nanmean(tch4[jnds]))
        ii += np.max((1,len(jnds)))

    t_avg = np.array(t_avg)
    lat_avg = np.array(lat_avg)
    lon_avg = np.array(lon_avg)
    xch4_avg = np.array(xch4_avg)

    fname = fi.split('/')[-1].split('.h5')[0]
    f = File(f'{fname}_{dx}m.h5','w')
    f.create_dataset('lat',data=lat_avg[:])
    f.create_dataset('lon',data=lon_avg[:])
    f.create_dataset('time',data=t_avg[:])
    f.create_dataset('xch4',data=xch4_avg[:])
    
    f.close()