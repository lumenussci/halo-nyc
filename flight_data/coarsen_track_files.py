import netCDF4 as nc
from h5py import File
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import pandas as pd
import glob,sys,pdb
from geopy.distance import geodesic

try:
    dx = float(sys.argv[1])
    halo_files = sorted(sys.argv[2:])
except:
    print('Usage: python coarse_track_files.py dx file_1 file_2 ...')
    sys.exit()

delt = np.max((1,int(dx/220)))#dt.timedelta(seconds=dx/220) #spatial averaging assumes a 220 m/s flight speed
print(delt)

t_lims = [[13.55,16.11],[18.2,20.75],[13.55,16.07],[18.2,20.75],[13.7,16.85],[16.98,20.15]]
for ifi,fi in zip(range(6),halo_files[:]):
    lon_avg,lat_avg,t_avg,xch4_avg = [],[],[],[]
    fid = nc.Dataset(fi,'r')
    xch4 = fid['CH4DataProducts']['XCH4_clear'][:].data.flatten()
    good_inds = range(len(xch4))#np.where(1-np.isnan(xch4))
    t = fid['Nav_Data']['gps_time'][:][good_inds].data.flatten()
    lat = fid['Nav_Data']['gps_lat'][:][good_inds].data.flatten()
    lon = fid['Nav_Data']['gps_lon'][:][good_inds].data.flatten()
    xch4 = xch4[good_inds]
    #pdb.set_trace()
    d = np.array([dt.datetime(2023,7,1) + dt.timedelta(seconds=ti*3600) for ti in t])
    df = pd.DataFrame({'t':t,'lat':lat,'lon':lon,'xch4':xch4},index=d)
    df_s = df.rolling(pd.Timedelta(seconds=dx/220)).mean()
    
    inds = np.where((t > t_lims[ifi][0])*(t < t_lims[ifi][1]))[0]
    tlat = df_s.iloc[inds]['lat']
    tlon = df_s.iloc[inds]['lon']
    tt = df_s.iloc[inds]['t']
    tch4 = df_s.iloc[inds]['xch4']
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

    fname = fi.split('.h5')[0]
    f = File(f'{fname}_{dx}m.h5','w')
    f.create_dataset('lat',data=lat_avg[:])
    f.create_dataset('lon',data=lon_avg[:])
    f.create_dataset('time',data=t_avg[:])
    f.create_dataset('xch4',data=xch4_avg[:])
    
    f.close()
