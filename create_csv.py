import pandas as pd
from h5py import File
import numpy as np
import datetime as dt
import glob

files = glob.glob('*300m*.h5')
for fi in files:
    flt = fi.split('_subCH4')[0]
    f = File(fi,'r')
    lat = f['lat'][:]
    lon = f['lon'][:]
    t = f['time'][:]
    obid = [f'{flt}_{i:05d}' for i in range(1,len(lat)+1)]
    #alt = np.array([max((hi,0)) for hi in h])

    d = fi.split('_')[0]
    datestr = f' {d[:4]}-{d[4:6]}-{d[6:8]}'
    nstep = len(t)
    i = 0
    while i < len(t):
        p2 = pd.DataFrame(columns=['obid',' lati',' long',' zagl',' UTC_date',' UTC_time'])
        p2['obid'] = obid[i:i+nstep]
        p2[' lati'] = lat[i:i+nstep]
        p2[' long'] = lon[i:i+nstep]
        p2[' zagl'] = np.zeros(nstep)
        p2[' UTC_date'][:] = [datestr for i in range(nstep)]
        d = np.array([dt.datetime(1970,1,1) + 3600*dt.timedelta(seconds=ti) for ti in t[i:i+nstep]])
        hrs = np.array([' '+d[i].strftime('%H:%M:%S') for i in range(len(d))])
        p2[' UTC_time'][:] = hrs[:]
        p2.round(4).to_csv(fi.split('.h5')[0]+'.csv',index=False)
        i += nstep
