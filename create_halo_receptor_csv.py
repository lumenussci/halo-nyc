import pandas as pd
from h5py import File
import numpy as np
import datetime as dt
from glob import glob
import sys

try:
    flight_file = sys.argv[1]
    nstep = sys.argv[2]
    out_dir = sys.argv[3]
except IndexError:
    print("Usage: python create_halo_receptor_csv.py flight_day_csv n_receptors_per_outfile save_dir")
    sys.exit()

p = pd.read_csv('../nyc_cf/receptors.csv')

nstep = 36
with p as pd.read_csv(flight_file):
    datestr = fi.split('/')[1][:11]
    i = 0
    while i < len(p[' lati'][:]):
        p2 = pd.DataFrame(columns=p.keys())
        p2['obid'] = p['obid'][i+nstep]
        p2[' lati'] = p[' lati'][i:i+nstep]
        p2[' long'] = p[' long'][i:i+nstep]
        p2[' zagl'] = p[' zagl'][i:i+nstep]
        p2[' UTC_date'][:] = p[' UTC_date'][i:i+nstep]
        p2[' UTC_time'][:] = p[' UTC_time'][i:i+nstep]
        p2.round(4).to_csv('receptors/%s_receptors_%04d_%04d.csv'%(datestr,i,i+len(p2['lati'])-1),index=False)
        i += nstep
