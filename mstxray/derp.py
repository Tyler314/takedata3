import os
import datetime
import h5py
import numpy as np

from mst import mdsplus as mds
from calibration import Calibration
import constants as const
import physical as mp
import dbutil

TD3_DB_XRAY_PATH = '/home/xray/takedata3/mstxray/db/derp.sqlite' #xRayTD3.sqlite'

shot = '1'
query  = 'SELECT * FROM xRayData WHERE sn=?'
       
row = dbutil.fetchall(TD3_DB_XRAY_PATH, query, (shot))[0]

print row        
