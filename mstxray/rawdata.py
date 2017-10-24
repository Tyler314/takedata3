from __future__ import division, print_function

import os
import datetime
import h5py
import numpy as np

from mst import mdsplus as mds
from calibration import Calibration
import constants as const
import physical as mp
import dbutil

TD3_DB_XRAY_PATH = '/home/xray/takedata3/mstxray/db/xRayTD3.sqlite'

def get_next_shot():
    """For saving non-mdsplus data, use this funtion to generate the next
    shot number to use. 
    """
    # Today's date.
    date_num = mds.date_to_date_num(datetime.date.today())
    
    # This is the directory the file will reside in. 
    dir_ = os.path.join(const.DATA_DIR, str(date_num))
    
    shot = date_num * 10000
    
    fn = os.path.join(dir_, '{0}.hdf5'.format(shot))
    while os.path.exists(fn):
        shot += 1
        fn = os.path.join(dir_, '{0}.hdf5'.format(shot))
        
    return shot


###############################################################################
# Class to represent raw x-ray data. 
###############################################################################
class RawData(object):
    """The RawData object is a container for raw x-ray data for a single 
    detector.

    Members:
    y           -- The data itself.
    shot        -- The shot number the data was taken from.
    detector    -- The detector used to collect the data.
    filter      -- The filter used.
    aperture    -- The aperture used.
    t0          -- The time of the first sample in seconds.
    timestamp   -- The (Unix) timestamp when the data was collected.
    calib       -- The calibration object for the data.
    physPort    -- The physical port.
    shapingTime -- The shaping time.
    gain        -- The gain.
    """
    def __init__(self, detector, shot):
        self.detector   = detector
        self.shot       = int(shot)
        
        self.calib      = None
        
        self.y          = None
        
        self.filter     = None
        self.aperture   = None
        self.t0         = None
        self.timestamp  = None
        
        self._time_axis   = None
        
        # New to takedata3
        self.physPort     = None
        self.shapingTime  = None
        self.gain         = None
        self.filename     = None


    def __hash__(self):
        return int(str(self.shot) + str(self.detector.sn))


    @property
    def time_axis(self):
        """The full time axis."""
        if self._time_axis is None:
            self._time_axis = np.arange(len(self.y)) * self.dt + self.t0
        return self._time_axis


    def save(self, hdf5_file):
        """Save data to the hdf5 file and sqlite3 database."""
        assert(not None in (self.shot, self.y, self.detector, self.dt, 
                            self.t0, self.timestamp))

        # Save raw data to the file.
        name = 'data_{0}'.format(self.detector.sn)
        hdf5_file.create_dataset(name, data=self.y, compression=1)
        
        # CREATE TABLE xRayData (sn INTEGER PRIMARY KEY, filename TEXT, shot INTEGER, detector INTEGER, port TEXT, shapingTime TEXT, gain TEXT, aperture INTEGER, filter INTEGER, t0 REAL, timestamp INTEGER);
        query = 'INSERT INTO xRayData '
        query += '(filename, shot, detector, port, shapingTime, gain,'
        query += ' aperture, filter, t0, timestamp) VALUES '
        query += '(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'

        aperture = None if self.aperture is None else self.aperture.sn
        filter   = None if self.filter is None else self.filter.sn
    
        dbutil.run_commit(TD3_DB_XRAY_PATH, 
                          query, 
                          (name, self.shot, self.detector.sn,
                           self.physPort, self.shapingTime, self.gain, 
                           aperture, filter, self.t0, self.timestamp))


    def load(self, hdf5_file):
        """Load the data from the given hdf5 file and the sqlite3 database."""
        # Load the calibration object associated with this data. 
        self.calib = Calibration(self.detector, _date_from_shot(self.shot))
        
        # Load the raw data. 
        data_name = 'data_{0}'.format(self.detector.sn)
        # Set self.y to a numpy array of the values of the detector with the specified serial number
        self.y = hdf5_file.get(data_name).value

        # Read data from the database. 
        query  = 'SELECT filename, shot, detector, port, shapingTime, '
        query += 'gain, aperture, filter, t0, timestamp, '
        query += 'FROM xRayData WHERE shot=? AND detector=?'
       
        row = dbutil.fetchall(TD3_DB_XRAY_PATH, 
                              query,
                              (self.shot, int(self.detector.sn)))[0]
        
        self.filename    = row[0]
        self.shot        = row[1]
        self.detector    = row[2]
        self.physPort    = row[3]
        self.shapingTime = row[4]
        self.gain        = row[5]
        self.aperture    = row[6]
        self.filter      = row[7]
        self.t0          = row[8]
        self.timestamp   = row[9]
        
        if self.aperture is not None:
            self.aperture = mp.XRayAperture.get(self.aperture)
        if self.filter is not None:
            self.filter = mp.XRayFilter.get(self.filter)
        if self.port is not None:
            self.port = mp.Port.get(self.port)


###############################################################################
# Loading and storing raw data. 
###############################################################################
def _date_from_shot(shot):
    """Given a shot, return the associated date as an integer."""
    if shot / 1e11 > 1:
        return int(shot / 10000)
    else:
        return mds.shot_to_date_num(shot)


def _get_data_path(shot):
    """Get the path to the raw data file for the given shot.
    
    Arguments:
    shot -- The shot for which we'd like the raw data path. Could be an 
            MDSplus shot like 1110216001, or a non-MDSplus shot like
            201102160020, or None to get a new non-MDSplus shot.
    
    Return: 
    dir      -- The directory.
    filename -- The filename. 
    """
    date = _date_from_shot(shot)
    dir_ = os.path.join(const.DATA_DIR, str(date))
    return dir_, '{0}.hdf5'.format(shot)


def save(data_list):
    """Save the list of RawData objects to a file (and the database). The 
    filename to use will depend on whether the shot member is set in the 
    RawData objects.
    
    data_list -- List of RawData objects to save. Each Data object 
                 represents the data from a single detector. 
                
    Return: Nothing.
    """
    # As a sanity check, make sure all shot numbers match. 
    shot = data_list[0].shot
    for data in data_list[1:]:
        assert(data.shot == shot)

    # This will get the correct save path.
    dir_, filename = _get_data_path(shot)
    path = os.path.join(dir_, filename)
    
    # We won't ever overwrite any files. 
    if os.path.exists(path):
        raise Exception('Refusing to overwrite file: {0}'.format(path))
    
    # Make sure the directory exists.
    try:
        os.makedirs(dir_)
    except:
        pass
        
    # Get a list of detector serial numbers. 
    sns = sorted([data.detector.sn for data in data_list])
    
    # Write out the data file.
    hdf5_file = h5py.File(path, 'w')
    hdf5_file.create_dataset('detectors', data=sns)
    
    # When the data saves itself to the file, it will also write 
    # associated data to the sqlite database. 
    for rawdata in data_list:
        rawdata.save(hdf5_file)
        
    # Close the hdf5 file.
    hdf5_file.close()


def load(shot, detector_sn=None):
    """Load data for the given shot.
    
    Arguments:
    shot        -- The shot number. Either an MDSplus shot number, or 
                   a non-MDSplus shot number. 
    detector_sn -- If the serial number is given, then return the data 
                   for the given serial number. Otherwise, return a 
                   list of all the data objects in the file. 
    """
    dir_, filename = _get_data_path(shot)
    path = os.path.join(dir_, filename)
    
    # Create a list of data objects.
    data_list = []
    
    # Open the file for reading. 
    hdf5_file = h5py.File(path, 'r')
    
    # Check for a single serial number. For some reason we need
    # to cast detector_sn and sn below as integers.
    if detector_sn is not None:
        rawdata = RawData(mp.XRayDetector.get(int(detector_sn)), shot)
        rawdata.load(hdf5_file)
        hdf5_file.close()
        return rawdata
    
    # Set sns to a numpy array of the detector serial numbers
    sns = hdf5_file.get('detectors').value
    # Load all data in the file. 
    for sn in sns:
        rawdata = RawData(mp.XRayDetector.get(int(sn)), shot)
        rawdata.load(hdf5_file)
        data_list.append(rawdata)
        
    return data_list
        

        
        
