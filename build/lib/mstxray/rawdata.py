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
    y         -- The data itself.
    shot      -- The shot number the data was taken from.
    detector  -- The detector used to collect the data.
    port      -- The port the detector was installed on.
    filter    -- The filter used.
    aperture  -- The aperture used.
    dt        -- The time step in seconds.
    t0        -- The time of the first sample in seconds.
    timestamp -- The (Unix) timestamp when the data was collected.
    calib     -- The calibration object for the data.
    """
    def __init__(self, detector, shot):
        self.detector   = detector
        self.shot       = int(shot)
        
        self.calib      = None
        
        self.y          = None
        
        self.port       = None
        self.filter     = None
        self.aperture   = None
        self.dt         = None
        self.t0         = None
        self.timestamp  = None
        
        self._time_axis = None


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
        
        query = 'INSERT INTO xRayData '
        query += '(shot, detector, aperture, filter, port, dt,'
        query += ' t0, timestamp) VALUES '
        query += '(?, ?, ?, ?, ?, ?, ?, ?)'

        aperture = None if self.aperture is None else self.aperture.sn
        filter = None if self.filter is None else self.filter.sn
        port = None if self.port is None else self.port.sn
        
        dbutil.run_commit(const.DB_XRAY_PATH, 
                          query, 
                          (self.shot, self.detector.sn,
                           aperture, filter, port, self.dt, self.t0,
                           self.timestamp))


    def load(self, hdf5_file):
        """Load the data from the given hdf5 file and the sqlite3 database."""
        # Load the calibration object associated with this data. 
        self.calib = Calibration(self.detector, _date_from_shot(self.shot))
        
        # Load the raw data. 
        data_name = 'data_{0}'.format(self.detector.sn)
        self.y = hdf5_file.get(data_name).value

        # Read data from the database. 
        query = 'SELECT shot, aperture, filter, port, dt, t0, timestamp '
        query += 'FROM xRayData WHERE shot=? AND detector=?'
        
        row = dbutil.fetchall(const.DB_XRAY_PATH, 
                              query,
                              (self.shot, int(self.detector.sn)))[0]
        
        self.shot      = row[0]
        self.aperture  = row[1]
        self.filter    = row[2]
        self.port      = row[3]
        self.dt        = row[4]
        self.t0        = row[5]
        self.timestamp = row[6]
        
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
    sns = hdf5_file.get('detectors').value
    
    # Check for a single serial number.  For some reasone we need
    # to cast detector_sn and sn below as integers.
    if detector_sn is not None:
        rawdata = RawData(mp.XRayDetector.get(int(detector_sn)), shot)
        rawdata.load(hdf5_file)
        hdf5_file.close()
        return rawdata
    
    # Load all data in the file. 
    for sn in sns:
        rawdata = RawData(mp.XRayDetector.get(int(sn)), shot)
        rawdata.load(hdf5_file)
        data_list.append(rawdata)
        
    return data_list
        

        
        
