from __future__ import print_function, division

import os
import h5py

import constants as const
import dbutil

###############################################################################
# Database code.
###############################################################################
class Calibration(object):
    """Calibration for a single detecotor.
    
    Members:
    date     -- The first date for which calibration is valid. Int: 20101224.
    detector -- The x-ray detector.
    eFactor  -- The conversion factor from raw data to energy in keV.
    m_sigma  -- The standard deviation an energy measurement varies linearly
    sigma_E  -- with energy as: sigma = E * m_sigma + sigma_E.
    y_pulse  -- The characteristic pulse. 
    filename -- The pulse filename. 
    """
    def __init__(self, detector, data_date):
        """Constructor."""
        self.detector  = detector
        self.data_date = data_date

        self.date     = None
        self.eFactor  = None
        self.m_sigma  = None
        self.sigma_E  = None
        self.y_pulse  = None
        self.filename = None 
        
        self._load()


    def _pulse_path(self):
        """Return the absolute path to the characteristic pulse hdf5 file
        given the filename. If the filename is None, return a path 
        corresponding to the date and detector in the Calibration object.
        
        Return: A path to the characteristic pulse. 
        """
        if self.filename is None:
            self.filename = '{0}_{1}.hdf5'.format(self.date, self.detector.sn)

        return os.path.join(const.PULSE_DIR, str(self.date)[:4], self.filename)
        
        
    def _load(self):
        # Load properties from the database. 
        query = 'SELECT date, eFactor, m_sigma, sigma_E, pulseFilename '
        query += 'FROM xRayCalib '
        query += 'WHERE detector = ? AND date <= ? '
        query += 'ORDER BY date DESC'
        
        try:
            row = dbutil.fetchall(const.DB_XRAY_PATH, query, 
                                  (self.detector.sn, self.data_date))[0]
        except:
            row = None
        
        if row is None:
            print('Warning: No calibration in database.',
                  '\n  Detector:', self.detector.sn,
                  '\n  Date:', self.data_date)
        else:
            (self.date, 
             self.eFactor, 
             self.m_sigma, 
             self.sigma_E, 
             self.filename) = row

        # Load the pulse shape from the HDF5 file. 
        path = self._pulse_path()
        if not os.path.exists(path):

            print('Warning: Characteristic pulse file not found: ', path)
            return 
        
        try:
            f = h5py.File(path, 'r')
            self.y_pulse = f.get('y_pulse').value
            f.close()
        except Exception, ex:
            msg = 'Failed to read pulse file.'
            msg += '\n    Pulse path: {0}'.format(path)
            msg += '\n    Exception: {0}'.format(ex)
            raise Exception(msg)
        

    def save(self):
        """Save the calibration data to the disk and database. All members
        must be set. 
        """
        # Check for any errors coming in here.
        assert(None not in (self.date, self.detector, self.eFactor, 
                            self.m_sigma, self.sigma_E, self.y_pulse))
        
        # The first order of business is the save the characteristic pulse.
        path = self._pulse_path()

        # Make sure pulse path exists.
        try:
            os.makedirs(os.path.split(path)[0])
        except:
            pass
        
        # Don't overwrite the pulse. 
        if os.path.exists(path):
            raise Exception("Refusing to overwrite pulse: {0}".format(path))

        h5py_file = h5py.File(path, 'w')
        h5py_file.create_dataset('y_pulse', data=self.y_pulse)
        h5py_file.close()
        
        # Next, we need to store the associated information in the database.
        query  = 'INSERT OR REPLACE INTO xRayCalib'
        query += '(date, detector, eFactor, m_sigma, sigma_E, pulseFilename) '
        query += 'VALUES(?, ?, ?, ?, ?, ?)'

        dbutil.run_commit(const.DB_XRAY_PATH, 
                          query, (self.date, self.detector.sn, self.eFactor, 
                                  self.m_sigma, self.sigma_E, self.filename))
        



