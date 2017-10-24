
from mst import mdsplus as mds
from rawdata import RawData
import physical as mp

def get_photons(shot, chan):
    """Get photons from mdsplus.
    shot -- The shot number.
    chan -- The detector channel or detector number. I'm not sure of how 
            the old system indexed things. 
    Return: t, y
    t -- The photon times in seconds.
    y -- The photon energies in keV. 
    """
    chan_str = '{0:02}'.format(chan)
    sig = r'\mst_hxr::top.hxr_array.fitdata.hxr' + chan_str
    t, y = mds.get_signal(shot, sig)
    return t, y


class OldData(RawData):
    def __init__(self, chan, shot):
        self.data_channel = chan
        RawData.__init__(self, mp.XRayDetector.get(2000), shot)
        self._time_axis = get_photons(shot, chan)[0]
        
    

