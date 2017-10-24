from __future__ import print_function, division
import os
import numpy as np
from scipy.optimize import curve_fit

import cache
from mst import signal


def crop_dead_time(dead_time, tmin, tmax):
    """Crop the dead_time array so that it only contains entries falling in the 
    time window [tmin, tmax].
    """
    if len(dead_time) == 0:
        return dead_time

    mask = np.ones_like(dead_time[:,0], dtype=np.bool)
    mask &= dead_time[:,1] > tmin
    mask &= dead_time[:,0] < tmax

    dead_time = dead_time[mask]

    if len(dead_time) > 0 and dead_time[0][0] < tmin:
        dead_time[0][0] = tmin
        
    if len(dead_time) > 0 and dead_time[-1][1] > tmax:
        dead_time[-1][1] = tmax
        
    return dead_time


def apply_mask(xrdata, mask):
    xrdata.time    = xrdata.time  [mask]
    xrdata.energy  = xrdata.energy[mask]
    xrdata.eff     = xrdata.eff   [mask]
    xrdata.flags   = xrdata.flags [mask]
    

def crop_time(xrdata, tmin, tmax):
    """Crops data arrays in xrdata to the time interval."""
    mask = np.ones(xrdata.time.shape[0], dtype=np.bool)
    mask &= (xrdata.time >= tmin)
    mask &= (xrdata.time < tmax)
        
    xrdata.dead_time = crop_dead_time(xrdata.dead_time, tmin, tmax)
    apply_mask(xrdata, mask)

    xrdata.t_start = tmin
    xrdata.t_end   = tmax


def calc_working_ratio(xrdata):
    time_ratio = 1.0
    if len(xrdata.dead_time) != 0:
        dt = xrdata.t_end - xrdata.t_start
        dt_err = (xrdata.dead_time[:,1] - xrdata.dead_time[:,0]).sum()
        time_ratio = (dt - dt_err) / dt
    
    photon_ratio = 1.0
    if len(xrdata.flags) != 0:
        N = 1.0 * xrdata.flags.shape[0] 
        Nerr = (xrdata.flags != 0).sum()
        photon_ratio = (N - Nerr) / N

    return time_ratio * photon_ratio


def spect(xrdata, Emin, Emax, NE, tmin, tmax, 
          correct_eff=False, mult_energy=False,
          divide_work=False, divide_time=False, 
          divide_energy=False, divide_etendue=False):
    """Generate a spectrum. 

    xrdata -- A mst.signal.XRayData object. Note that xrdata will be modified
              during the spectrum generation. You may use copy.copy in order
              to pass in a copy. 
    Emin -- Emin, Emax and NE are used to construct the energy bins using
    Emax -- np.linspace(Emin, Emax, NE). 
    NE   -- ...
    tmin -- Select photons in a given time interval. Time interval is half-open,
    tmax -- [tmin, tmax).
    
    The following change the weighting of the photons being histogramed. 
    correct_eff -- Correct histogram and error bars for collection efficiency. 
    mult_energy -- Multiply counts by their energy to create an energy spectrum. 

    The following change the scale of the returned result. 
    divide_work    -- Divide by the working ratio (fraction of time 
                      detector is working). 
    divide_time    -- Divide by the total time in the window. 
    divide_energy  -- Divide by width of energy bins. 
    divide_etendue -- Divide by aperature etendue. 
    
    Return: E, H, err
    E   -- Energy bin centers. 
    H   -- Histogram of photons. 
    err -- Error bars. 
    """
    # Energy bin centers. 
    Ebins = np.linspace(Emin, Emax, NE)
    dE = Ebins[1] - Ebins[0]
    E = Ebins[:-1] + (dE / 2.0)

    # First we mask out the time interval. 
    crop_time(xrdata, tmin, tmax)

    # Get working ratio. 
    working_ratio = calc_working_ratio(xrdata)
    
    # Keep good only good photons. 
    apply_mask(xrdata, xrdata.flags == 0)

    # Construct weights. 
    weights = np.ones_like(xrdata.energy)
    
    if correct_eff:
        weights /= xrdata.eff
    if mult_energy:
        weights *= xrdata.energy
        
    # Construct histogram. 
    H = np.histogram(xrdata.energy, Ebins, weights=weights)[0]

    # Calculate error bars.
    err = np.histogram(xrdata.energy, Ebins, weights=weights**2)[0]
    err = np.sqrt(err)
    
    # Scale results. 
    if divide_work:
        H /= working_ratio
        err /= working_ratio
        
    if divide_energy:
        H /= dE
        err /= dE
        
    if divide_time:
        dt = tmax - tmin
        H /= dt
        err /= dt
        
    if divide_etendue:
        H /= xrdata.etendue
        err /= xrdata.etendue
        
    return E, H, err
    

def shot_spect(shot, Emin, Emax, NE, tmin, tmax, **kwargs):
    """Return a list of spectra (signal.XRayData, H, err) for each detector 
    in the given shot.
    
    shot -- The shot number. 
    Emin -- See spect. 
    Emax -- "
    NE   -- "
    tmin -- "
    tmax -- "
    
    Return: E, spect_list
    E       -- Energy bin centers. 
    spectra -- List of spectra. Each element is of the form (signal.XRayData, H, err). 
    """
    E = None
    ret = []
    for xrdata in signal.xray_data(shot):
        try:
            E, H, err = spect(xrdata, Emin, Emax, NE, tmin, tmax, **kwargs)
            ret.append((xrdata, H, err))
        except:
            print('Failed to produce spectrum. Shot: {0} Detector: {1}'.format(
                    shot, xrdata.det_sn))

    return E, ret


def multi_shot_spect(shot_list, Emin, Emax, NE, tmin, tmax, **kwargs):
    """Exactly like shot_spect, but analyses multiple shots."""
    E = None
    ret = []
    for shot in shot_list:
        try:
            E, spectra = shot_spect(shot, Emin, Emax, NE, tmin, tmax, **kwargs)
            ret.extend(spectra)
        except:
            print('Failed to get spectra for shot: {0}'.format(shot))
            
    return E, ret

