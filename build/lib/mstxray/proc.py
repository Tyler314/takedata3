
from __future__ import print_function, division

import numpy as np
import MDSplus as mds

import rawdata
import physical as mp

from mst.signal import XRayData  # Struct for processed data.


PROC_VERSION = 201408


def _config_eff(energies, data):
    """Calculate the total x-ray collection efficiency at the detector for the
    given configuration.
    """
    detector = data.detector
    port     = data.port
    filt     = data.filter
    aper     = data.aperture

    # We start with perfect transmission. We also keep track of the
    # distance the photons have traveled in air.
    trans = np.ones_like(energies)
    air_len = 0

    # Check for a port window.
    if port is not None:
        thickness = port.winThk
        trans *= port.winMat.get_transmission(energies, thickness)

    # Check for a filter.
    if filt is not None:
        thickness = filt.thk
        trans *= filt.mat.get_transmission(energies, thickness)
        air_len += filt.tubeLen - thickness

    # Check for aperture tube.
    if aper is not None:
        air_len += aper.tubeLen + aper.standoff

    # There is some standoff in the detector mounting as well (very small).
    air_len += detector.standoff

    # Now the transmission of the air. The serial numbers of materials
    # (or anything else) should never change. This is dry air near sea
    # level. I think it should be close enough.
    air = mp.Material.get(2)
    trans *= air.get_transmission(energies, air_len)

    # Finally, the detector has an intrinsic efficiency. This information
    # comes from the manufacturer.
    trans *= detector.efficiency(energies)

    return trans


def _energy_lims(value, data):
    """Given some efficiency value, 0 < value < 1, return two
    energies. The first is the lowest energy with an efficiency
    higher than value, the second is the highest energy with a
    efficiency higher than value.

    This function is limited to the range [1, 100] in keV.

    Resolution is about 0.1 keV.
    """
    trans_max = 100.0
    energies = np.linspace(1.0, trans_max, 2000)
    trans = _config_eff(energies, data)
    energies = energies[trans > value]

    E_min = energies[0]
    E_max = energies[-1]

    return E_min, E_max


def _geometry(data):
    """Return detector geometry:
    R, Z in m.
    p_tilt in radians.
    """
    # Calculate standoff.
    standoff = data.port.standoff

    if data.filter is not None:
        standoff += data.filter.tubeLen

    if data.aperture is not None:
        standoff += data.aperture.tubeLen

    standoff += data.detector.standoff

    # We start at the outer shell center point.
    R, Z = data.port.outer_shell_center_pt()

    p_tilt = data.port.pTilt # Tilt in radians.
    R += standoff * np.sin(np.pi - p_tilt)
    Z += standoff * np.cos(np.pi - p_tilt)

    # Add major radius to R.
    R += mp.R0

    return float(R), float(Z), float(p_tilt)


def _get_dead_time_amptek_ev1(data):
    """This flags dead time around the large negative spikes that occur
    when the preamplifier resets.
    """
    # When a reset occurs, the digitizer is driven to zero during the spike.
    # The recovery time varies between the two.
    n_pre = None
    n_post = None

    if data.detector.type == 'Amptek1':
        n_pre  = int(0.002e-3 / data.dt)
        n_post = int(0.045e-3 / data.dt)
    elif data.detector.type == 'eV1':
        n_pre  = int(0.015e-3 / data.dt)
        n_post = int(0.010e-3 / data.dt)
    elif data.detector.type == 'Amptek2':
	n_pre = int(0.002e-3 / data.dt)
	n_post = int(0.4e-3 / data.dt)
    else:
	n_pre = int(0 / data.dt)
	n_post = int(0 / data.dt)
	print('Warning...deadtime calculation. Unrecognized detector type.')

    mask = data.y == 0
    inds = np.where((mask[1:]) & (~mask[:-1]))[0] + 1

    for ind in inds:
        i0 = np.max((0, ind - n_pre))
        i1 = np.min((data.y.shape[0], ind + n_post))
        mask[i0:i1] = 1

    # For simplicity, set first and last mask elements to zero.
    mask[0] = mask[-1] = 0

    # Find start and end points.
    i0s = np.where((mask[1:] == 1) & (mask[:-1] == 0))[0]
    i1s = np.where((mask[1:] == 0) & (mask[:-1] == 1))[0]

    # Create dead times array.
    dead_time = np.empty(shape=(i0s.shape[0], 2), dtype=np.int)
    dead_time[:,0] = i0s
    dead_time[:,1] = i1s

    return dead_time


def _fit_data(data):
    """Fit the data and return times, energies, flags, and dead_time."""
    from pulsefit_block.fit import fit_mpoc_mle as pulsefit

    # Set the pulse fitting threshold.
    th = _energy_lims(0.001, data)[0]

    if data.detector.type == 'eV1' and th < 4.5:
        th = 4.5
    elif data.detector.type == 'Amptek1' and th < 2.5:
        th = 2.5
    elif data.detector.type == 'Amptek2' and th < 2.5:
        th = 2.5

    r = data.y.astype(np.float64) * data.calib.eFactor
    s = data.calib.y_pulse.astype(np.float64)

    inds, amps, offsets, flags, dead_time = pulsefit(r, s, th)

    # Extend dead-time for detector models, and sort according to
    # dead-time start-time.
    dead_time_ = np.concatenate((dead_time, _get_dead_time_amptek_ev1(data)))

    # Create unified dead-time list.
    dead_mask = np.zeros(r.shape[0], dtype=np.int)
    for i0, i1 in dead_time_:
        dead_mask[i0:i1] = 1
    dead_mask[0] = dead_mask[-1] = 0

    d_start = np.where(dead_mask[1:] - dead_mask[:-1] ==  1)[0] + 1
    d_end   = np.where(dead_mask[1:] - dead_mask[:-1] == -1)[0] + 1

    dead_time = np.empty((d_start.size, 2), dtype=np.int)
    dead_time[:, 0] = d_start
    dead_time[:, 1] = d_end

    # Remove pulses occuring during dead periods.
    mask = np.ones(inds.shape[0], dtype=np.bool)
    for i0, i1 in dead_time:
        mask &= (inds < i0) | (inds > i1)

    # Remve amps below threshold.
    mask &= amps > th

    inds    = inds   [mask]
    amps    = amps   [mask]
    offsets = offsets[mask]
    flags   = flags  [mask]

    return (data.time_axis[inds].astype(np.float32),
            amps.astype(np.float32),
            flags.astype(np.int8),
            data.time_axis[dead_time].reshape(-1).astype(np.float32))


def _create_proc_data(data):
    pd = XRayData()

    pd.det_sn = int(data.detector.sn)
    pd.det_type = str(data.detector.type)
    pd.sigma_E = float(data.calib.sigma_E)

    E_min, E_max = _energy_lims(0.1, data)
    pd.E_min = float(E_min)
    pd.E_max = float(E_max)

    pd.etendue = -1 if data.aperture is None else float(data.aperture.etendue)
    pd.port_sn = int(data.port.sn)

    pd.R, pd.Z, pd.p_tilt = _geometry(data)

    pd.t_start = float(data.time_axis[0])
    pd.t_end = float(data.time_axis[-1] + data.dt)

    pd.time, pd.energy, pd.flags, pd.dead_time = _fit_data(data)

    pd.eff = _config_eff(pd.energy, data)

    pd.efficiency_E = np.linspace(1, 200, 2001)
    pd.efficiency = _config_eff(pd.efficiency_E, data)

    return pd


def create_mds_tree(shot, pd_list):
    """Create the MDSplus tree for the given shot.

    shot         : The shot number.
    detector_sns : The detector serial numbers.
    """
    print('#' * 79)
    print('Creating tree for shot:', shot)
    tree = mds.Tree('MST_HXR', shot, 'NEW')
    tree.write()

    print('Opening tree for editing...')
    tree = mds.Tree('MST_HXR', shot, 'EDIT')

    print('proc_version'); tree.addNode(':PROC_VERSION', 'NUMERIC')
    print('det_in_use'  ); tree.addNode(':DET_IN_USE', 'NUMERIC')

    for pd in pd_list:
        print('#' * 59)
        print('Detector:', pd.det_sn, '--', shot)

        node = tree.addNode(':DETECTOR_{0}'.format(pd.det_sn), 'STRUCTURE')
        old_default = tree.setDefault(node)

        print('Adding nodes...')
        print('    det_type'    ); tree.addNode(':DET_TYPE',     'TEXT')
        print('    sigma_E'     ); tree.addNode(':SIGMA_E',      'NUMERIC')
        print('    E_min'       ); tree.addNode(':E_MIN',        'NUMERIC')
        print('    E_max'       ); tree.addNode(':E_MAX',        'NUMERIC')
        print('    etendue'     ); tree.addNode(':ETENDUE',      'NUMERIC')
        print('    port_sn'     ); tree.addNode(':PORT_SN',      'NUMERIC')
        print('    R'           ); tree.addNode(':R',            'NUMERIC')
        print('    Z'           ); tree.addNode(':Z',            'NUMERIC')
        print('    p_tilt'      ); tree.addNode(':P_TILT',       'NUMERIC')
        print('    t_start'     ); tree.addNode(':T_START',      'NUMERIC')
        print('    t_end'       ); tree.addNode(':T_END',        'NUMERIC')
        print('    efficiency_E'); tree.addNode(':EFFICIENCY_E', 'NUMERIC')
        print('    efficiency'  ); tree.addNode(':EFFICIENCY',   'NUMERIC')
        print('    time'        ); tree.addNode(':TIME',         'NUMERIC')
        print('    energy'      ); tree.addNode(':ENERGY',       'NUMERIC')
        print('    eff'         ); tree.addNode(':EFF',          'NUMERIC')
        print('    flags'       ); tree.addNode(':FLAGS',        'NUMERIC')
        print('    dead_time'   ); tree.addNode(':DEAD_TIME',    'NUMERIC')

        tree.setDefault(old_default)

    print('Writing tree...')
    tree.write()


def save_mds_data(shot, pd_list):
    """Save the list of processed data to the HXR tree."""
    print('#' * 79)
    print('Saving data for shot:', shot)
    print('Opening tree for editing...')
    tree = mds.Tree('MST_HXR', shot, 'EDIT')

    print('Processing code version:', PROC_VERSION)
    tree.getNode(':PROC_VERSION').putData(PROC_VERSION)

    det_in_use = np.asarray([pd.det_sn for pd in pd_list], dtype=np.int32)
    print('Detectors in use:')
    print(det_in_use)
    tree.getNode(':DET_IN_USE').putData(det_in_use)

    for pd in pd_list:
        print('#' * 59)
        print('Detector:', pd.det_sn, '--', shot)

        node = tree.getNode(':DETECTOR_{0}'.format(pd.det_sn))
        old_default = tree.setDefault(node)

        print('det_type', pd.det_type);
        tree.getNode(':DET_TYPE').putData(pd.det_type)
        print('sigma_E', pd.sigma_E)
        tree.getNode(':SIGMA_E').putData(pd.sigma_E)
        print('E_min', pd.E_min)
        tree.getNode(':E_MIN').putData(pd.E_min)
        print('E_max', pd.E_max)
        tree.getNode(':E_MAX').putData(pd.E_max)
        print('etendue', pd.etendue)
        tree.getNode(':ETENDUE').putData(pd.etendue)
        print('port_sn', pd.port_sn)
        tree.getNode(':PORT_SN').putData(pd.port_sn)
        print('R', pd.R)
        tree.getNode(':R').putData(pd.R)
        print('Z', pd.Z)
        tree.getNode(':Z').putData(pd.Z)
        print('p_tilt', pd.p_tilt)
        tree.getNode(':P_TILT').putData(pd.p_tilt)
        print('t_start', pd.t_start)
        tree.getNode(':T_START').putData(pd.t_start)
        print('t_end', pd.t_end)
        tree.getNode(':T_END').putData(pd.t_end)

        print('Writing efficiency_E...')
        tree.getNode(':EFFICIENCY_E').putData(pd.efficiency_E)
        print('Writing efficiency...')
        tree.getNode(':EFFICIENCY').putData(pd.efficiency)
        tree.write()

        print('Writing time...')
        tree.getNode(':TIME').putData(pd.time); tree.write()
        print('Writing energy...')
        tree.getNode(':ENERGY').putData(pd.energy); tree.write()
        print('Writing eff...')
        tree.getNode(':EFF').putData(pd.eff); tree.write()
        print('Writing flags...')
        tree.getNode(':FLAGS').putData(pd.flags); tree.write()
        print('Writing dead_time...')
        tree.getNode(':DEAD_TIME').putData(pd.dead_time); tree.write()

        tree.setDefault(old_default)

    print('Writing tree...')
    tree.write()
    print('Done.')


def mds_save_shot(shot, pd_list):
    """Save the given processed-data list to MDSplus for the given shot."""
    create_mds_tree(shot, pd_list)
    save_mds_data(shot, pd_list)


def proc_shot(shot):
    try:
        data_list = rawdata.load(shot)
    except Exception, ex:
        print('Failed to load raw data for shot:', shot)
        print('Exception', ex)
        return None

    pd_list = []
    for data in data_list:
        try:
            pd_list.append(_create_proc_data(data))
        except Exception, ex:
            print("Failed to process data for sn:", data.detector.sn)
            print('Exception:', ex)

    return pd_list
