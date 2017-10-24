from __future__ import division, print_function

import numpy as np
import bottleneck as bn

import np_inline
#import simplegrid
import sshfn

from mstxray import rawdata
from pulsefit_block import fit_mpoc_mle as _fit


def calc_threshold(y, window=1024):
    """Calculate the pulse-height threshold.

    Arguments:
    y      -- The raw data.
    window -- The window size.

    The raw data is broken into segments with length given by window. The
    median of the standard deviations of the segments is used to determine a
    reasonable threshold.

    Return: The estimated pulse-height threshold.
    """
    # Trim the input array so it can be reshaped.
    excess = y.size % window
    y = y[:y.size - excess]

    # Reshape the array (window columns).
    y = y.reshape((-1, window))

    # Warn if there are too few rows.
    if y.shape[0] < 1000:
        print("Warning:", y.shape[0], "segments in y. More is better.")

    # Take the median standard deviation.
    std = np.median(y.std(1))

    # We choose a threshold that is some multiple of this standard deviation.
    # This seems to give a reasonable value for calibration data from all of
    # our detector types. -jdl 2012-01-20.
    return float(5 * std)


def calc_y_norm(y, threshold, window=1024):
    """Create a "normalized" data array. This is accomplished by subtracting
    a rolling median filter from the original data.

    Arguments:
    y         -- The raw data.
    threshold -- The pulse detection threshold.
    window    -- The window size.

    Return: A raw data array with non-zero offsets removed.
    """
    y = y.astype(np.float64)
    y_med = bn.move_median(y, window)
    y_med[:window - 1] = 0
    y_norm = y - y_med
    y_norm[:window - 1] = 0
    return y_norm


_find_peak_indices_c_code = r"""
npy_int64 i, i_peak, num_peaks = 0;
npy_float64 peak_value;

i = 0;
while(1) {
  // Skip samples below threshold.
  while(i < y_norm_shape(0) && y_norm(i) < threshold) {
    ++i;
  }

  // Exit conditions.
  if(i == y_norm_shape(0)) {
    break;
  }

  // Find the next peak.
  i_peak = i;
  y_norm_assert(i);
  peak_value = y_norm(i);
  ++i;
  y_norm_assert(i);
  while(i < y_norm_shape(0) && y_norm(i) >= threshold) {
    if(y_norm(i) > peak_value) {
      y_norm_assert(i);
      peak_value = y_norm(i);
      i_peak = i;
    }
    ++i;
  }

  // Save the peak.
  peak_inds_assert(num_peaks);
  peak_inds(num_peaks) = i_peak;
  ++num_peaks;
}

return_val = num_peaks;
"""

_find_peak_indices_c = np_inline.inline(
    py_types=((float, 'threshold'),),
    np_types=((np.float64, 1, 'y_norm'),
              (np.int64,   1, 'peak_inds')),
    code=_find_peak_indices_c_code,
    return_types=((int, 'return_val'),))


def find_peak_indices(y_norm, threshold):
    """Calculate the indices of the peaks in the normalized y data.

    Arguments:
    y_norm    -- The normalized y data returned by calc_y_norm.
    threshold -- The pulse detection threshold.

    Return: An array of indices where peaks are located.
    """
    peak_inds = np.zeros(y_norm.shape, dtype=np.int64);
    num_peaks = _find_peak_indices_c(threshold, y_norm, peak_inds)
    peak_inds.resize([num_peaks])
    return peak_inds


_calc_pulse_fwhm_c_code = r"""
npy_int64 i, i_peak, i_left, i_right;
npy_float threshold;

for(i = 0; i < peak_inds_shape(0); ++i) {
  peak_inds_assert(i);
  i_peak = peak_inds(i);

  y_norm_assert(i_peak);
  threshold = y_norm(i_peak) / 2.0;

  i_left = i_peak - 1;
  y_norm_assert(i_left);
  while(i_left > 0 && y_norm(i_left) > threshold) {
    --i_left;
    y_norm_assert(i_left);
  }

  i_right = i_peak + 1;
  y_norm_assert(i_right);
  while(i_right < y_norm_shape(0) && y_norm(i_right) > threshold) {
    ++i_right;
    y_norm_assert(i_right);
  }
  widths_assert(i);
  widths(i) = i_right - i_left;
}
"""

_calc_pulse_fwhm_c = np_inline.inline(
    np_types=((np.float64, 1, 'y_norm'),
              (np.int64, 1, 'peak_inds'),
              (np.int64, 1, 'widths')),
    code=_calc_pulse_fwhm_c_code)


def calc_pulse_fwhm(y_norm, peak_inds, threshold):
    """Calculate the approximate pulse full-width at half-maximum.

    Arguments:
    y_norm    -- The normalized y data returned by calc_y_norm.
    peak_inds -- Indices of the peaks returned by find_peak_indices.
    threshold -- The pulse detection threshold.

    Return: The pulse width in samples.
    """
    widths = np.zeros(peak_inds.shape, dtype=np.int64)
    _calc_pulse_fwhm_c(y_norm, peak_inds, widths)
    return np.median(widths)


def fit_pulse_to_data(y, p):
    """Find the best fit parameters for fitting p to y. This is a
    two parameter fit, A * p + b.

    Arguments:
    y -- The raw data (same shape as p).
    p -- The pulse shape.

    Return: A, b
    A   -- The amplitude.
    b   -- The offset.
    """
    N = p.size
    P = p.sum()
    Q = (p * p).sum()
    D = y.sum()
    F = (y * p).sum()

    # Best fit parameters. "A" is the amplitude, and "b" is the offset.
    b = (D * Q - F * P) / (N * Q - P**2)
    A = (F - P * b) / Q

    return A, b


def find_best_fit_index_shift(y, p, fwhm):
    """Find the shift in index at which p fits y with the least error.

    Arguments:
    y    -- The raw data. Same shape as p.
    p    -- The pulse shape.
    fwhm -- The pulse full-width at half-maximum.

    Return: i_shift, A, b
    i_shift -- The shift in index at which p fits y with the least error.
               This is an integer in the range [-fwhm/4, fwhm/4].
    A       -- The best fit amplitude.
    b       -- The best fit offset.
    """
    i_start = int(p.argmax() - 1.5 * fwhm)
    i_end   = i_start + 2 * fwhm
    p = p[i_start:i_end]

    err_best = None
    shift_best = None
    A_best = None
    b_best = None

    for i_shift in range(int(-fwhm/4), int(fwhm/4)):
        y_test = y[i_start + i_shift:i_end + i_shift]
        A, b = fit_pulse_to_data(y_test, p)
        res = (A * p + b - y_test)
        err = np.sum(res**2)
        if err_best is None or err < err_best:
            err_best = err
            shift_best = i_shift
            A_best = A
            b_best = b

    return shift_best, A_best, b_best


def get_pulse_sum(shot, detector_sn, win_pre, win_post, p=None):
    """Generate a sum of pulses found in the given shot.

    Arguments:
    shotlist     -- List of shots to use for creating the characteristic pulse.
    detector_sn  -- The detector serial number.
    win_pre      -- Number of pulse half-widths to keep before peak.
    win_post     -- Number of pulse half-widths to keep after peak.
    p            -- Previously created pulse shape.

    Return: ps, n_pulses
    ps       -- The sum of pulses.
    n_pulses -- The number of pulses in the sum.
    """
    threshold = None # The pulse-height threshold.
    y         = None # The raw data.
    y_norm    = None # The "normalized" raw data.
    fwhm      = None # The pulse full-width at half-maximum.
    i_dist    = None # fwhm * half_win.
    peak_inds = None # The indices of the peaks.
    pulse_sum = None # The sum of good pulses (for averaging).
    n_pulses  = 0    # The number of pulses.

    data = rawdata.load(shot, detector_sn)
    y = data.y.astype(np.float64)

    # Threshold.
    threshold = calc_threshold(y)

    y_norm    = calc_y_norm(y, threshold)
    peak_inds = find_peak_indices(y_norm, threshold)

    # Pulse FWHM.
    if p is None:
        fwhm = calc_pulse_fwhm(y_norm, peak_inds, threshold)
    else:
        fwhm = int(p.shape[0] / (win_pre + win_post))

    i_pre  = fwhm * win_pre
    i_post = fwhm * win_post

    # Find indices of peaks that have the appropriate distances to their
    # nearest neighbors.
    i_dist = max(i_pre, i_post)
    mask = np.ones(peak_inds.shape, np.bool)
    mask[0]  = peak_inds[0] > i_dist
    mask[-1] = (y.shape[0] - peak_inds[-1]) > i_dist
    mask[1:]  &= (peak_inds[1:] - peak_inds[:-1]) > i_dist
    mask[:-1] &= (peak_inds[1:] - peak_inds[:-1]) > i_dist

    peak_inds = peak_inds[mask]

    # Create characteristic and sum pulse.
    if p is None:
        i = peak_inds[0]
        pulse_sum = y[i - i_pre:i + i_post].copy()
        pulse_sum -= pulse_sum[0]
        pulse_sum /= pulse_sum.max()
        n_pulses = 1
        p = pulse_sum / n_pulses
        p /= p.max()
        peak_inds = peak_inds[1:]
    else:
        pulse_sum = np.zeros_like(p)

    # Loop through peaks.
    for i_peak in peak_inds:
        y_ = y[i_peak - i_pre:i_peak + i_post]

        # Shift the peak.
        i_shift, A, b = find_best_fit_index_shift(y_, p, fwhm)
        i_peak += i_shift
        y_ = y[i_peak - i_pre:i_peak + i_post]

        # If the residual is within the threshold, we can update the
        # characteristic pulse.
        res = np.absolute(A * p + b - y_)
        if res.max() < threshold:
            y_ = y_ - b
            y_ /= y_.max()
            pulse_sum += y_
            pulse_sum -= pulse_sum[0]
            n_pulses += 1
            p = pulse_sum / n_pulses
            p /= p.max()

    # Remove offset.
    return pulse_sum, n_pulses


def get_char_pulse(shotlist, detector_sn, win_pre=4, win_post=20):
    """Generate a characteristic pulse from the given shotlist.

    Arguments:
    shotlist     -- List of shots to use for creating the characteristic pulse.
    detector_sn  -- The detector serial number.
    win_pre      -- Number of pulse half-widths to keep before peak.
    win_post     -- Number of pulse half-widths to keep after peak.

    Return: p, n_pulses
    p        -- A characteristic pulse. The pulse will need to be cropped
                and normalized.
    n_pulses -- The number of pulses averaged to make up the characteristic
                pulse.
    """
    pulse_sum = None # The sum of good pulses (for averaging).
    n_pulses  = None # The number of pulses.
    p         = None # The current best characteristic pulse.

    # Create initial pulse to use for the rest of the data.
    shot = shotlist[0]
    pulse_sum, n_pulses = get_pulse_sum(shot, detector_sn, win_pre, win_post)
    p = pulse_sum / n_pulses
    p -= p[0]
    p /= p.max()
    fwhm = (p >= 0.5).sum()

    shotlist = shotlist[1:]

    client = sshfn.Client()

    for shot in shotlist:
        client.add_job(get_pulse_sum, shot, detector_sn, win_pre, win_post, p)

    results = client.run_jobs()

    for result in results:
        if isinstance(result, Exception):
            raise result

        (pc, nc) = result
        i_shift, A, b = find_best_fit_index_shift(pulse_sum, pc, fwhm)

        if i_shift < 0:
            pulse_sum[:i_shift] += pc[-i_shift:]
        elif i_shift > 0:
            pulse_sum[i_shift:] += pc[:-i_shift]
        else:
            pulse_sum += pc

        n_pulses += nc

    return pulse_sum, n_pulses


def crop_normalize_pulse(p):
    """Crop and normalize the characteristic pulse.

    Arguments:
    p -- The characteristic pulse returned by get_char_pulse.

    Return: A cropped and normalized version of the pulse.
    """
    # Find the approximate start of the pulse.
    i_start = p.argmax()
    while i_start > 0 and p[i_start - 1] < p[i_start]:
        i_start -= 1

    # Crop from the front of the pulse and reset the zero level.
    p = p[i_start:].copy()
    p -= p[0]
    p /= p.max()

    return p


def get_amps(shot, detector_sn, p, th):
    """Get pulse amplitudes in the given shot."""
    data = rawdata.load(shot, detector_sn)
    y = data.y.astype(np.float64)
    inds, amps, offset, flags, dead_time = _fit(y, p, th)
    mask = (flags == 0) & (amps > th)
    return amps[mask]


def get_amp_hist(shotlist, detector_sn, p, th, n_bins=200):
    """Get an amplitude histogram for the given detector.

    Arguments:
    shotlist    -- List of shots to pull data from.
    detector_sn -- Detector's serial number.
    p           -- The characteristic pulse shape.
    th          -- The amplitude threshold.
    n_bins      -- Number of bins to use.

    Return: A histogram of the amplitudes found in the given list of shots,
            and the center amplitude of each bin.
    """
    client = sshfn.Client()

    for shot in shotlist:
        client.add_job(get_amps, shot, detector_sn, p, th)

    amps_list = client.run_jobs()

    amps = np.concatenate(amps_list)

    hist, edges = np.histogram(amps, n_bins)
    E = edges[:-1] + (edges[1:] - edges[:-1]) / 2
    return hist, E
