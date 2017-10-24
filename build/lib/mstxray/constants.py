"""Some constants used throughout the x-ray code."""

import os

_this_dir = os.path.split(__file__)[0]

XRAY_PATH  = os.path.expanduser('~/x-ray')
DATA_DIR   = os.path.join(XRAY_PATH, 'xRayData')
PULSE_DIR  = os.path.join(XRAY_PATH, 'xRayCharPulse')

DB_PHYSICAL_PATH = os.path.join(_this_dir, 'db', 'mstPhysical.sqlite')
DB_XRAY_PATH     = os.path.join(XRAY_PATH, 'xRay.sqlite')

ADDR_XRAY_SRV    = 'dave.physics.wisc.edu'
PORT_DIGITIZER   = 18495
PORT_DB_PHYSICAL = 18496
PORT_DB_XRAY     = 18497
PORT_FILE_SRV    = 18498

