#!/usr/bin/env python

from distutils.core import setup

scripts = [
    'mstxray-digi-server',
    'mstxray-takedata',
    'mstxray-takedata2',
    'mstxray-takedata-basic',
    'mstxray-procshot',
    'mstxray-procday',
    'mstxray-takedata_CINDE',
    ]

setup(name          = "mstxray",
      version       = "1.0",
      description   = "X-ray processing package.",
      author        = "J. David Lee",
      author_email  = "johnl@cs.wisc.edu",
      maintainer    = "johnl@cs.wisc.edu",
      url           = "none",
      packages      = ["mstxray"],
      package_data  = {'mstxray': ['c_code/*', 'db/*']},
      scripts       = scripts
     )

