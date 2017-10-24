# coding=UTF-8
"""Physical attributes of equipment relating to the mst x-ray diagnostic.

Constants defined here:
    R0    -- The major radius.
    a0    -- The minor radius.
    a0out -- The minor radius of the outside of the torus.
"""
###############################################################################
# Imports
###############################################################################
import numpy as np
from scipy.interpolate import interp1d

import dbutil
import constants as const


###############################################################################
# Constants.
###############################################################################
R0    = 1.5  # Major radius.
a0    = 0.52 # Minor radius.
a0out = 0.57 # Minor radius of the outside of the torus.


###############################################################################
# Base database object. 
###############################################################################
class DBObject(object):
    """Not for public use. This is a base class for objects stored in the 
    database allowing them to be used as row factories (see sqlite3 docs).
    """
    __OBJECT_CACHE = {}

    def __init__(self, cursor, row):
        """For use as a row_factory (sqlite3)."""
        
        # Every object has a serial number that uniquely identifies it.
        # Every object also has a name.
        self.sn = None
        self.name = None

        # Other collumns are added from the cursor. 
        cols = [d[0] for d in cursor.description]
        for i, col in enumerate(cols):
            setattr(self, col, row[i])


    def __str__(self):
        return '{0} (sn: {1})'.format(self.name, self.sn)


    def __repr__(self):
        return '<{0} (sn: {1})>'.format(self.name, self.sn)


    @classmethod
    def get(cls, sn):
        """Get an item with the given serial number. If you pass an object of 
        type cls, it will be returned.
        """
        # If the user passed in an object of this class, return it. 
        if isinstance(sn, cls):
            return sn

        try:
            # We first attempt to return the object from cache. 
            return cls.__OBJECT_CACHE[cls][sn]

        except:
            # The object wasn't in cache, so we read it from the database.
            table_name = cls.__name__ + 's'
            query = 'SELECT * FROM {0} WHERE sn=?'.format(table_name)
            obj = dbutil.fetchall(const.DB_PHYSICAL_PATH, query, (sn,), cls)[0]

            # And of course we cache it. 
            if cls not in cls.__OBJECT_CACHE:
                cls.__OBJECT_CACHE[cls] = {}
            cls.__OBJECT_CACHE[cls][sn] = obj

            return obj
    

    @classmethod
    def get_all(cls, in_use=None):
        """Get all items. If in_use is not None, then we only want those items
        that are in use. This condition isn't defined for all objects. 
        """
        table_name = cls.__name__ + 's'
        query = 'SELECT * FROM {0}'.format(table_name)
        query_args = []
        if in_use is not None:
            query += ' WHERE inUse=?'
            query_args.append(in_use)
        query += ' ORDER BY sn'
        objs = dbutil.fetchall(const.DB_PHYSICAL_PATH, query, query_args, cls)

        if cls not in cls.__OBJECT_CACHE:
            cls.__OBJECT_CACHE[cls] = {}

        for obj in objs:
            cls.__OBJECT_CACHE[cls][obj.sn] = obj
            
        return objs
        


###############################################################################
# Public objects.
###############################################################################
class Material(DBObject):
    """Materials are defined in the database. 

    Members:
        sn     -- The material's unique serial number. 
        symbol -- The symbol for the material. 
        name   -- The full name of the material. 
        rho    -- The material's density in kg/m^3.
    """
    def __init__(self, cursor, row):
        DBObject.__init__(self, cursor, row)

        # _mu_o_rho will be used to calculate x-ray transmission efficiency.
        self._mu_o_rho = None
        self._energies = None
        self._load_mu_o_rho()
        

    def _load_mu_o_rho(self):
        """Load the mass attenuation coefficients."""
        query = 'SELECT energy, mu_o_rho FROM xRayMassAttenuationCoeffs'
        query += ' WHERE material=? ORDER BY energy ASC'
        rows = np.array(dbutil.fetchall(const.DB_PHYSICAL_PATH, 
                                        query, (self.sn,)))
        self._energies = rows[:, 0]
        self._mu_o_rho = rows[:, 1]

        
    def get_transmission(self, energies, thickness):
        """Return a the x-ray transmission efficiency at the given energies 
        for the given thickness. (See Beer-Lambert law.)
        
        Arguments:
        energies  -- The energies in keV.
        thickness -- The thickness of the material in m. 
        
        Return: The x-ray transmission efficiency at each energy in 
                energies. 
        """
        trans = np.exp(-self._mu_o_rho * thickness * self.rho)
        func = interp1d(self._energies, trans, bounds_error=False, 
                        fill_value=0, kind=1)
        return func(energies)


class Port(DBObject):
    """Physical port specifications.
    
    Members:
        sn           -- The port's serial number. Must be unique and fixed. 
        name         -- Human readable name. 
        radius       -- The port radius in m. 
        radiusUnc    -- The uncertainty in the radius. 
        winMat       -- Port window material. See material.py.
        winThk       -- Port window thickness in m.
        winThkUnc    -- Port window thickness uncertainty in m.
        tAng         -- Toroidal angle in radians.
        tAngUnc      -- Toroidal angle uncertainty in radians. 
        pAng         -- Poloidal angle in radians.
        pAngUnc      -- Poloidal angle uncertainty in radians.
        pTilt        -- Poloidal tilt angle in radians.
        pTiltUnc     -- Poloidal tilt angle uncertainty in radians.
        standoff     -- Distance from outside wall to mounting surface on
                        center in m. 
        standoffUnc  -- Standoff height uncertainty in m.
        inUse        -- True if port is in use.    
    """
    def __init__(self, cursor, row):
        DBObject.__init__(self, cursor, row)
        self.winMat = Material.get(self.winMat)


    def _x_y(self, r, theta):
        """Given r and theta, return x and y."""
        return r * np.array((np.cos(theta), np.sin(theta)))


    def impact_pt(self):
        """Return the (x,y) coordinates of the impact point. That is, the 
        point on the port center-line closest to the geometric center (0,0) 
        of the poloidal cross-section.
        """
        pTilt = self.pTilt
        pAng  = self.pAng
        r_i = a0out * np.sin(pTilt - pAng - np.pi / 2)
        theta_i = pTilt - np.pi
        return self._x_y(r_i, theta_i)


    def impact_param(self):
        """Return the impact parameter of the port in meters.
        This is a signed quantity -- it is negative if the impact parameter
        is inboard, and positive if it is outboard.
        """
        x, y = self.impact_pt()
        return np.sqrt(x**2 + y**2) * x / np.absolute(x)

    
    def outer_shell_center_pt(self):
        """Return the (x, y) coordinates of the outer shell center point"""
        return self._x_y(a0out, self.pAng)
    
    
    def drill_hole_center_pt(self, standoff):
        """Return the (x, y) coordinates of the a hole's center point at a 
        given standoff from the outer shell. Positive is outward. 
        """
        r_out = self.outer_shell_center_pt()
        r = r_out - self.impact_pt()
        rHat = r / np.linalg.norm(r)
        return r_out + rHat * standoff

    
    def drill_hole_edge_pts(self, standoff, radius = None):
        """Return a set of points, (x1, y1), (x2, y2) that are the edges of a
        drill hole at the standoff through a perpendicular plane. 
        """
        if radius is None:
            radius = self.radius

        # Center point.
        rc = self.drill_hole_center_pt(standoff)
        r_imp = self.impact_pt()
        r_imp_hat = r_imp / np.linalg.norm(r_imp)
        return (rc - radius * r_imp_hat, rc + radius * r_imp_hat)
    
    
    def inner_shell_center_pt(self):
        """Return the (x, y) coordinates of the center point of the port
        on the inner shell.
        """
        norm = np.linalg.norm

        r_out = self.outer_shell_center_pt() 
        r_imp = self.impact_pt()

        r1   = r_imp - r_out
        r1_hat = r1 / norm(r1)

        x = np.sqrt(a0**2 - norm(r_imp)**2)
        r2_mag = norm(r1) - x
        r2 = r2_mag * r1_hat 
        
        return r_out + r2


    def exit_pt(self):
        """Return the point on the outer shell where the port center-line
        would exit the vessel.
        """
        # Construct a vector that crosses the far outer wall.
        r_in  = self.inner_shell_center_pt()
        r_imp = self.impact_pt()
        
        return 2 * r_imp - r_in


class XRayDetector(DBObject):
    """Physical x-ray detector specifications. 
    
    Members:
        sn          -- The detector serial number. Must be unique and 
                       fixed.
        name        -- The human readable name.
        type        -- The detector type. 
        standoff    -- Distance between front of detector and mounting 
                       tube.
        standoffUnc -- Uncertainty in standoff
        winMat      -- Window material on front of detector.
        winThk      -- Window thickness.
        winThkUnc   -- Window thickness uncertainty.
        colOffset   -- Offset from front of detector to front of 
                       columnator.
        colThk      -- Thickness of columnator
        colRadius   -- Radius of the columnator's aperature.
        detOffset   -- Offset from front of detector to chip front.
        detRadius   -- Radius of the detector chip.
        detThk      -- Thickness of detector chip. 
        digiCard    -- Address of digitizer card (VME).
        digiChan    -- Digitizer card channel.
        inUse       -- True if the detector is in use. 

        efficiency  -- Function taking photon energies and returning the 
                       detector efficiency at those energies.
    """
    def __init__(self, cursor, row):
        DBObject.__init__(self, cursor, row)
        self.winMat = Material.get(self.winMat)
        self._efficiency = None
        
        self._energies     = None
        self._efficiencies = None

        self._load_efficiency()

        
        
    def efficiency(self, energies):
        if self._energies is not None:
            fn = interp1d(self._energies, self._efficiencies, 
                          bounds_error=False, fill_value=0, kind=1)
            return fn(energies)
        else:
            return 1.0


    def _load_efficiency(self):
        try:
            query = 'SELECT energy, efficiency FROM '
            query += 'XRayDetectorEfficiencies '
            query += 'WHERE detector=? ORDER BY energy ASC'

            rows = np.array(dbutil.fetchall(
                    const.DB_PHYSICAL_PATH, query, (self.sn,)))
            self._energies = rows[:, 0]
            self._efficiencies = rows[:, 1]

        except Exception:
            print('No detector efficiency for SN {0}'.format(self.sn))
            self._energies     = None
            self._efficiencies = None


class XRayAperture(DBObject):
    """Physical x-ray aperture specifications.
    
    Members:
        sn           -- The filter's unique serial number. Unique and 
                        fixed. 
        name         -- The filter's human readable name.
        topRadius    -- Top plate aperture radius.
        topRadiusUnc -- Top plate aperture radius uncertainty
        botRadius    -- Bottom plate aperture radius.
        botRadiusUnc -- Bottom plate aperture radius uncertainty
        plateThk     -- Plate thickness.
        plateThkUnc  -- Plate thickness uncertainty.
        spLen        -- Spacer length.
        spLenUnc     -- Specer length uncertainty.
        standoff     -- Distance between tube bottom and bottom aperture
                        plate.
        standoffUnc  -- Uncertainty in distance.
        tubeLen      -- Length of the containing tube.
        tubeLenUnc   -- Uncertainty in the tube length. 
        inUse        -- True if the aperature set is in use. 
    """

    @property
    def solidAngle(self):
        """The solid angle seen by the aperture in sr."""
        length = self.spLen + 2 * self.plateThk
        r_bot = self.botRadius
        cos_theta = length / np.sqrt(length**2 + r_bot**2)
        return 2 * np.pi * (1 - cos_theta)
    
    @property
    def etendue(self):
        """Etendue in m^2 * sr."""
        return self.solidAngle * np.pi * self.topRadius**2


class XRayFilter(DBObject):
    """Physical x-ray filter specifications.

    Members:
        sn         -- The filter's unique serial number. Unique and 
                      fixed. 
        name       -- The filter's human readable name.
        radius     -- The radius of the filter. 
        radiusUnc  -- The uncertainty in the radius. 
        mat        -- The filter material's sn (see mst.material).
        thk        -- The thickness of the filter. 
        thkUnc     -- The uncertainty in the filter's thickness. 
        tubeLen    -- The length of the filter's tube.
        tubeLenUnc -- The uncertainty in the filter's tube length. 
        inUse      -- True if the filter is in use. 
    """
    def __init__(self, cursor, row):
        DBObject.__init__(self, cursor, row)
        self.mat = Material.get(self.mat)

