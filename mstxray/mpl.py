# coding=UTF-8
"""This file contains functions to assist in plotting mst-related stuff."""
import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import tkSimpleDialog

from mst import mdsplus as mds
from mstxray import physical as mp


def shapely_shell(num_points = 8192):
    """Return a shapely Polygon object representing the shell of MST."""
    from shapely.geometry import Polygon
    
    theta = np.linspace(0, 2 * np.pi, num_points)
    
    # Points defining the outside of the shell.
    ptsOut = np.ndarray(shape=(num_points, 2), dtype=np.double)
    ptsOut[:,0] = mp.a0out * np.cos(theta)
    ptsOut[:,1] = mp.a0out * np.sin(theta)
    
    # Points defining the inside of the shell.
    ptsIn = np.ndarray(shape=(num_points, 2), dtype=np.double)
    ptsIn[:,0] = mp.a0 * np.cos(theta)
    ptsIn[:,1] = mp.a0 * np.sin(theta)
    
    # Return the polygon. 
    return Polygon(ptsOut, (ptsIn,))


def shapely_port_cutter(port):
    """Return a shapely Polygon object that can be subtracted from the shell
    to make a port hole in the shell.
    """
    from shapely.geometry import Polygon
    
    standoffIn = - np.linalg.norm(port.outer_shell_center_pt() - 
                                  port.impact_pt())
    
    pt1, pt2 = port.drill_hole_edge_pts(port.standoff)
    pt4, pt3 = port.drill_hole_edge_pts(standoffIn)

    return Polygon((pt1, pt2, pt3, pt4))


def plot_poloidal_cross_section(ax, geo_axis = True, mag_axis = True, 
                                x_axis = True, y_axis = True, ports = ()):
    """Plot a poloidal cross-section on the given axes. 
    
    Arguments:
    ax       -- The matplotlib axes object to plot into. 
    geo_axis -- Draw the geometric axis. 
    mag_axis -- Draw the magnetic axis. 
    x_axis   -- Draw the x-axis.
    y_axis   -- Draw the y-axis.
    ports    -- A list of ports to draw (port objects).

    Note: This code will draw the cross-section incorrectly if the poloidal 
    angle is zero. 
    """
    # Equal aspect ratio.
    ax.set_aspect('equal')

    from shapely.geometry import MultiPolygon
    polys = MultiPolygon((shapely_shell(), ))

    # Carve out each port. 
    for port in ports:
        cutter = shapely_port_cutter(port)
        polys = polys.difference(cutter)

    if polys.type != 'MultiPolygon':
        polys = MultiPolygon((polys, ))

    # Plot each polygon.
    for poly in polys:
        x, y = poly.exterior.xy
        ax.plot(x, y, 'k-')
    
    if geo_axis:
        ax.plot(0, 0, 'ko')
        
    if mag_axis:
        ax.plot(0.06, 0, 'ko')

    if x_axis:
        ax.plot((-mp.a0, mp.a0), (0, 0), 'k--')
        
    if y_axis:
        ax.plot((0, 0), (-mp.a0, mp.a0), 'k--')



def plot_port(ax, port, outer_shell = True, inner_shell = True, 
              standoff = True, line_of_sight = True, impact_point = True,
              exit_point = True):
    """Plot port geometry.
    
    Arguments:
    ax            -- The matplotlib axis object to plot onto. 
    port          -- The port to draw. 
    outer_shell   -- Draw the point on the outer shell.
    inner_shell   -- Draw the point on the inner shell.
    standoff      -- Draw the standoff. 
    line_of_sight -- Draw the line of sight of the port. 
    impact_point  -- Draw the impact point. 
    exit_point    -- Draw the exit point on the inner shell. 
    """
    if outer_shell:
        x, y = port.outer_shell_center_pt()
        ax.plot(x, y, 'bo')
        
    if inner_shell:
        x, y = port.inner_shell_center_pt()
        ax.plot(x, y, 'bo')
        
    if line_of_sight:
        x1, y1 = port.outer_shell_center_pt()
        if standoff:
            x1, y1 = port.drill_hole_center_pt(port.standoff)
        x2, y2 = port.exit_pt()
        
        ax.plot((x1, x2), (y1, y2), 'b-')
    
    if impact_point:
        x, y = port.impact_pt()
        ax.plot(x, y, 'bo')
        
    if exit_point:
        x, y = port.exit_pt()
        ax.plot(x, y, 'bo')
        
    if standoff:
        from shapely.geometry import LineString, Polygon

        shell = Polygon(shapely_shell().exterior)
        
        pt1, pt3 = port.drill_hole_edge_pts(-(mp.a0out - mp.a0))
        pt2, pt4 = port.drill_hole_edge_pts(port.standoff)
        
        l1 = LineString((pt1, pt2))
        l2 = LineString((pt3, pt4))
        l1 = l1.difference(shell)
        l2 = l2.difference(shell)
       
        x, y = l1.xy
        ax.plot(x, y, 'k-')
        x, y = l2.xy
        ax.plot(x, y, 'k-')


def plot_filter(ax, filter, port, standoff = 0):
    """Plot a filter.
    
    Arguments:
    ax       -- The matplotlib axis object to plot onto. 
    filter   -- The filter to draw. 
    port     -- The port to place the filter on. 
    standoff -- The standoff from the port mounting surface. 
    """
    
    f = filter
    p = port
    standoff += p.standoff

    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(standoff, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(standoff + f.tubeLen, p.radius)
    ((x5, y5), (x6, y6)) = p.drill_hole_edge_pts(standoff + f.tubeLen, f.radius)
    ((x7, y7), (x8, y8)) = p.drill_hole_edge_pts(standoff + f.thk, f.radius)
    
    ax.plot((x1, x2, x4, x6, x8, x7, x5, x3, x1),
            (y1, y2, y4, y6, y8, y7, y5, y3, y1),
            'c-')


def plot_aperture(ax, ap, port, standoff = 0):
    """Plot an aperature.
    ax       -- The matplotlib axis object to plot onto. 
    ap       -- The aperture to draw. 
    port     -- The port to place the filter on. 
    standoff -- The standoff from the port mounting surface. 
    """
    a = ap
    p = port

    standoff += p.standoff

    # Plot the tube walls.
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(standoff, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(standoff + a.tubeLen, p.radius)
    
    ax.plot((x1, x3), (y1, y3), 'k-')
    ax.plot((x2, x4), (y2, y4), 'k-')

    # Plot the bottom aperture plate.
    soBotBot = standoff + a.standoff
    soBotTop = soBotBot + a.plateThk
    
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(soBotBot, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(soBotBot, a.botRadius)
    ((x5, y5), (x6, y6)) = p.drill_hole_edge_pts(soBotTop, a.botRadius)
    ((x7, y7), (x8, y8)) = p.drill_hole_edge_pts(soBotTop, p.radius)
    ax.plot((x1, x3, x5, x7), (y1, y3, y5, y7), 'k-')
    ax.plot((x2, x4, x6, x8), (y2, y4, y6, y8), 'k-')
    
    # Plot the top aperture plate.
    soTopBot = standoff + a.standoff + a.plateThk + a.spLen
    soTopTop = soTopBot + a.plateThk

    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(soTopBot, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(soTopBot, a.topRadius)
    ((x5, y5), (x6, y6)) = p.drill_hole_edge_pts(soTopTop, a.topRadius)
    ((x7, y7), (x8, y8)) = p.drill_hole_edge_pts(soTopTop, p.radius)
    ax.plot((x1, x3, x5, x7), (y1, y3, y5, y7), 'k-')
    ax.plot((x2, x4, x6, x8), (y2, y4, y6, y8), 'k-')



def plot_detector(ax, det, port, standoff = 0, label=True):
    """Plot a detector.
    ax       -- The matplotlib axis object to plot onto. 
    det      -- The detector to draw. 
    port     -- The port to place the detector on. 
    standoff -- The standoff from the port mounting surface.
    """
    d = det
    p = port
    
    standoff += p.standoff
    
    # Draw detector standoff.
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(standoff, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(standoff + d.standoff, 
                                                 p.radius)
    ax.plot((x1, x3), (y1, y3), 'k-')
    ax.plot((x2, x4), (y2, y4), 'k-')
    
    standoff += d.standoff

    # Plot detector body.
    soTop = standoff + 0.08
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(standoff, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(standoff, p.radius * 1.25)
    ((x5, y5), (x6, y6)) = p.drill_hole_edge_pts(soTop, p.radius * 1.25)

    ax.plot((x1, x3, x5, x6, x4, x2),
            (y1, y3, y5, y6, y4, y2), 'k-')

    # Plot the label.
    if label:
        x, y = p.drill_hole_center_pt(standoff + 0.04)
        ax.text(x, y, str(d.sn), ha='center', va='center',
                rotation=p.pTilt * 180 / np.pi - 180)

    # plot the window.
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(standoff, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(standoff + d.winThk, p.radius)
    ax.plot((x1, x2), (y1, y2), 'b-')
    ax.plot((x3, x4), (y3, y4), 'b-')

    # Plot the internal columnator.
    soBot = standoff + d.colOffset
    soTop = soBot + d.colThk
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(soBot, p.radius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(soBot, d.colRadius)
    ((x5, y5), (x6, y6)) = p.drill_hole_edge_pts(soTop, d.colRadius)
    ((x7, y7), (x8, y8)) = p.drill_hole_edge_pts(soTop, p.radius)

    ax.plot((x1, x3, x5, x7, x1), (y1, y3, y5, y7, y1), 'k-')
    ax.plot((x2, x4, x6, x8, x2), (y2, y4, y6, y8, y2), 'k-')

    # Plot the detector.
    detBot = standoff + d.detOffset
    detTop = detBot + d.detThk
    ((x1, y1), (x2, y2)) = p.drill_hole_edge_pts(detBot, d.detRadius)
    ((x3, y3), (x4, y4)) = p.drill_hole_edge_pts(detTop, d.detRadius)
    
    ax.plot((x1, x2, x4, x3, x1), (y1, y2, y4, y3, y1), 'r-')
    

def plot_view_cone(ax, port, filter, ap, det):
    from shapely.geometry import LineString

    p = port
    f = filter
    a = ap
    d = det

    # Calculate the bottom aperature offset.
    soApBot = p.standoff
    if f is not None:
        soApBot += f.tubeLen

    soApBot += a.standoff
    
    # Get two points on bottom of aperature.
    (ptBot1, ptBot2) = p.drill_hole_edge_pts(soApBot, a.botRadius)
    
    # Get two point on the top aperature.
    soApTop = soApBot + a.plateThk * 2 + a.spLen

    (ptTop1, ptTop2) = p.drill_hole_edge_pts(soApTop, a.topRadius)
    
    # Get two points defining the detector plane. 
    soDetBot = p.standoff
    if f is not None:
        soDetBot += f.tubeLen
        
    soDetBot += a.tubeLen
    
    soDetBot += d.standoff + d.detOffset

    (ptDet1, ptDet2) = p.drill_hole_edge_pts(soDetBot, d.detRadius)

    # First line.
    # Intersection with interior shell.
    (pt1, pt2) = lineCircIntersect(ptTop1, ptBot2,  R_POL_IN) 
    # Intersection with detector.
    pt3 = lineLineIntersect(ptTop1, ptBot2, ptDet1, ptDet2)
    
    line = LineString((pt1, pt2, pt3))
    x, y = line.xy
    ax.plot(x, y, 'r-')

    # Second line.
    (pt1, pt2) = lineCircIntersect(ptTop2, ptBot1,  R_POL_IN) 
    # Intersection with detector.
    pt3 = lineLineIntersect(ptTop2, ptBot1, ptDet1, ptDet2)
    
    line = LineString((pt1, pt2, pt3))
    x, y = line.xy
    ax.plot(x, y, 'r-')
    
