# Copyright (C) 2002-2006
# Max-Planck-Institut fuer Radioastronomie Bonn
# Argelander Institut fuer Astronomie
# Astronomisches Institut der Ruhr-Universitaet Bochum
#
# Produced for the LABOCA project
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Library General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
# details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#

"""
DES: module used to plot different symbols and forms
"""

__version__ = '$Revision: 2314 $'
__data__    = '$Date: 2007-08-08 10:43:50 +0200 (Wed, 08 Aug 2007) $'

import Plot
from ppgplot import  *
from Numeric import *


def ellipse(x,y,xFWHM, yFWHM, tilt,overplot=1, style='l', ci=1):
    """
    DES: draw a ellipse
    INP: (f) x/y      : the position offset
         (f) x/y FWHM : the minor and major axis
         (f) tilt     : the angle of the ellipse (radian)
         (l) overplot : should we overplot the ellipse (default: yes)
         (s) style    : the style of the plot (default: line)
         (i) ci       : the color index (default white)
    """

    theta = arange(31.)/30*2*pi

    xp = xFWHM/2.*cos(theta)
    yp = yFWHM/2.*sin(theta)

    ix = xp*cos(-1.*tilt)-yp*sin(-1.*tilt)
    iy = xp*sin(-1.*tilt)+yp*cos(-1.*tilt)

    ix = ix+x
    iy = iy+y

    Plot.plot(ix,iy,overplot=overplot,style=style,ci=ci)

def circle(x, y, radius, fill=0, ci=1):
    """
    DES: draw a (filled) circle using pgplot primitive
    INP: (f) x/y      : the position offset
         (f) radius   : the radisu of the circle
         (s) fill     : should we fill the circle  (default: no)
         (i) ci       : the color index (default white)
    """
    savedFillStyle = pgqfs()
    
    if fill:
        pgsfs(1)
    else:
        pgsfs(2)
        
    pgsci(ci)
    pgcirc(x,y,radius)

    pgsfs(savedFillStyle)
    
def shadeY(Ymin,Ymax,style=3,ci=1):
    """
    DES: shape an area between two Y values
    INP: (f) Y min/max : the limits (could be off plot limits)
         (i) style     : the shading style (see pgplot help)
         (i)    ci     : the color of the shade
    """
    savedFillStyle = pgqfs()

    pgsfs(style)
    pgsci(ci)
    xAxis = Plot.xAxis['limits']
    pgpoly(array([xAxis[0],xAxis[0],xAxis[1],xAxis[1]]),\
           array([Ymin,Ymax,Ymax,Ymin]))
    pgsfs(savedFillStyle)

def shadeX(Xmin,Xmax,style=3,ci=1):
    """
    DES: shape an area between two X values
    INP: (f) X min/max : the limits (could be off plot limits)
         (i) style     : the shading style (see pgplot help)
         (i)    ci     : the color of the shade
    """
    savedFillStyle = pgqfs()

    pgsfs(style)
    pgsci(ci)
    yAxis = Plot.yAxis['limits']
    pgpoly(array([Xmin,Xmax,Xmax,Xmin]),\
           array([yAxis[0],yAxis[0],yAxis[1],yAxis[1]]))
    pgsfs(savedFillStyle)
