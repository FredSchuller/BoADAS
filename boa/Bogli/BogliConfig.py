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
NAM: BogliConfig.py (file)
DES: contains Bogli configuration parameters
"""
__version__ = '$Revision: 2109 $'
__date__    = '$Date: 2007-06-03 14:33:01 +0200 (Sun, 03 Jun 2007) $'

from Numeric import *
import os.path

# lut directory
lutDir =  os.path.dirname(__file__)+'/lut/'

# Description of the Label along the X axis
xAxis={'charheight'     : 1.0, \
       'color'          : 1,   \
       'displacement'   : 3.0, \
       'coordinate'     : 0.5, \
       'justification'  : 0.5, \
       'side'           : 'b', \
       'text'           : ' ', \
       'log'            : 0, \
       'draw_it'        : 1, \
       'limits'         : array([0., 0.])}

# Description of the Label along the Y axis, same definition as
# LabelX
yAxis={'charheight'     : 1.0, \
       'color'          : 1,   \
       'displacement'   : 6.0, \
       'coordinate'     : 0.5, \
       'justification'  : 0.5, \
       'side'           : 'l', \
       'text'           : ' ', \
       'log'            : 0, \
       'draw_it'        : 1, \
       'limits'         : array([0., 0.])}

# Description of the Caption, same definition as LabelX
cLabel = {'charheight'      : 1.0, \
         'color'           : 1,   \
         'displacement'    : 1.0, \
         'coordinate'      : 0.5, \
         'justification'   : 0.5, \
         'side'            : 't', \
         'text'            : ' ', \
         'draw_it'         : 1}

# Description of the wedge, same definition as LabelX
zAxis={'charheight'   : 2.0, \
       'color'        : 1,   \
       'displacement' : 0.5, \
       'side'         : 'ri',\
       'text'         : ' ', \
       'draw_it'      : 1,\
       'width'        : 2.0, \
       'linestyle'    : 1,   \
       'linewidth'    : 1,   \
       'limits'       : array([0.,0.])}


# Description for the Point style
point={'symbol': 1,\
       'size': 4.0 }

# Description for the default line style
line={'color'     : 1, \
      'linestyle' : 1, \
      'linewidth' :1 }

# Description for the contour lines not used for the moment...
contour={'charheight': 0.5, \
         'color'     : 1, \
         'linestyle' : 1, \
         'linewidth' : 1 }

# Description for the boxes
box = { 'color'     : 1,   \
        'charheight': 1.0, \
        'linestyle' : 1,   \
        'linewidth' : 1,\
        'draw_it'   : 1}

# Description for the output text
xyouttext = {'charheight'     : 0.8,  \
             'color'          : 5,    \
             'angle'          : 0., \
             'justification'  : 0.5}

# Define the margin around the global viepoint : viewpoint are
# defined in a device-independent manner, using a coordinate
# system whose coordinates run from 0 to 1 in both x and y.

globalViewPoint = { 'marging_left'  : 0.15, \
                    'marging_right' : 1-0.02, \
                    'marging_bottom': 0.1,  \
                    'marging_top'   : 1-0.05, \
                    'wedge_size'    : 0.07}

# marging for clearer plots...                     
addFracLim=0.05


def pointSize(size=1):
    """
    DES:    change plotting point size
    INP:    (int) size = point size (default = 1)
    USAGE:  BogliConfig.pointSize(3)
            pointSize(3)
    """
    # 20070603FB created
    
    point['size']=size

