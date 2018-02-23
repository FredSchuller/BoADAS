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
NAM: BoaConfig (module)
DES: Set the path/filename for i/o operation
"""

import os

inDir   = '../fits'                 # input directory
outDir  = './'                      # output directory

projID  = 'T-79.F-0002-2007'        # Your ESO project Id
# Date of obs - in file name since 2007/8/13 (optional)
date = '2007-08-14'

inFile   = ''                       # input file 
outFile  = ''                       # output file
rcpPath  = os.getenv('BOA_HOME_RCP') # Path for the rcp file when needed
if rcpPath is None:
    rcpPath = ""
    
online            = 0               # Is boa used online or not
maxMessHandWeight = 4

goodFeedList = ['AC','DC']
DEBUG = 2
