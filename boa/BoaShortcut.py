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
####################################################################
# update started 08/12/06 by C.Vlahakis
# update continued 09/03/07 by C. Vlahakis
# more updating and re-organisation on 20.03.07 by C. Vlahakis
####################################################################

"""
DES: Little module to define Boa Shortcut, updated
USE: execfile('BoaShortcut.py') to reset the current shortcut
"""

#from boa import BoaReduce
from boa import BoaConfig, BoaDir
from boa import BoaMessageHandler
from boa import BoaMapping ,BoaPointing, BoaFocus
from boa.Bogli import *

# Baseline, noise and statistics related
####################################################################
polynomialBaseline = baseline    = base       = data.polynomialBaseline
medianNoiseRemoval = mediannoise              = data.medianNoiseRemoval
statistics         = stats       = stat       = data._DataAna__statistics
computeRms         = maprms                   = data.computeRmsFromMap
medianBaseline     = medianBase  = medianbase = data.medianBaseline


# Mapping related
####################################################################
mapping         = doMap         = domap       = data.doMap
smoothBy        = smooth                      = data.smoothMap
display         = mapdisp       = mapdisplay  = data.displayMap
chanMap         = ChanMap       = chanmap     = data.chanMap
computeWeight   = computeweight = weight      = data.computeWeight
showMap                                       = data.showMap
mapsum          = mapSum                      = BoaMapping.mapSum
mapsum2         = mapSum2                     = BoaMapping.mapSum2
mapsumfast      = mapSumFast                  = BoaMapping.mapsumfast
zoom                                          = data.zoom
setValPoly      = setValuesPolygon            = BoaMapping.setValuesPolygon

# Plotting related
####################################################################
plotArray            = plotarray                      = data.BolometerArray.plotArray
plotAzimuth          = azimuth   = azim = az          = data.ScanParam.plotAzimuth
plotElevation        = elevation = elev = el          = data.ScanParam.plotElevation
plotAzEl             = azel                           = data.ScanParam.plotAzEl
plotAzimuthOffset    = azimuthOffset = azimoff = azo  = data.ScanParam.plotAzimuthOffset
plotElevationOffset  = elevationOffset = eleoff = elo = data.ScanParam.plotElevationOffset
plotAzElOffset       = azeloff = azelo                = data.ScanParam.plotAzElOffset
plotAzElSpeed        = azelspeed = azelsp             = data.ScanParam.plotAzElSpeed
plotAzElAcceleration = azelaccel = azelac             = data.ScanParam.plotAzElAcceleration
plotCorrel           = plotcorrel = plotcor = plotCor = data.plotCorrel
plotFFT                                               = data.plotFFT
plotMean             = plotmean                       = data.plotMean
plotRms              = plotrms                        = data.plotRms
plotMeanChan         = plotmeanchan                   = data.plotMeanChan
plotRmsChan          = plotrmschan                    = data.plotRmsChan
plotSubscan          = plotSub                        = data.ScanParam.plotSubscan
plotSubscanOffsets   = plotSubOff                     = data.ScanParam.plotSubscanOffsets
signal               = signa = sign = sig             = data.signal


# Pointing and Focus related
####################################################################
solvePointingOnMap = solvepointing = solvepoint = solvepoin = solvepoi = data.solvePointingOnMap
solveFocus         = solvefocus = solveFoc = solvefoc                  = data.solveFocus


# Data handling related
####################################################################
setCurrChanList   = channels = channel = chan = data.BolometerArray.setCurrChanList
checkChanList     = checkChannels = checkChan = data.BolometerArray.checkChanList
printCurrChanList = printChannels = printChan = data.BolometerArray.printCurrChanList
getPixel          = getPix                    = data.getPixelFromMap


# Bogli Related
####################################################################
selectDev = device = dev      = DeviceHandler.selectDev
clear     = cle = cl          = Plot.clear
openDev   = op                = DeviceHandler.openDev
closeDev  = close = clo = cls = DeviceHandler.closeDev
pointSize = BogliConfig.pointSize
poly      = defPoly = defPolygon = Plot.defPolygon
savePoly  = savePolygon       = Plot.savePolygon
loadPoly  = loadPolygon       = Plot.loadPolygon
plotPoly  = plotPolygon       = Plot.plotPolygon

# SNF related
####################################################################
correlatedNoiseRemoval = cnr    = CNR  = data.correlatedNoiseRemoval
pca = PCA = corrpca = corrPCA = data.corrPCA


# Save & restore
####################################################################
dumpData    = dump    = data.dumpData
restoreData = restore = data.restoreData


# Flag related
####################################################################
flag                                       = data.flag
flagMJD                                    = data.flagMJD
flagLon                                    = data.flagLon
flagRms                                    = data.flagRms
flagSpeed      = flagS                     = data.flagSpeed
flagAccel      = flagA                     = data.flagAccel
flagFractionRms = flagFRms                 = data.flagFractionRms
flagPosition   = flagPos                   = data.flagPosition
flagSubscan    = flagSub                   = data.flagSubscan  # do we want this???    
flagChannels   =   flagC =   flagCh =  fCh = data.flagChannels
unflag                                     = data.unflag
unflagChannels = unflagC = unflagCh = ufCh = data.unflagChannels
despike        = dspike                    = data.despike
itDespike      = iterativeDespike          = data.iterativeDespike


# I/O related
####################################################################
findInDir        = find    = fd        = BoaDir.findInDir
removeScans      = remove  = rs        = BoaDir.removeScans 
selectInDir      = select  = slt       = BoaDir.selectInDir
listInDir        = ils     = inls      = BoaDir.listInDir
listInDir2       = ils2    = inls2     = BoaDir.listInDirNsamp
setInDir         = indir   = ind       = BoaDir.setInDir
setInFile        = infile  = inf       = BoaDir.setInFile
setOutFile       = outfile = outf      = BoaDir.setOutFile
setOutDir        = outdir  = outd      = BoaDir.setOutDir
setProjectID     = setproj = proj      = BoaDir.setProjectID
setDate          = setdate = date      = BoaDir.setDate
resetCurrentList = resetCurrList = rls = BoaDir.resetCurrentList
read                                   = data.read


# Misc
####################################################################
rcp       = readRCP     = readRCPfile = data.BolometerArray.readRCPfile
restore   = restoreData               = data.restoreData
updateRCP                             = data.BolometerArray.updateRCP
flagRCP                               = data.flagRCP
flatfield = flat                      = data.flatfield
correctOpacity = opacity   = opac     = data.correctOpacity
flattenFreq                           = data.flattenFreq
blankFreq                             = data.blankFreq
addSource = addSourceModel            = data.addSource  # defined in BoaMapping
#saveMambo = mambo = data.saveMambo


# MessageHandler related
####################################################################
#mess=messages.setMess
#messweight=messweigh=messweig=messwei=messwe=messw=BoaB.MessHand.setMaxWeight
messweight  =  messw  = data.MessHand.setMaxWeight
