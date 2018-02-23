import copy
from Numeric import *
import math
import BoaMapping
import BoaPointing
import Utilities
import BoaFlagHandler
from Bogli import DeviceHandler, Plot, MultiPlot
from BoaDataEntity import Telescope, BolometerArray, ScanParameter
from slalib import *
        
def pointSource():
    
    modelImage = BoaMapping.Image()

    # Source position
    azimuth_source   = 0.0
    elevation_source = 0.0
    size_source      = 3.0 # Point source = 3 ; Jupiter = 40
    flux_source      = 1.0
    size_map         = 2*400.0 # arcsec
    
    pixelsize_model = 3 # arcsec
    
    modelImage.computeWCS(pixelsize_model,\
                          minmax=array([-0.5,0.5,-0.5,0.5],Float)*size_map)

    p_source = array([0.0,0.0,0.0,\
                              flux_source,azimuth_source,elevation_source,\
                              size_source,size_source,0.0])
    
    dimX,dimY       = modelImage.WCS['NAXIS1'], modelImage.WCS['NAXIS2']
    
    azimuth_model, elevation_model = modelImage.physicalCoordinates()
    modelImage.data = Utilities.modelBaseEllipticalGaussian(p_source,[azimuth_model,elevation_model])
    

    return modelImage

def giorgioImage():
    import Image
    modelImage = BoaMapping.Image()
    size_map         = 2*400.0 # arcsec
    pixelsize_model = 3 # arcsec
    modelImage.computeWCS(pixelsize_model,\
                          minmax=array([-0.5,0.5,-0.5,0.5],Float)*size_map,\
                          )
    dimX,dimY       = modelImage.WCS['NAXIS1'], modelImage.WCS['NAXIS2']
    
    im = Image.open("foto4.jpg")
    im = im.convert('L')
    im = im.rotate(-90)
    im = im.resize((dimY,dimX))
    im = reshape(list(im.getdata()),(dimX,dimY))
    modelImage.data = im

    return modelImage

def rasterPositions(azimuth_raster,elevation_raster,\
                    azimuth_source,elevation_source):
    global timeSubscan,nSubscans,sampling_rate,lst_0
    scan_time = nSubscans*timeSubscan
    
    print '-------------------'
    print 'Simulating raster :'
    print '    # subscans : ', nSubscans
    print ' scan duration : ', scan_time,' s'
    print ' sampling rate : ', sampling_rate, ' Hz'

    # The output struture
    outParam              = ScanParameter()

    Nint                  = int(sampling_rate*scan_time)
    ndumps                = int(sampling_rate*timeSubscan) # total number of samples per subscan
    integnum              = arrayrange(ndumps*nSubscans)
                                       
    outParam.NInt         = Nint
    outParam.NObs         = nSubscans
    outParam.SubscanNum   = arrayrange(nSubscans)+1
    outParam.SubscanNum   = outParam.SubscanNum.tolist()
    outParam.SubscanTime  = zeros(nSubscans,'f')
    outParam.SubscanIndex = zeros((2,nSubscans))

    outParam.LST     = arrayrange(Nint,typecode=Float)/\
                       (Nint+1)*scan_time + lst_0
    mjd0 = 54331.5
    outParam.MJD     = arrayrange(Nint,typecode=Float)/\
                       (Nint+1)*scan_time/86400. + mjd0
    
    # Origin of the circle
    outParam.Az     = zeros(Nint,'f')
    outParam.El     = zeros(Nint,'f')
    outParam.LonOff = zeros(Nint,'f')
    outParam.LatOff = zeros(Nint,'f')
    outParam.Baslon = zeros(Nint,'f')
    outParam.Baslat = zeros(Nint,'f')

    outParam.FlagHandler = BoaFlagHandler.createFlagHandler(zeros(Nint,Int32))

    # For each subscan, fix the position of the offsets at the raster position
    for sub in range(nSubscans):

        # Az/El in deg
        outParam.Az[sub*ndumps:(sub+1)*ndumps] = azimuth_source   + azimuth_raster[sub]/3600.
        outParam.El[sub*ndumps:(sub+1)*ndumps] = elevation_source + elevation_raster[sub]/3600.

        # Lon/Lat in deg
        outParam.LonOff[sub*ndumps:(sub+1)*ndumps] = azimuth_raster[sub]/3600.
        outParam.LatOff[sub*ndumps:(sub+1)*ndumps] = elevation_raster[sub]/3600.
        
        outParam.SubscanTime[sub]    = outParam.LST[sub*ndumps]
        outParam.SubscanIndex[:,sub] = array([sub*ndumps,(sub+1)*ndumps])

    return outParam
    

def raster_otf():
    global nSubscans,azimuth_raster,elevation_raster,specificParam
    azimuth_raster = zeros((nSubscans),'f')
    elevation_raster = (arrayrange(nSubscans,typecode='f') -
                        (nSubscans-1)/2.) * specificParam['Ystep']
    
def telescopePositions_otf(outParam,specificParam):
    global lst_0,ra0,dec0,azimuth_raster,elevation_raster
    print '--------------------'
    print 'Simulating OTF :'
    print ' Xlen, Xstep   :', specificParam['Xlen'],specificParam['Xstep'], ' arcsec'
    print ' Ylen, Ystep   :', specificParam['Ylen'],specificParam['Ystep'], ' arcsec'
    print ' linear speed  :', specificParam['Xlen']/timeSubscan/60.,' arcmin/s'

    ra0rad = ra0 * pi/180.
    de0rad = dec0 * pi/180.
    # Precess to date of obs.
    epoch = 2007.64  # sometime late August 2007
    ra0rad,de0rad = sla_preces ('FK5',2000.,epoch,ra0rad,de0rad)
    phi = (defineAPEX().Latitude)*pi/180.
    lst = outParam.LST

    for sub in range(outParam.NObs):

        # Find out time in sec for this subscan
        startSubscan = outParam.SubscanIndex[0,sub]
        endSubscan   = outParam.SubscanIndex[1,sub]
        subscan_time = outParam.LST[startSubscan:endSubscan]-outParam.LST[startSubscan]
        ndumps       = len(subscan_time)
        azimuth_array   = (arrayrange(ndumps,typecode='f') * specificParam['Xstep']
                           - specificParam['Xlen'] /2.) / 3600.
        elevation_array = zeros((ndumps),'f')

        # Add the offset to the output structure
        for j in range(startSubscan,endSubscan):
            lst = outParam.LST[j]
            ha = lst / 86400. * 2.*pi - ra0rad
            az0, el0 = sla_e2h(ha,de0rad,phi)
            outParam.Az[j] = az0 * 180./pi + azimuth_array[j-startSubscan]
            outParam.El[j] = el0 * 180./pi + elevation_array[j-startSubscan]
        outParam.Az[startSubscan:endSubscan] += azimuth_array.astype('f')
        #outParam.El[startSubscan:endSubscan] += elevation_array.astype('f')
        outParam.El[startSubscan:endSubscan] += array(elevation_raster[sub]/3600.,'f')
        
        outParam.LonOff[sub*ndumps:(sub+1)*ndumps] += azimuth_array.astype('f')
        outParam.LatOff[sub*ndumps:(sub+1)*ndumps] += elevation_array.astype('f')

    return outParam

def telescopePositions_circle(outParam,specificParam):

    print '--------------------'
    print 'Simulating circles :'
    print ' circle radius   :', specificParam['radius'], ' arcsec'
    print 'rotation speed   :', specificParam['rotationSpeed'], ' deg/s'
            
    for sub in range(outParam.NObs):
        
        # Find out time in sec for this subscan
        startSubscan = outParam.SubscanIndex[0,sub]
        endSubscan   = outParam.SubscanIndex[1,sub]
        subscan_time = outParam.LST[startSubscan:endSubscan]-outParam.LST[startSubscan]

        dumps = subscan_time*2*pi*specificParam['rotationSpeed']/360.
        ndumps = len(dumps)
        
        azimuth_array   = specificParam['radius']*cos(dumps) / 3600.
        elevation_array = specificParam['radius']*sin(dumps) / 3600.

        # Add the offset to the output structure
        outParam.Az[startSubscan:endSubscan] += azimuth_array.astype('f')
        outParam.El[startSubscan:endSubscan] += elevation_array.astype('f')

        outParam.LonOff[sub*ndumps:(sub+1)*ndumps] += azimuth_array.astype('f')
        outParam.LatOff[sub*ndumps:(sub+1)*ndumps] += elevation_array.astype('f')

    return outParam

def telescopePositions_spiral(outParam,specificParam):

    print '-------------------'
    print 'Simulating spiral :'
    print ' start radius  :', specificParam['radius'], ' arcsec'
    print '  rotat speed  :', specificParam['rotationSpeed'], ' deg/sec'
    print ' radial speed  :', specificParam['radiusSpeed'], ' arcsec /sec'
    rotSpeed = specificParam['rotationSpeed'] *2*pi/360.  # in rad/sec
    print 'min lin speed  :', specificParam['radius']*rotSpeed/60.,' arcmin /sec'
    print 'max lin speed  :', (specificParam['radius']+timeSubscan*specificParam['radiusSpeed']) \
                                *rotSpeed/60.,' arcmin /sec'
    
    
    for sub in range(outParam.NObs):

        # Find out time in sec for this subscan
        startSubscan = outParam.SubscanIndex[0,sub]
        endSubscan   = outParam.SubscanIndex[1,sub]
        subscan_time = outParam.LST[startSubscan:endSubscan]-outParam.LST[startSubscan]

        # the parameter of the spiral (in degree)
        theta = specificParam['rotationSpeed']/180*pi*subscan_time
        radius = (specificParam['radius']+specificParam['radiusSpeed']*subscan_time)/3600 

        azimuth_array   = radius*cos(theta)
        elevation_array = radius*sin(theta)
        
        # Add the offset to the output structure
        outParam.Az[startSubscan:endSubscan] += azimuth_array.astype('f')
        outParam.El[startSubscan:endSubscan] += elevation_array.astype('f')
        
        outParam.LonOff[startSubscan:endSubscan] += azimuth_array.astype('f')
        outParam.LatOff[startSubscan:endSubscan] += elevation_array.astype('f')

    return outParam


def telescopePositions_lissajou(outParam,specificParam):

    print '---------------------'
    print 'Simulating lissajou :'
    print ' az amplitude  : ',specificParam['az_amplitude'], ' arcsec'
    print ' az period     : ',specificParam['az_period'], ' s'
    print ' el amplitude  : ',specificParam['el_amplitude'], ' arcsec'
    print ' el period     : ',specificParam['el_period'], ' s'

    for sub in range(outParam.NObs):
        
        # Find out time in sec for this subscan
        startSubscan = outParam.SubscanIndex[0,sub]
        endSubscan   = outParam.SubscanIndex[1,sub]
        subscan_time = outParam.LST[startSubscan:endSubscan]-outParam.LST[startSubscan]

        # Lissajou parameters in degree
        theta_azimuth   = 2*pi*subscan_time/specificParam['az_period']
        theta_elevation = 2*pi*subscan_time/specificParam['el_period']        
        azimuth_array   = specificParam['az_amplitude']*cos(theta_azimuth)/3600
        elevation_array = specificParam['el_amplitude']*sin(theta_elevation)/3600

        # Add the offset to the output structure 
        outParam.Az[startSubscan:endSubscan] += azimuth_array.astype('f')
        outParam.El[startSubscan:endSubscan] += elevation_array.astype('f')
        
        outParam.LonOff[startSubscan:endSubscan] += azimuth_array.astype('f')
        outParam.LatOff[startSubscan:endSubscan] += elevation_array.astype('f')

    return outParam


def timeSeries(modelImage, BolometerArray, ScanParam):

    print '------------------------'
    print 'Simulating Observation :'
    print BolometerArray
    print ScanParam
    print ''

    beamSize = BolometerArray.BeamSize
    nBeams   = BolometerArray.NChannels
    px,py    = BolometerArray.Offsets
    nDumps   = shape(ScanParam.LST)[0]
    
    signal   = zeros((nDumps,nBeams), Float)

    # Convolve with the beam
    # Same beam for every channel so out of the loop, otherwhise one
    # should convolve with a different beam for each channel

    simulated_image = copy.copy(modelImage)
    simulated_image.smoothBy(beamSize)
    
    myTiming = Utilities.Timing()
    myTiming.setIter(nBeams)
    myTiming.resetTime()
    
    for chan in range(nBeams):
        
        azimuth_beam   = ScanParam.LonOff*3600 - px[chan]
        elevation_beam = ScanParam.LatOff*3600 - py[chan]
        
        i,j            = simulated_image.wcs2pix(azimuth_beam,elevation_beam)

        for dump in range(nDumps):
            if 0 < i[dump] < modelImage.WCS['NAXIS1'] and \
                   0 < j[dump] < modelImage.WCS['NAXIS2']:
                signal[dump,chan] = simulated_image.data[int(i[dump]),int(j[dump])]

        myTiming.timeLeft(iter=chan)

    return signal


def pack(signal,BolometerArray,ScanParameters):
    # Put everything into one Boa object
    
    outPut = BoaPointing.Point()
    outPut.BolometerArray = Laboca
    outPut.ScanParam = SimulationParam

    nDumps = signal.shape[0]/nSubscans

    outPut.Data        = signal
    outPut.DataWeights = ones((nDumps*nSubscans,Laboca.NChannels),'f')    
    for i in xrange(SimulationParam.NInt):
        outPut.DataWeights[i,::] = BolometerArray.Weights
    outPut.FlagHandler = \
        BoaFlagHandler.createFlagHandler(zeros((nDumps*nSubscans,Laboca.NChannels),'i'))
    
    outPut.ChanRms  = zeros((3,Laboca.NChannels),Float)
    outPut.ChanMean = zeros((3,Laboca.NChannels),Float)
    outPut.ChanMed  = zeros((3,Laboca.NChannels),Float)
    
    outPut.ChanRms_s  = zeros((3,Laboca.NChannels,nSubscans),Float)
    outPut.ChanMean_s = zeros((3,Laboca.NChannels,nSubscans),Float)
    outPut.ChanMed_s  = zeros((3,Laboca.NChannels,nSubscans),Float)
    
    return outPut



# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

def defineAPEX():
    APEX = Telescope()
    APEX.Name = 'APEX'
    APEX.Diameter = 12.0
    APEX.Longitude =  67.+45/60+33.0/60/60
    APEX.Latitude  = -23.+00/60+20.8/60/60
    APEX.Elevation = 5100.
    return APEX

def defineArray():
    # Define the array
    APEX = defineAPEX()

    Laboca = BolometerArray()
    Laboca.Telescope = APEX
    Laboca.FeBe = 'LABOCA-ABBA2'
    Laboca.readRCPfile('LABOCA-centred.rcp')     # 320 channels
    Laboca.updateRCP('master-laboca-may07.rcp')  # 257 good channels
    Laboca.EffectiveFrequency = 345.0e9
    Laboca._BolometerArray__computeBeamSize()
    Laboca.flag([201,227])
    cross = [87,101,103,106,110,139,143,149,163,183,241,255] 
    Laboca.flag(cross)

    # Read in channel weights
    Laboca.Weights = zeros((Laboca.NChannels),'f')
    Laboca.Gain_bak = copy.copy(Laboca.Gain)
    f = file('rcp/pixel-11121.dat')
    rl = f.readlines()
    for l in rl[1::]:
        tmp = l.split()
        Laboca.Gain[int(tmp[0])-1] = float(tmp[2])
        Laboca.Weights[int(tmp[0])-1] = float(tmp[3])
    f.close()

    bad = []
    for i in range(Laboca.NChannels):
        if Laboca.Weights[i] == 0.:
            bad.append(i+1)
    Laboca.flag(bad)
    return(Laboca)

def defineArraySZ():
    # Define the array
    APEX = defineAPEX()

    ASZCa = BolometerArray()
    ASZCa.Telescope = APEX
    ASZCa.FeBe = 'BOLOSZ-SZACBE'
    ASZCa.readRCPfile('ASZCa_April2007_positions.rcp')     # 331 channels, numberes are 1-based
    ASZCa.EffectiveFrequency = 150.0e9
    ASZCa.Telescope.Diameter *= 2./3.  # they use only 2/3 of the dish...
    ASZCa._BolometerArray__computeBeamSize()

    # Read in channel weights
    ASZCa.Weights = zeros((ASZCa.NChannels),'f')
    f = file('rcp/aszca_nefd.txt')   # Here numbers are zero-based
    rl = f.readlines()
    f.close()
    all = range(1,332)
    for l in rl[3::]:
        tmp = l.split()
        ASZCa.Weights[int(tmp[0])] = 1./float(tmp[2])**2
        all.remove(int(tmp[0])+1)
    ASZCa.flag(all)
    return(ASZCa)


################################################################################
# Functions to define some Raster patterns

def raster1():
    # only one point, no offset
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 30
    azimuth_raster   = zeros((30),'f')
    elevation_raster = zeros((30),'f')
    #azimuth_raster   = array([0],'f')
    #elevation_raster = array([0],'f')

def raster4():
    # 4 circles centered
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 4
    azimuth_raster   = array([-1,1,1,-1],'f')*18
    elevation_raster = array([-1,1,-1,1],'f')*18

def hexagon1():
    # hexagon with 7 positions
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 7
    azimuth_raster   = array([0,  0, 0, sqrt(3)/2, -sqrt(3)/2,sqrt(3)/2, -sqrt(3)/2],'f')*18*25
    elevation_raster = array([0, -1, 1,       0.5,        0.5,     -0.5,      -0.5 ],'f')*18*25

def hexagon2():
    # hexagon with 14 subscans
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 14
    azimuth_raster1   = array([0,  0, 0, sqrt(3)/2, -sqrt(3)/2,sqrt(3)/2, -sqrt(3)/2],'f')*18*36
    elevation_raster1 = array([0, -1, 1,       0.5,        0.5,     -0.5,      -0.5 ],'f')*18*36
    azimuth_raster2= array([0,  0, 0, sqrt(3)/2, -sqrt(3)/2,sqrt(3)/2, -sqrt(3)/2],'f')*18*36+18*7
    elevation_raster2 = array([0, -1, 1,       0.5,        0.5,     -0.5,      -0.5 ],'f')*18*36+18*7
    angle=30.
    rot_mat = array([[cos(angle*pi/180), sin(angle*pi/180)],[-sin(angle*pi/180),cos(angle*pi/180)]])
    new = matrixmultiply(rot_mat,array([azimuth_raster1,elevation_raster1])/2)
    azimuth_raster2   = new[0]
    elevation_raster2 = new[1]
    azimuth_raster = ravel([azimuth_raster1,azimuth_raster2])
    elevation_raster = ravel([elevation_raster1,elevation_raster2])

def hexagon3():
    # hexagon Deep Field
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 14
    azimuth_raster1   = array([0,  0, 0, sqrt(3)/2, -sqrt(3)/2,sqrt(3)/2, -sqrt(3)/2],'f')*18*35
    elevation_raster1 = array([0, -1, 1,       0.5,        0.5,     -0.5,      -0.5 ],'f')*18*35
    azimuth_raster2   = array([sqrt(3)/2, -sqrt(3)/2, sqrt(3)/2, -sqrt(3)/2, sqrt(3), sqrt(3), sqrt(3)],'f')*18*35
    elevation_raster2 = array([1.5,       1.5,        -1.5,     -1.5, 1,      -1, 0 ],'f')*18*35
    azimuth_raster = ravel([azimuth_raster1,azimuth_raster2])
    elevation_raster = ravel([elevation_raster1,elevation_raster2])

def square1(spacing=180.):
    # 3x3 raster
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 9
    azimuth_raster   = array([-1,0,1,-1,0,1,-1,0,1],'f')*spacing
    elevation_raster = array([-1,-1,-1,0,0,0,1,1,1],'f')*spacing
    
def square2(spacing=180.):
    # 3x3 raster + 2x2 raster
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 30
    nSubscans        = 13
    azimuth1         = array([-1,0,1,-1,0,1,-1,0,1],'f')*spacing
    elevation1       = array([-1,-1,-1,0,0,0,1,1,1],'f')*spacing
    azimuth2         = array([-1,1,1,-1],'f')*spacing/2.
    elevation2       = array([-1,-1,1,1],'f')*spacing/2.
    azimuth_raster   = concatenate((azimuth1,azimuth2))
    elevation_raster = concatenate((elevation1,elevation2))
    
def square3(spacing=180.):
    # 4x4 raster + 3x3 raster
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40
    nSubscans        = 25
    azimuth1         = array([-1.5,-0.5,0.5,1.5])*spacing
    elev1            = array(4*[-1.5])*spacing
    elev2            = array(4*[-0.5])*spacing
    elev3            = array(4*[0.5])*spacing
    elev4            = array(4*[1.5])*spacing
    azimuth2         = array([-1,0,1,-1,0,1,-1,0,1],'f')*spacing
    elevation2       = array([-1,-1,-1,0,0,0,1,1,1],'f')*spacing
    azimuth_raster   = concatenate((azimuth1,azimuth1,azimuth1,azimuth1,azimuth2))
    elevation_raster = concatenate((elev1,elev2,elev3,elev4,elevation2))

def strip1():
    # Long strip - 2deg x 15'
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 25     # sec
    nSubscans        = 27
    azimuth_raster   = []
    elevation_raster = []
    for i in range(nSubscans):
        azimuth_raster.append(-3900.+300.*i)
        elevation_raster.append((-1)**(i+1) * 300.)

def strip2():
    # Long strip - 2deg x 0.5deg
    global nSubscans,timeSubscan,azimuth_raster,elevation_raster
    timeSubscan      = 40     # sec
    nSubscans        = 19
    azimuth_raster   = []
    elevation_raster = []
    for i in range(nSubscans):
        azimuth_raster.append(-4200.+460.*i)
        elevation_raster.append((-1)**(i+1) * 400.)

################################################################################
# Definition of specific parameters
def param_otf1():
    global specificParam,timeSubscan,sampling_rate,nSubscans
    specificParam = {
        'Xlen': 300.,\
        'Xstep': 6.,\
        'Ylen': 240.,\
        'Ystep': 12. }
    timeSubscan = int(specificParam['Xlen'] / specificParam['Xstep'] + 1) * 1./sampling_rate
    nSubscans   = int(specificParam['Ylen'] / specificParam['Ystep'] + 1)

def param_otf2():
    global specificParam,timeSubscan,sampling_rate,nSubscans
    specificParam = {
        'Xlen': 1800. * 1.6,\
        'Xstep': 3.,\
        'Ylen': 1800.,\
        'Ystep': 20. }
    timeSubscan = int(specificParam['Xlen'] / specificParam['Xstep'] + 1) * 1./sampling_rate
    nSubscans   = int(specificParam['Ylen'] / specificParam['Ystep'] + 1)

def param_otf3():
    global specificParam,timeSubscan,sampling_rate,nSubscans
    specificParam = {
        'Xlen': 24,\
        'Xstep': 12.,\
        'Ylen': 20.,\
        'Ystep': 20. }
    timeSubscan = int(specificParam['Xlen'] / specificParam['Xstep'] + 1) * 1./sampling_rate
    nSubscans   = int(specificParam['Ylen'] / specificParam['Ystep'] + 1)

def param_circle():
    global specificParam
    # laboca beam is 18.2" on APEX
    specificParam = {'radius':  5.*18.2,
                     'rotationSpeed': 90.}

def param_spiral1():
    global specificParam
    specificParam = {
        'rotationSpeed': 120.,\
        'radius':  30.,\
        'radiusSpeed': 2. }
    
def param_spiral2():
    global specificParam
    specificParam = {
        'rotationSpeed': 35.,\
        'radius':  240.,\
        'radiusSpeed': 2. }

def param_spiral3():
    global specificParam
    specificParam = {
        'rotationSpeed': 60.,\
        'radius':  120.,\
        'radiusSpeed': 1.5 }

def param_spiral4():
    global specificParam
    specificParam = {
        'rotationSpeed': 120.,\
        'radius':  18./2,\
        'speed': 18./6 }

def param_lissajou():
    global specificParam
    specificParam = {
        'az_amplitude': 1000.,\
        'el_amplitude': 1000.,\
        'az_period': 2.,\
        'el_period': 2.*sqrt(2) }


################################################################################
#
# Main program
#

# Global parameters
#sampling_rate    = 25.0   # sample rate Hz
sampling_rate    =  50.0   # sample rate Hz

# Choose the input image
#sourceMap = giorgioImage
sourceMap = pointSource

# Source Position :
#lst_0            = 80721.
ra0,dec0         = 35.9,-4.
#lst_0            = 8640.
lst_0            = 8640. - 15.*60.  - 1.5* 3600.

# compute Az, El at start
ra0rad = ra0 * pi/180.
de0rad = dec0 * pi/180.
# Precess to date of obs.
epoch = 2007.64  # sometime late August 2007
ra0rad,de0rad = sla_preces ('FK5',2000.,epoch,ra0rad,de0rad)
phi = (defineAPEX().Latitude)*pi/180.
# convert RA0, Dec0 to Az0, El0
ha = lst_0 / 86400. * 2.*pi - ra0rad
az0, el0 = sla_e2h(ha,de0rad,phi)

azimuth_source   = az0 * 180./pi
elevation_source = el0 * 180./pi

# print " ***** RA,Dec = ",ra0,dec0
# print " ***** Az,El  = ",azimuth_source,elevation_source
# print " ***** LST(0) = ",lst_0

########################################
# Here we choose the scan pattern and parameters
# Example: OTF
telescopePositions = telescopePositions_otf
param_otf2()
raster_otf()

# # Example: raster of circles
#telescopePositions = telescopePositions_circle
#raster1()
#param_circle()

# # Example: raster of spirals
#telescopePositions = telescopePositions_spiral
#hexagon1()
#param_spiral1()

# Example: raster of spirals
#telescopePositions = telescopePositions_spiral
#strip1()
#param_spiral2()

# Example: raster of spirals
#telescopePositions = telescopePositions_spiral
#square2(500.)
#param_spiral3()

# Example: raster of spirals
#telescopePositions = telescopePositions_spiral
#square3(600.)
#param_spiral3()


########################################

# Global Scan parameters
timeScan = timeSubscan*nSubscans   # seconds
#Laboca   = defineArray()
Laboca   = defineArraySZ()

# Doing the simulation :
print ''
print 'Simulation started'

# Compute the model image
modelImage = sourceMap()

# Compute the telescope positions :
SimulationParam = rasterPositions(azimuth_raster,elevation_raster,\
                                  azimuth_source,elevation_source)
SimulationParam = telescopePositions(SimulationParam, specificParam)

# print " **** Az = ",SimulationParam.Az
# print " **** El = ",SimulationParam.El

SimulationParam.Object  = 'simulation'
SimulationParam.ScanNum = 1
SimulationParam.Frames = 'EQEQHO'
SimulationParam.Telescope = Laboca.Telescope
SimulationParam.AzOff = SimulationParam.LonOff
SimulationParam.ElOff = SimulationParam.LatOff
 
# for i in range(SimulationParam.NInt):
#     SimulationParam.Baslon[i], SimulationParam.Baslat[i] = sla_dh2e(SimulationParam.Az[i]/180*pi,\
#                                                                     SimulationParam.El[i]/180*pi,\
#                                                                     Laboca.Telescope.Latitude/180*pi)

# SimulationParam.Baslon = SimulationParam.Baslon/pi*180+SimulationParam.LST*360/86400
# SimulationParam.Baslat = SimulationParam.Baslat/pi*180

# Simulate the observation
signal = timeSeries(modelImage, Laboca, SimulationParam)

# Put everything into one BoaPoint object
myMap = pack(signal,Laboca,SimulationParam)

# Do some plotting
# myMap.doMap(oversamp=3)
# myMap.Map.smoothBy(myMap.BolometerArray.BeamSize/2.)
# myMap.Map.display(coverage=1,style='heat',caption='Coverage')
# #myMap.BolometerArray.plotArray(overplot=1)
# myMap.ScanParam.plotAzElOffset(overplot=1)
# raw_input()
# myMap.Map.display(weight=1,style='heat',caption='Weights')
# raw_input()

#####
# some more parameters for ASZCa
myMap.ScanParam.Coord = (ra0,dec0)
myMap.ScanParam.computeRa0De0()


myMap.ScanParam.computeRaDec() 

# print " *** all RA,Dec = ",myMap.ScanParam.RA,myMap.ScanParam.Dec
# print " *** all LST = ",myMap.ScanParam.LST

myMap.ScanParam.computeParAngle()
myMap.doMap(system='EQ')
myMap.Map.display(weight=1,style='heat',caption='Weights')
