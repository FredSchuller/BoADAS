#######################################################################################
#
# computeRCP.py
#
# script to determine array parameters (RCP) from a beam map
#
# FS - 2007/3/07 - 2007/5/20
#    - 2007/10/04: remove call to ppgplot
#
#######################################################################################
from fortran import fStat

def iterat(chans=[],scanNum=0,oversamp=3,sizeX=[-300,300],sizeY=[-300,300],\
           data=None,flux=1.,radius=50):
    """
    This function solves for pointing (i.e. fits a 2D Gaussian) on each individual
    channel map, and stores the results in arrayParamOffsets, as well as in a text
    file (offsets, widths, angle, integrated and peak fluxes)
    """
    newBadChan = []
    newGoodChan = []
    ox=[]
    oy=[]
    wx=[]
    wy=[]
    gain=[]
    pa=[]
    arrayParamOffsets = []

    dummy = -1
    beam = data.BolometerArray.BeamSize
    gainnorm = 0.0
    
    

    for c in chans:
        print "Channel = ",c
        num = data.BolometerArray.getChanIndex(c)[0]
        Plot.nextpage()
        data.doMap(c,noPlot=1,oversamp=oversamp,sizeX=sizeX,sizeY=sizeY)
        try:
            #data.solvePointingOnMap(radius=radius,plot=0)
            data.solvePointing(c,radius=radius,plot=1)
            ##raw_input()
            #data.showPointing(limitsZ=[-0.1,0.7])
        except:
            #output.write("%3i no fit found\n"%(c))
            newBadChan.append(c)
            ox.append(dummy)
            oy.append(dummy)
            wx.append(dummy)
            wy.append(dummy)
            gain.append(dummy)
            pa.append(dummy)
            
        Plot.xyout(40,40,str(c))
        result = data.PointingResult
        if result == -1:
            #output.write("%3i no fit found\n"%(c))
            newBadChan.append(c)
            ox.append(dummy)
            oy.append(dummy)
            wx.append(dummy)
            wy.append(dummy)
            gain.append(dummy)
            pa.append(dummy)
        else:
            offX = result['gauss_x_offset']['value']
            offY = result['gauss_y_offset']['value']
            widX = result['gauss_x_fwhm']['value']
            widY = result['gauss_y_fwhm']['value']
            ang  = result['gauss_tilt']['value']
            peak = result['gauss_peak']['value']
            # relative gain (we will have to divide by g for flat-fielding)
            g = peak/1.0
            result['gain'] = g 
            offX = data.BolometerArray.Offsets[0,num] - offX
            offY = data.BolometerArray.Offsets[1,num] - offY

            if g>0.05 and widX<6*beam and widY<6*beam:
                
                arrayParamOffsets.append({'channel' : c,\
                                      'result'  : result})
                
                ox.append(offX)
                oy.append(offY)
                wx.append(widX)
                wy.append(widY)
                gain.append(g)
                pa.append(ang)
                newGoodChan.append(c)
                gainnorm = gainnorm+g
            else:
                newBadChan.append(c)
                ox.append(dummy)
                oy.append(dummy)
                wx.append(dummy)
                wy.append(dummy)
                gain.append(dummy)
                pa.append(dummy)

    nrgood = len(newGoodChan)
    if nrgood > 0:
        gainnorm = gainnorm/len(newGoodChan)
        for i in range(len(gain)):
            gain[i]=gain[i]/gainnorm
    
    filename = "beams_%i_%s.txt"%(oversamp,scanNum)
    output = file(filename,'w')
    for i in range(len(chans)):
        if gain[i]>0:
            text = "%3i %10.3f %10.3f %5.2f %5.2f %5.2f %5.3f\n" % \
                       (chans[i],ox[i],oy[i],wx[i],wy[i],pa[i],gain[i])
            output.write(text)
        else:
            output.write("%3i no fit found\n"%(chans[i]))
    output.close()
    if newBadChan:
    	flagC(newBadChan)
    print gainnorm
    if nrgood > 0:
        for param in arrayParamOffsets:
            param['result']['gain'] /= gainnorm

    data.arrayParamOffsets = arrayParamOffsets   # store results in current object
    
    return nrgood
#######################################################################################
def plotRCP(chans,scanNum,num=1,sizeX=[500,-500],sizeY=[-500,500],data=None):
    """
    Plot ellipses showing the beams shapes, using the fits results stored
    in data.arrayParamOffsets
    """
    cap = "LABOCA Beam map - Scan %s"%(scanNum)
    plot([0],[0],limitsX=sizeX,limitsY=sizeY,nodata=1,aspect=1,\
         labelX='Az offset ["]',labelY='El offset ["]',caption=cap)

    nb = len(chans)
    chanList = []
    for i in range(len(data.arrayParamOffsets)):
        chanList.append(data.arrayParamOffsets[i]['channel'])
        
    for i in range(nb):
        num = chanList.index(chans[i])
        theParam = data.arrayParamOffsets[num]
        x = data.BolometerArray.Offsets[0,chans[i]-1]
        y = data.BolometerArray.Offsets[1,chans[i]-1]
        wx = theParam['result']['gauss_x_fwhm']['value']
        wy = theParam['result']['gauss_y_fwhm']['value']
        tilt = theParam['result']['gauss_tilt']['value']
        Forms.ellipse(x,y,wx,wy, tilt*pi/180.,overplot=1)
    if num:
        BogliConfig.xyouttext['color'] = 2
        for i in range(nb):
            x = data.BolometerArray.Offsets[0,chans[i]-1]
            y = data.BolometerArray.Offsets[1,chans[i]-1]
            Plot.xyout(x,y,str(chans[i]))

def plotFromFile(infile='',sizeX=[500,-500],sizeY=[-500,500],cap='',noerase=0):
    """
    Plot beam shapes, as they were stored in a beam*.txt file
    """
    if not infile:
        print "Need a valid input file name... exiting"
        return
    
    try:
        f = file(infile)
    except IOError:
        print "Need a valid input file name... exiting"
        return

    num,Xpos,Ypos = [],[],[]
    fw1,fw2,ang   = [],[],[]
    rl = f.readlines()
    f.close()
    for l in rl:
        tmp = l.split()
        try:
            n = int(tmp[0])
            x = float(tmp[1])
            y = float(tmp[2])
            w1= float(tmp[3])
            w2= float(tmp[4])
            a = float(tmp[5])
            num.append(n)
            Xpos.append(x)
            Ypos.append(y)
            fw1.append(w1)
            fw2.append(w2)
            ang.append(a)
        except:
            pass
    print "Successfully read in %i channels"%(len(num))
    plot([0],[0],limitsX=sizeX,limitsY=sizeY,nodata=1,aspect=1,\
         labelX='Az offset ["]',labelY='El offset ["]',caption=cap,noerase=noerase)
    BogliConfig.xyouttext['color'] = 2
    for i in range(len(num)):
        Forms.ellipse(Xpos[i],Ypos[i],fw1[i],fw2[i],ang[i]*pi/180.,overplot=1)
        Plot.xyout(Xpos[i],Ypos[i],str(num[i]))
    
#######################################################################################
def example():
    # One complete example of use
    #
    # Read a beam map
    scan = 3383
    read(str(scan))

    # Use a good RCP file as starting point, if available
    updateRCP('saturn_3249.rcp')
    
    # Do basic reduction steps
    medianbase()
    flagRms(below=1e-4)  # or whatever corresponds to dead channels
    flagC(317)           # the blind bolometer
    flagC(257)           # a very noisy one
    flagPos(radius=50)   # flag the source
    mediannoise(chanRef=108)
    despike()
    unflag(flag=5)

    # now do a solvepoi() on an individual channel map. Let's say that the result peak is
    # 0.00566, and the source (Saturn) is 2427 Jy. The next command scales all signals to
    # a first reasonable value
    data.Data *= array(2427./0.00566,'f')

    # We will now fit gaussians to all individual channel maps
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-100,100],\
           data=data,flux=2427.,radius=50.)
    close()

    # This will update the RCP in the current object, and write the results out
    # in saturn_3383.rcp
    data.updateArrayParameters("saturn_3383.rcp")

    # Let's see how they look like
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[400,-400],sizeY=[-400,400],data=data)


def oneSZ(scan):
    """
    Process one beam map observed with ASZCa
    INP: (int) scan = scan number
    OUT: 3 files:
           - rcp as text file
           - beam shapes as one ps file
           - individual channel maps + fits in one big ps file
    """

    # Planet fluxes from Astro
    calib = {'Jupiter':2670.,
             'Saturn':612.,
             'Mars':61.}
    
    read(str(scan),febe='BOLOSZ-SZACBE')
    updateRCP('rcpBoa.rcp')
    data.Data *= array(15.,'f')  # 1st approx. to Jy
    data.rebin()    # don't need 100 Hz resolution
    medianbase(subscan=0)
    flagPos(radius=60)   # flag the source
    stat()
    rms = data.getChanListData('rms')
    med = fStat.f_median(rms)
    flagRms(below=med/10.)
    mediannoise(chanRef=102)
    stat()
    rms = data.getChanListData('rms')
    med = fStat.f_median(rms)
    flagRms(above=med*10.)
    #despike()  # can't despike if the source is completely off
    unflag(flag=5)

    # now compute RCPs
    Fnu = calib[data.ScanParam.Object]  # expected flux
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('Fred/beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-250,250],sizeY=[-200,200],\
           data=data,flux=Fnu,radius=100.)
    close()

    # This will update the RCP in the current object, and write the results out
    # in sz_????.rcp
    data.updateArrayParameters("sz_%i.rcp"%(scan))

    # Let's see how they look like
    op('Fred/beamshapes_%i.ps/CPS'%(scan))
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[700,-700],sizeY=[-700,700],data=data)
    close()
    

def beam1():
    scan = 10050
    read(str(scan))
    
    # Use a good RCP file as starting point, if available
    updateRCP('saturn_3288.rcp')
    correctDelay (data,delay=3)
    medianbase()
    mediannoise()
    mediannoise()

    data.Data *= array(1500./0.03,'f')
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-150,150],\
           data=data,flux=1500.,radius=50.)
    close()

    data.updateArrayParameters("laboca-10050.rcp")
    op('beams-10050.ps/CPS')
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[400,-400],sizeY=[-400,400],data=data)
    close()
    
def beam1b():
    # still scan 10050, a few bad fits
    scan = 10050
    read(str(scan))
    bad = [86,139,95,183,204,218,309]
    updateRCP('saturn_3288.rcp')
    correctDelay (data,delay=3)
    medianbase()
    mediannoise()
    mediannoise()
    
    data.Data *= array(1500./0.03,'f')
    op('beams_%i-b.ps/CPS'%(scan))
    iterat(chans=bad,scanNum=scan,oversamp=5,sizeX=[-250,250],sizeY=[-250,250],\
           data=data,flux=1500.,radius=250.)
    close()
    data.updateArrayParameters("laboca-10050-b.rcp")

def beam2():
    scan = 10143
    read(str(scan))
    
    # Use a good RCP file as starting point, if available
    updateRCP('laboca-10050.rcp')
    correctDelay (data,delay=1)
    medianbase()
    mediannoise()
    mediannoise()

    data.Data *= array(400./0.008,'f')
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-150,150],\
           data=data,flux=400.,radius=50.)
    close()

    data.updateArrayParameters("laboca-%s.rcp"%(scan))
    op('beam_shape_%i.ps/CPS'%(scan))
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[400,-400],sizeY=[-400,400],data=data)
    close()

def beam2b():
    # scan 10143, a few bad fits
    scan = 10143
    read(str(scan))
    bad = [132,220,255]
    updateRCP('laboca-10050.rcp')
    correctDelay (data,delay=1)
    medianbase()
    mediannoise()
    mediannoise()
    
    data.Data *= array(400./0.008,'f')
    op('beams_%i-b.ps/CPS'%(scan))
    iterat(chans=bad,scanNum=scan,oversamp=5,sizeX=[-250,250],sizeY=[-250,250],\
           data=data,flux=400.,radius=250.)
    close()
    data.updateArrayParameters("laboca-%s-b.rcp"%(scan))

def beam3():
    # older scan, but use newer rcp as starting point
    scan = 10051
    read(str(scan))
    
    updateRCP('laboca-10143.rcp')
    correctDelay (data,delay=5)
    medianbase()
    mediannoise()
    mediannoise()

    data.Data *= array(1500./0.03,'f')
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-150,150],\
           data=data,flux=1500.,radius=50.)
    close()

    data.updateArrayParameters("laboca-%i.rcp"%(scan))
    op('beams-shape_%i.ps/CPS'%(scan))
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[400,-400],sizeY=[-400,400],data=data)
    close()

def beam4():
    scan = 10203
    read(str(scan))
    
    updateRCP('laboca-10051.rcp')
    correctDelay (data,delay=1)
    medianbase()
    mediannoise()
    mediannoise()

    data.Data *= array(1500./0.03,'f')
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-150,150],\
           data=data,flux=1500.,radius=50.)
    close()

    data.updateArrayParameters("laboca-%i.rcp"%(scan))
    op('beams-shape_%i.ps/CPS'%(scan))
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[400,-400],sizeY=[-400,400],data=data)
    close()

def beam5():
    scan = 11207
    read(str(scan))
    updateRCP('master-laboca-may07.rcp')
    
    medianbase()
    mediannoise()
    data.Data *= array(1500./0.03,'f')
    cc = data.BolometerArray.checkChanList([]) # list of usable channels
    op('beams_%i.ps/CPS'%(scan))               # output plots, for ckecking the fit results
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-150,150],\
           data=data,flux=1500.,radius=50.)
    close()

    data.updateArrayParameters("laboca-%i.rcp"%(scan))
    op('beams-shape_%i.ps/CPS'%(scan))
    plotRCP(chans=cc,scanNum=scan,num=1,sizeX=[400,-400],sizeY=[-400,400],data=data)
    close()

def read5():
    f = file('beams_5_11102.txt')
    rl = f.readlines()
    f.close

    execfile('Test/cabling.py')
    bad = [85,86,96,136,173,191,192,201,205,220,202,218]
    bad.append(163)
    bad.append(4)
    bad.extend([101,103,106,110])
    posX = []
    posY = []
    wid  = []
    num  = []
    for l in rl:
        tmp = l.split()
        n = int(tmp[0])
        if n not in resistor and n not in bad:
            num.append(n)
            posX.append(float(tmp[1]))
            posY.append(float(tmp[2]))
            major = max([float(tmp[3]),float(tmp[4])])
            wid.append(major-22.)

    plot([0],[0],limitsX=[-400,350],limitsY=[-350,400],aspect=1)
    for i in range(len(num)):
        plot([posX[i]-wid[i],posX[i]+wid[i]],[posY[i],posY[i]],overplot=1,style='l')
        #Plot.xyout(posX[i],posY[i],str(num[i]))
                     
def beam6():
    #scan = 11540
    scan = 11924
    read(str(scan))
    updateRCP('master-laboca-may07.rcp')
    
    medianbase()
    mediannoise(nbloop=2)
    data.Data *= array(data.BolometerArray.BEGain / 270. * 7.15E6,'f') # 350 mV bias
    cc = data.BolometerArray.checkChanList([])
    op('beams_%i.ps/CPS'%(scan))
    iterat(chans=cc,scanNum=scan,oversamp=5,sizeX=[-150,150],sizeY=[-150,150],\
           data=data,flux=400.,radius=50.)
    close()

    data.updateArrayParameters("laboca-%i.rcp"%(scan))
    chanList = []
    bad = copy.copy(resistor)
    bad.extend([4,163,86,90,96,136,173,183,191,192,201,202,204,205])
    bad.extend([218,220,257,309])
    bad.extend([101,103,106,110,241])
    for i in range(len(data.arrayParamOffsets)):
        if data.arrayParamOffsets[i]['channel'] not in bad:
            chanList.append(data.arrayParamOffsets[i]['channel'])
    op('beams-shape_%i.ps/CPS'%(scan))
    plotRCP(chans=chanList,scanNum=scan,num=1,
            sizeX=[350,-400],sizeY=[-350,400],data=data)
    close()

def read6():
    global num,posX,posY,wid1,wid2
    #f = file('beams_5_11540.txt')
    f = file('beams_5_11924.txt')
    rl = f.readlines()
    f.close

    execfile('Test/cabling.py')
    bad = copy.copy(resistor)
    bad.extend([4,163,86,90,93,96,136,173,183,191,192,201,202,204,205])
    bad.extend([218,220,257,309])
    bad.extend([101,103,106,110,241])
    posX = []
    posY = []
    wid1  = []
    wid2  = []
    num  = []
    for l in rl:
        tmp = l.split()
        n = int(tmp[0])
        if n not in bad:
            num.append(n)
            posX.append(float(tmp[1]))
            posY.append(float(tmp[2]))
            wid1.append(max([float(tmp[3]),float(tmp[4])]))
            wid2.append(min([float(tmp[3]),float(tmp[4])]))

