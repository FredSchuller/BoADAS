#!/usr/bin/env python

# APEX - Atacama Pathfinder EXperiment Project
#
# who       when        what
# --------  ----------  ----------------------------------------------
# mdumke    2010-05-01  Added flagging & inverting as f(MJD)
# alundgre  2010-04-11  Sign switch commented out
# aweiss    2010-03-14  Created
#
'''
SABOCA On-Off Reduction
'''
import random

saboca = os.getenv('BOA_HOME_SABOCA')
execfile(saboca+'/saboca-cabling.py')
execfile(saboca+'/boaSaboca.py')

def GetOnBolometers(tol=1.0):
    
    tolerance = data.BolometerArray.BeamSize/tol

    #GET AZ,EL OFFSETS OF REFERENCE PIXEL
    refchan = data.BolometerArray.RefChannel
    refaz = -1*data.BolometerArray.Offsets[0,refchan-1]
    refel = -1*data.BolometerArray.Offsets[1,refchan-1]

    #CALCULATE INDIX OF POSITION CHANGES
    
    azoff=data.ScanParam.get('azoff')
    eloff =data.ScanParam.get('eloff')
    SubIndex = []
    oldpos=[azoff[0],eloff[0]]

    for i in range(len(azoff)):
        if i == 0:
            SubIndex.append(0)
        else:
            newpos=[azoff[i],eloff[i]]
            dist=sqrt((oldpos[0]-newpos[0])**2+(oldpos[1]-newpos[1])**2)
            if (dist > tolerance):
                SubIndex.append(i)
                oldpos=copy.deepcopy(newpos)
                
                

    SubIndex.append(len(azoff))

    #DETERMINE ON BOLOMETERS BASED ON AZ EL OFFSETS FOR EACH INTERVAL

    onbolos = []
    index = []
    
    for i in range(len(SubIndex)-1):
        medaz=-1*fStat.f_median(azoff[SubIndex[i]:SubIndex[i+1]])-refaz
        medel=-1*fStat.f_median(eloff[SubIndex[i]:SubIndex[i+1]])-refel
        found = 0
        for j in range(data.BolometerArray.NChannels):
            dist = sqrt((medaz-data.BolometerArray.Offsets[0,j])*(medaz-data.BolometerArray.Offsets[0,j])+(medel-data.BolometerArray.Offsets[1,j])*(medel-data.BolometerArray.Offsets[1,j]))
            if dist < tolerance and found == 0:
                if i == 0:
                    onbolos.append(j+1)
                    lastbolo = j+1
                    index.append(SubIndex[i])
                    found = 1
                else:
                    test = j+1
                    if test != lastbolo:
                        onbolos.append(j+1)
                        index.append(SubIndex[i])
                        found = 1
                        lastbolo = j+1
                    
                
                
    index.append(len(azoff))
    return onbolos,index

####################################

def flagPhaseOn():
	chref = data.BolometerArray.RefChannel
	flagPosition(channel=[chref],radius=20,flag=8)
	mask = data.FlagHandler.isSetMask([8],dim=1,index=chref-1)
	chanList = data.BolometerArray.checkChanList([])
	chanListIndexes = data.BolometerArray.getChanIndex(chanList)
	for chan in chanListIndexes:
        	data.FlagHandler.setOnMask(mask,8,dim=1, index=chan)

def flagPhaseOff():
        chref = data.BolometerArray.RefChannel
        flagPosition(channel=[chref],radius=20,flag=8,outer=1)
        mask = data.FlagHandler.isSetMask([8],dim=1,index=chref-1)
        chanList = data.BolometerArray.checkChanList([])
        chanListIndexes = data.BolometerArray.getChanIndex(chanList)
        for chan in chanListIndexes:
                data.FlagHandler.setOnMask(mask,8,dim=1, index=chan)

def flagOnePhase(side=0):
    # flags all bolometers when the wobbler is on one side
    # INP: (int) side, 0 for left, 1 for right
    #
    # this has to be done per subscan
    nbSub = data.ScanParam.NObs
    mask = array([],'i')
    for i in range(nbSub):
        subNum = data.ScanParam.SubscanNum[i]
        # get all Az offsets for that subscan, and compute mean
        azo = data.ScanParam.get('azoff',subscans=[subNum])
        meanAz = fStat.f_mean(azo)
        # now determine what is left, what is right
        # But we need to work on all datapoints, including flagged ones
        azo = data.ScanParam.get('azoff',subscans=[subNum],flag='None')
        maskR = greater(azo,meanAz)
        maskL = less(azo,meanAz)
        if side:
            mask = concatenate((mask,maskR))
        else:
            mask = concatenate((mask,maskL))

    chanList = data.BolometerArray.checkChanList([])
    chanListIndexes = data.BolometerArray.getChanIndex(chanList)
    for chan in chanListIndexes:
        data.FlagHandler.setOnMask(mask,8,dim=1, index=chan)
        
def flagPhaseLeft():
    # flags all bolometers when the wobbler is on left position
    flagOnePhase(side=0)
    
def flagPhaseRight():
    # flags all bolometers when the wobbler is on right position
    flagOnePhase(side=1)
    

def calibrate(data,tau):
    mjdref = data.ScanParam.MJD[0]
#     ###CntstoJy(data)
    
    rcp=getSabocaRCP(mjdref)
    print rcp
    data.BolometerArray.updateRCP(rcp)
    
    flagRCP(rcp)
    cross=getSabocaCross(mjdref)
    print cross
    #data.flagChannels(cross)
#     data.flatfield()
    data.correctOpacity(tau)
    scanel  = fStat.f_mean(data.ScanParam.get('El'))
    taucorr = exp(tau/sin(scanel * pi / 180.))
    return taucorr

def flagCable(data,cable='M'):
    # flag channels corresponding to one cable
    filename = os.path.join(os.getenv('BOA_HOME_RCP'),
                            'channels_full_set.20070509.txt')
    data.BolometerArray.readAdditionnalIndexFile(filename,0,5,'#')
    chanList = data.BolometerArray.selectAdditionnalIndex(cable)
    data.flagChannels(chanList)

def process(data,weak=0):
    # Standard reduction of data, before computing on-off

    data.zeroStart()
    data.Data *= array(VtoJy,'f')

    mjdref = data.ScanParam.MJD[0]

    # RCP:
    rcp=getSabocaRCP(mjdref)
    data.MessHand.info("RCP file : %s " % rcp)
    updateRCP(rcp)
    # Bad channels:
    bad=getSabocaCross(mjdref)
    data.MessHand.info("Bad channels : %s " % bad)
    data.flagChannels(bad)
    # Inverted channels:
    inverted=getSabocaInvert(mjdref)
    data.MessHand.info("Inverted channels : %s " % inverted)
    invertSomeChannels(inverted)

    flagRCP(rcp)
    flat()
    # flag all data taken while wobbler is moving
    refch=data.BolometerArray.RefChannel
    data.flagSpeed(above=30.)
    # zero-order baseline per subscan
    data.polynomialBaseline(order=0,subscan=1)

    if weak == 1:
        data.flagFractionRms(ratio=5)
        data.despike()
        data.medianNoiseRemoval(computeFF=1,chanRef=-1,nbloop=1,factor=0.2)
        data.despike()
        data.flagFractionRms(ratio=3)
        data.medianNoiseRemoval(computeFF=1,chanRef=-1,nbloop=3,factor=0.7)
        data.flagFractionRms(ratio=3)
        data.despike()
        data.polynomialBaseline(order=1,subscan=1)
        data.flattenFreq(below=0.5/data.ScanParam.WobCycle,
                         hiref=0.7/data.ScanParam.WobCycle)
        data.polynomialBaseline(order=0,subscan=1)
        data.despike()
    else:
        data.flagPosition(radius=20)
        data.flagFractionRms(ratio=5)
        data.despike()
        data.medianNoiseRemoval(computeFF=1,chanRef=-1,nbloop=1,factor=0.2)
        data.despike()
        data.flagFractionRms(ratio=3)
        data.medianNoiseRemoval(computeFF=1,chanRef=-1,nbloop=3,factor=0.7)
        data.despike()
        data.flagFractionRms(ratio=3)
        data.polynomialBaseline(order=1,subscan=1)
        data.unflag(flag=8)

    data.unflagChannels([refch])
    
def reduceOO(ScanNr=0,tau=1.0,subscans=[],weak=0):
    tst = read(ScanNr,subscans=subscans)
    if tst:
        return -999,-999,-999,-999,-999
    
    # define result arrays: phase_dif + dispersion for each subscan
    result,sigma,inttime = [],[],[]
    count = 0
    refch=data.BolometerArray.RefChannel
    rmson = []
    rmsoff = []
    rmsdiff = []
    
    # calibrate the data
    taucorr = calibrate(data,tau)
    scantime = (data.ScanParam.MJD[-1]-data.ScanParam.MJD[0])*24*3600
    process(data,weak=weak)
    data.signal(refch)

    # Store list of "off" bolo
    offChanList = data.BolometerArray.checkChanList([])  # all valid channels
    offChanList = offChanList.tolist()
    offChanList.remove(refch)
    nbChan = len(offChanList) + 1  # including the "on" bolo

    # build arrays with flux, time and flags for the full scan
    flux = data.getChanData('flux',refch,flag='None')
    offFluxes = data.getChanListData('flux',offChanList,dataFlag='None')
    allFluxes = [flux]
    allFluxes.extend(offFluxes)
    allFluxes = array(allFluxes)
    tt = data.getChanData('mjd',refch,flag='None')

    # compute arrays of flags corresponding to each phase
    flagPhaseOn()
    flagOff = copy.copy(data.getChanData('flag',refch,flag='None'))
    data.signal(refch,overplot=1,ci=2,style='p')
    data.unflag(flag=8)
    
    flagPhaseOff()
    flagOn = copy.copy(data.getChanData('flag',refch,flag='None'))
    data.signal(refch,overplot=1,ci=3,style='p')
    data.unflag(flag=8)
    
    # now starts loop on subscans
    nbSub = len(data.ScanParam.SubscanNum)
    subInd = data.ScanParam.SubscanIndex
    subNum = data.ScanParam.SubscanNum
    for sub in range(nbSub):
        # is it nodding A or B?
        subtype = data.ScanParam.WobblerSta[subInd[0,sub]]
        #            print "Subscan: %i  - type: %s"%(subNum[sub],str(subtype))
            
        if subNum[sub]%2 == 1:   # if subscans is odd, reset list of phase_diff and errors:
            tmpRes,dtmpRes  = [[]],[[]] # to store phase_dif and error for this subscan
            for c in range(nbChan-1):
                tmpRes.append([])
                dtmpRes.append([])
            tmptime = []
            count = count+1
            
        # here starts the loop that computes each on-off
        nMax = subInd[1,sub]
        index = subInd[0,sub]

        while(index < nMax-1):
            timeOn  = []
            timeOff = []
            while(flagOn[index] and index < nMax-1):
                index += 1
            # here is the 1st value for "on"
            onStart = index
            while(index < nMax and not flagOn[index]):
                timeOn.append(tt[index])
                index += 1
            # that was the end of this "on" serie
            onEnd = index
            
            if index < nMax:
                # not yet at the end, go to beginning of next "off"
                while(flagOff[index] and index < nMax-1):
                    index += 1
                # here is the 1st value for "off"
                offStart = index
                while(index < nMax and not flagOff[index]):
                    timeOff.append(tt[index])
                    index += 1
                # that was the end of this "off" serie
                offEnd = index

            # if at least two values in each phase, compute mean values
            if len(timeOn) > 1 and len(timeOff) > 1:
                valOn = allFluxes[:,onStart:onEnd]
                valOff = allFluxes[:,offStart:offEnd]
                intOn = timeOn[-1]-timeOn[0]
                tOn = fStat.f_mean(timeOn)
                tOff = fStat.f_mean(timeOff)
                for num in range(nbChan):
                    onChan  = fStat.f_median(valOn[num])
                    offChan = fStat.f_median(valOff[num])
                    rmsOn   = fStat.f_rms(valOn[num],onChan)
                    rmsOff  = fStat.f_rms(valOff[num],offChan)
                    # compute and store phase diff + error for all channels
                    theDiff = onChan - offChan
                    dtheDiff = sqrt(rmsOn**2+rmsOff**2)
                    tmpRes[num].append(theDiff)
                    dtmpRes[num].append(dtheDiff/sqrt(len(valOn[num])))
                    # and only for ref. chan also stores rms on each phase
                    if num == 0:
                        sigOn  = onChan
                        sigOff = offChan
                        rmson.append(rmsOn)
                        rmsoff.append(rmsOff)
                        rmsdiff.append(dtheDiff)
                       
                tmptime.append(intOn)

                # plot all results
                plot([tOn],[sigOn],overplot=1,ci=2)
                plot([tOn,tOn],[sigOn-rmson[-1],sigOn+rmson[-1]],
                     overplot=1,ci=2,style='l')
                plot([tOff],[sigOff],overplot=1,ci=3)
                plot([tOff,tOff],[sigOff-rmsoff[-1],sigOff+rmsoff[-1]],
                     overplot=1,ci=3,style='l')
                
        if  subNum[sub]%2 == 0:
            meanRes = []
            sigRes  = []
            sumtime = sum(tmptime)
            for num in range(nbChan):
                w = 1./array(dtmpRes[num])**2
                wsum = sum(w)
                meanRes.append(sum(array(tmpRes[num])*w) / wsum)
                sigRes.append(sqrt(1./wsum))
                           
            result.append(meanRes)
            sigma.append(sigRes)
            inttime.append(sumtime)
            print "Nod-cycle: %i :"%(count)
            print " ... %6.3f +/- %5.3f Jy/b"%(meanRes[0],sigRes[0])
            
    ontime = sum(inttime)
    scaneff = ontime/scantime
    print "Scan time efficiency: %6.1f %% "%(100.*scaneff)
    source = data.ScanParam.Object
    mrmson   = fStat.f_mean(rmson)
    mrmsoff  = fStat.f_mean(rmsoff)
    mrmsdiff = fStat.f_mean(rmsdiff)

    # compute NEFD for each phase, on the full scan
    adt=(data.ScanParam.get('deltat'))
    dt=fStat.f_median(adt)
    print 1/dt
    sendiff=mrmsdiff/taucorr*sqrt(dt)*1000.
    senon=mrmson/taucorr*sqrt(dt)*1000.
    senoff=mrmsoff/taucorr*sqrt(dt)*1000.
    print "Scan NEFD: On: %6.1f, Off: %6.1f mJy sqrt(s), "%(senon,senoff)
    print "          phase difference: %6.1f mJy sqrt(s)"%(sendiff)

    return result,sigma,inttime,taucorr,source
    
#####################################################################

def doOO(scanlist=[],tau=1.0,subscans=[],weak=1,
         filename='',update=0,source=''):
    # to do: use a list of taus as input

    # Initialise or lists and dictionary to store results
    fluxlist  = []
    dfluxlist = []
    nefdlist  = []
    timelist  = []
    rmsOff,meanOff = [],[]

    # reading previous results if required
    if update:
        if os.path.exists(filename):
            fin = file(filename)
            rl = fin.readlines()
            fin.close()
            for l in rl:
                tmp = l.split()
                fluxlist.append(float(tmp[0]))
                dfluxlist.append(float(tmp[1]))
                nefdlist.append(float(tmp[2]))
                timelist.append(float(tmp[3]))
                meanOff.append(float(tmp[4]))
                rmsOff.append(float(tmp[5]))
        
    fluxOff   = {}
    dfluxOff  = {}
    cumOff    = {}
    dcumOff   = {}
    src = ''

    for snum in scanlist:
        flux,dflux,tOn,taucorr,src = reduceOO(snum,tau,subscans,weak=weak)
        if flux != -999:
            nb = len(flux)
            nbChan = len(flux[0])
            timelist.extend(tOn)
            for i in range(nb):
                fluxlist.append(flux[i][0])
                dfluxlist.append(dflux[i][0])
                nefdlist.append(1000.*dflux[0][i]/taucorr*sqrt(tOn[i]))
                # compute weighted mean and rms of OFF signals
                offSig = flux[i][1:]
                nbOff = float(len(offSig))
                offWei = 1./array(dflux[i][1:])**2
                w_offSig = offWei * array(offSig)
                meanOff.append(sum(w_offSig) / sum(offWei))
                rmsOff.append(sum(1./array(dflux[i][1:])) / sum(offWei))
                
    # Export results, if required
    entries = len(fluxlist)
    if filename and scanlist:
        # don't write file if no input scan
        fout = file(filename,'w')
        for i in range(entries):
            oneLine = str("%10.6f %10.6f "%(fluxlist[i],dfluxlist[i]))
            oneLine += str("%6.2f %6.2f "%(nefdlist[i],timelist[i]))
            oneLine += str("%10.6f %10.6f\n"%(meanOff[i],rmsOff[i]))
            fout.write(oneLine)
        fout.close()

    # Compute NEFD and statistical errors
    meannefd = fStat.f_mean(nefdlist)
    mean = fStat.f_mean(fluxlist)
    sysrms=fStat.f_rms(fluxlist,mean)
    syserror=sysrms/sqrt(len(fluxlist))

    # Compute cumulated fluxes (ON and OFF), and their errors
    cumflux = zeros(entries,'f')
    dcumflux= zeros(entries,'f')
    cumoff  = zeros(entries,'f')
    dcumoff = zeros(entries,'f')
    cumtime = zeros(entries,'f')
    inlist  = range(entries)

    # Shuffle all measurements many times...
    nbIter = 100
    for k in range(nbIter):
        outlist = reorder(inlist)
        sumtime = 0.0
        sumflux = 0.0
        sumoff  = 0.0
        wonsum  = 0.0
        woffsum = 0.0
        for i in range(entries):
            num = outlist[i]
            won    = 1./dfluxlist[num]**2
            wonsum = wonsum+won
            woff   = 1./rmsOff[num]**2
            woffsum= woffsum+woff
            sumtime = sumtime+timelist[num]
            sumflux = sumflux+fluxlist[num]*won
            sumoff  = sumoff+meanOff[num]*woff
            cumflux[i] = cumflux[i]+sumflux/wonsum
            dcumflux[i]= dcumflux[i]+sqrt(1./wonsum)
            cumoff[i]  = cumoff[i]+sumoff/woffsum
            dcumoff[i] = dcumoff[i]+sqrt(1./woffsum)
            cumtime[i] = cumtime[i]+sumtime
            
    cumflux  = cumflux/float(nbIter)
    dcumflux = dcumflux/float(nbIter)
    cumoff   = cumoff/float(nbIter)
    dcumoff  = dcumoff/float(nbIter)
    cumtime  = cumtime/float(nbIter)
    
    # Finally, plot the results
    toplot = cumflux.tolist()
    toplot.extend(cumflux+dcumflux)
    toplot.extend(cumflux-dcumflux)
    toplot.extend(cumoff+dcumoff)
    toplot.extend(cumoff-dcumoff)
    mini,maxi = fStat.minmax(toplot)
    xmin,xmax=fStat.minmax(cumtime)

    Plot.panels(1,2)
    BogliConfig.box['charheight'] = 1.5
    BogliConfig.cLabel['charheight'] = 1.5
    BogliConfig.xAxis['charheight'] = 1.5
    BogliConfig.yAxis['charheight'] = 1.5
    BogliConfig.globalViewPoint['marging_top'] = 0.9
    BogliConfig.globalViewPoint['marging_bottom'] = 0.15

    # First plot: down-integrating results
    Plot.nextpage()
    if not source:
        source = src
    if maxi > 1.0:
        plot(cumtime,cumflux,limitsX=[xmin,xmax],limitsY=[mini,maxi],
             style='l',ci=2,labelX='ON-Time [s]',labelY='Flux [Jy]',
             caption='Cumulative on-off flux on '+source,noerase=1)
        for i in range(len(cumflux)):
            up = (cumflux[i]+dcumflux[i])
            low = (cumflux[i]-dcumflux[i])
            plot([cumtime[i],cumtime[i]],[low,up],overplot=1,ci=2,style='l')
        plot(cumtime,cumoff+dcumoff,style='l',ci=3,overplot=1)
        plot(cumtime,cumoff-dcumoff,style='l',ci=3,overplot=1)
        plot(cumtime,cumoff,style='l',ci=1,overplot=1)
    else:
        plot(cumtime,1000.0*array(cumflux),limitsX=[xmin,xmax],
             limitsY=[1000.0*mini,1000.0*maxi],style='l',ci=2,
             labelX='ON-Time [s]',labelY='Flux [mJy]',
             caption='Cumulative on-off flux on '+source,noerase=1)
        for i in range(len(cumflux)):
            up = 1000.0*(cumflux[i]+dcumflux[i])
            low = 1000.0*(cumflux[i]-dcumflux[i])
            plot([cumtime[i],cumtime[i]],[low,up],overplot=1,ci=2,style='l')
        plot(cumtime,1000.*cumoff+1000.*dcumoff,style='l',ci=3,overplot=1)
        plot(cumtime,1000.*cumoff-1000.*dcumoff,style='l',ci=3,overplot=1)
        plot(cumtime,1000.*cumoff,style='l',ci=1,overplot=1)

    # and print the final numbers
    print "###################################################"
    if cumflux[-1] > 1.:
        print "Final result = %6.3f +/- %6.3f [+/- %6.3f] Jy/b"%(cumflux[-1],dcumflux[-1],syserror)
    else:
        print "Final result = %6.2f +/- %6.2f [+/- %6.2f] mJy/b"%(1000.*cumflux[-1],1000.*dcumflux[-1],1000.*syserror)
        
    print "Average OFF bolo = %6.2f +/- %5.1f mJy/b"%(1000.*cumoff[-1],1000.*dcumoff[-1])
    print "Effective NEFD: %6.1f mJy sqrt(s)"%(meannefd)
    print "###################################################" 

    # Second plot: show individual data points with error bars
    Plot.nextpage()
    valmin = array(fluxlist) - array(dfluxlist)
    valmax = array(fluxlist) + array(dfluxlist)
    sumtime = array(timelist)
    for i in range(entries):
        sumtime[i] = sum(timelist[:i+1])
    ymin,ymax = fStat.minmax(concatenate([valmin,valmax]))
    pointSize(6)
    if ymax > 1.0:
        plot(sumtime,fluxlist,noerase=1,labelX='ON-Time [s]',
             labelY='Flux [Jy]',limitsY=[ymin,ymax],
             caption='Individual ON-OFF measurements',ci=2)
        for i in range(entries):
            plot([sumtime[i],sumtime[i]],[valmin[i],valmax[i]],
                 overplot=1,style='l',ci=2)
    else:
        plot(sumtime,1000.*array(fluxlist),noerase=1,labelX='ON-Time [s]',
             labelY='Flux [mJy]',limitsY=[1000.*ymin,1000.*ymax],
             caption='Individual ON-OFF measurements',ci=2)
        valmin *= 1000.
        valmax *= 1000.
        for i in range(entries):
            plot([sumtime[i],sumtime[i]],[valmin[i],valmax[i]],
                 overplot=1,style='l',ci=2)

    #return sumtime,fluxlist

#-------------------------------------------------------------

def reorder(inlist):
    nb = len(inlist)
    for i in range(nb):
        j = random.randrange(i,nb)
        inlist[i],inlist[j] = inlist[j],inlist[i]

    for i in range(nb):
        j = random.randrange(i,nb)
        inlist[i],inlist[j] = inlist[j],inlist[i]

    for i in range(nb):
        j = random.randrange(i,nb)
        inlist[i],inlist[j] = inlist[j],inlist[i]

    return inlist




def flagElJitter(scan,deltaEl=data.BolometerArray.BeamSize/5.0):
    read(scan)
    data.flagSpeed(above=30.)
    flagPosition(radius=20,outer=1,flag=8)
    ondumps=len(data.getChanData('eloff'))
    melof=fStat.f_mean(data.getChanData('eloff'))
    elrms=fStat.f_rms(data.getChanData('eloff'),melof)
    low = melof-deltaEl
    up=melof+deltaEl
    data.flag(dataType='eloff',channel=data.BolometerArray.RefChannel,below=low,flag=3)
    data.flag(dataType='eloff',channel=data.BolometerArray.RefChannel,above=up,flag=3)
    remaindumps=len(data.getChanData('eloff'))
    elrms=fStat.f_rms(data.getChanData('eloff'),melof)
    fracflag=100.0*(1.0-1.0*remaindumps/ondumps)
    unflag(flag=8)
    print "%4.1f %% of on dumps flagged"%(fracflag)
    print "remaining rms in elevation tracing: %4.1f arcsec"%(elrms)
