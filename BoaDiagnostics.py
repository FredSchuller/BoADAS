#
# diagno.py
#
# Boa utilities with GUI
#
from Bogli.Interface import *
from Bogli import Plot

# external packages
from ppgplot import *
import string
from Numeric import array
import copy

def openGui(data):
    """ Open a new window to display graphic interface,
    and start waiting loop. Various diagnostic plots
    are displayed when the user clicks on buttons.

    INP data: object of class BoaDataAnalyser.DataAna on
              which processing is performed """

    global chanList,chanRef,integRms

    devId = pgopen('/XW')
    pgask(0)   # no prompt when closing window
    pgpap(10,0.3)             # window size (X in inch, ratio)
    pgscr(0, 0.85, 0.9, 0.85)  # background color
    pgsvp(0,1,0,1)           # defines world: 100% window
    pgenv(0,1000,300,0,0,-2) # world coordinate limits
    # window title
    pgsch(5)  # font size
    pgscf(2)  # font
    pgsci(4)  # font color
    pgptxt(500,0,0,0.5,"BoA Diagnostic Tools")

    pgsci(16)  # black
    # defines 7 buttons: times, tel. pattern, signal, FFT, rms, diff, correl
    center = []
    sizes = []
    x0 = 80
    dx = 140
    y0 = 80
    for i in range(7):
        center.append((x0+i*dx, y0))
        sizes.append((110,80))
    label = ["Times|MJD,LST","Telesc.|pattern","Signal","FFT",\
             "RMS","Diff.","Correl."]

    # Add two big boutons for saving / exit
    center.append((125,220))
    sizes.append((200,80))
    label.append('Save last')
    center.append((350,220))
    sizes.append((200,80))
    label.append('Close')

    # add indications
    # button for selection of channels
    xx=array([300,300,530,540,550,980,980])
    yy=array([130,145,145,205,145,145,130])
    pgsci(4)
    pgline(xx,yy)
    center.append((540,235))
    sizes.append((130,50))
    label.append('channels')
    # button for reference channel
    xx=array([720,720,840,850,860,980,980])
    yy=array([155,170,170,205,170,170,155])
    pgline(xx,yy)
    center.append((850,235))
    sizes.append((130,50))
    label.append('Ref chan')
    # RMS: number of integrations
    pgline(array([670,670]),array([125,160]))
    center.append((660,180))
    sizes.append((80,50))
    label.append('integ.')

    # default values
    pgsfs(1)
    pgsci(16)
    pgrect(475,605,270,300)
    pgrect(785,915,270,300)
    pgrect(620,700,215,245)
    pgsci(2)
    pgscf(1)
    pgsch(2)
    chanList = 'all'
    chanRef = data.BolometerArray.RefChannel
    integRms = 10
    pgptxt(485,290,0,0,str(chanList))
    pgptxt(795,290,0,0,str(chanRef))
    pgptxt(630,235,0,0,str(integRms))
    
    # display buttons and wait for clicks
    boutons = Fenetre(1,center,sizes,label,(0,10),\
                      font=3.,coltxt=1,colfond=14,family=2)
    boutons.dessine(new=0)

    again = 1
    while(again):
        numChoix = boutons.saisie()
        if numChoix == 8:
            again = 0
        else:
            processPlot(numChoix,data,devId)
        
    pgclos()


def processPlot(num,data,prevDev):
    """
    DES: do one of the diagnostic plots
    INP: (int) num = type of diagnostic
         (DataAna object) data = the data
         (int) prevdev = device number of main window
    """
    if num == 0:
        plotTimes(data)
    elif num == 1:
        plotPattern(data)
    elif num == 2:
        plotSignal(data)
    elif num == 4:
        plotRMS(data)
    elif num == 5:
        plotDiff(data)
    elif num == 6:
        plotCorrel(data)

    elif num == 9:
        # select channels
        selectChan(data)
    elif num == 10:
        # select ref. channel
        selectRef(data)
    elif num == 11:
        # set nb. of integrations to compute rms
        setInteg(data)
    
    else:
        print num,'not implemented'

    pgslct(prevDev)

##############
# individual plotting routines
#
def plotTimes(data):

    pgopen('/XW')
    pgask(0)

    lst=data.ScanParam.LST
    MJD_renorm=(data.ScanParam.MJD-data.ScanParam.MJD[0])*86400.

    Plot.plot(lst,style='l',labelX='Integration number',\
              labelY='LST [sec]',\
              caption='Scan '+str(data.ScanParam.ScanNum)+\
              ' - '+data.ScanParam.Object + \
              ' - '+data.BolometerArray.FeBe +\
              ' - '+data.BolometerArray.Telescope.Name)
    print "Press return for next plot"
    raw_input()

    Plot.plot(MJD_renorm,style='l',labelX='Integration number',\
              labelY='MJD - MJD[0] [sec]',\
              caption='Scan '+str(data.ScanParam.ScanNum)+\
              ' - '+data.ScanParam.Object + \
              ' - '+data.BolometerArray.FeBe +\
              ' - '+data.BolometerArray.Telescope.Name)

    dt=[]
    nbInteg = len(lst)
    for i in range(nbInteg-1):
        dt.append(MJD_renorm[i+1]-MJD_renorm[i])

    print "Press return for next plot"
    raw_input()
    Plot.plot(1./array(dt),style='p',labelX='Integration number',\
              labelY='Sampling rate [Hz]',\
              caption='Scan '+str(data.ScanParam.ScanNum)+\
              ' - '+data.ScanParam.Object + \
              ' - '+data.BolometerArray.FeBe +\
              ' - '+data.BolometerArray.Telescope.Name)

    print "Press return for next plot"
    raw_input()

    Plot.plot(lst,MJD_renorm,style='p',labelX='LST [sec]',\
              labelY='MJD - MJD[0] [sec]',\
              caption='Scan '+str(data.ScanParam.ScanNum)+\
              ' - '+data.ScanParam.Object + \
              ' - '+data.BolometerArray.FeBe +\
              ' - '+data.BolometerArray.Telescope.Name)
    print "Press return to close window"
    raw_input()
    pgclos()

def plotPattern(data):

    pgopen('/XW')
    pgask(0)

    data.ScanParam.plotAzEl()
    data.ScanParam.plotAzEl(style='p',overplot=1,ci=2)

    print "Press return for next plot"
    raw_input()

    data.ScanParam.plotAzElOffset()
    data.ScanParam.plotAzElOffset(style='p',overplot=1,ci=2)
    print "Press return to close window"
    raw_input()
    pgclos()

def plotSignal(data):

    # simply use current list of channels
    pgopen('/XW')
    pgask(0)

    data.signal()
    print "Press return for next plot"
    raw_input()
    data.signal(plotMap=1)
    print "Press return to close window"
    raw_input()
    pgclos()
    
def plotFFT(data):
    # to be done...
    #data.fft([17],mjd=1,style='l')
    pass

def plotRMS(data):

    global integRms
    slideRms = data.slidingRms(nbInteg = integRms)
    pgopen('/XW')
    pgask(0)
    Plot.draw(slideRms,wedge=1,caption='Signal RMS',labelY = 'Channel #')
    print "Press return to close window"
    raw_input()
    pgclos()
    
def plotDiff(data):

    global chanRef

    backup = copy.copy(data.Data)
    # subtract signal of chanRef
    for n in data.BolometerArray.CurrChanList:
        data.Data[:,n-1] = data.Data[:,n-1]-backup[:,chanRef-1]

    pgopen('/XW')
    pgask(0)
    data.signal()
    print "Press return for next plot"
    raw_input()
    data.signal(plotMap=1)
    print "Press return to close window"
    raw_input()
    pgclos()

    # restore data to initial state
    data.Data = copy.copy(backup)

    
def plotCorrel(data):

    global chanRef

    pgopen('/XW')
    pgask(0)
    # simply use current list of channels
    data.plotCorrel(chanRef = chanRef)
    print "Press return to close window"
    raw_input()
    pgclos()
    

##############
# routines to set paramaters
#
def selectChan(data):
    """ open a box to fill by the user with valid
    list of channels. Return list when valid

    INP: (object) data = data to process
    """

    global chanList
    
    answer = fenetreInteractive(x0=475,y0=270,x1=605,y1=300,\
                                      prevText=chanList,bgCol=16)
    chan = answer
    # interpret answer like '1-20'
    if '-' in chan:
        words = string.split(chan,'-')
        try:
            toList = range(int(words[0]),int(words[1])+1)
            chan = toList
        except ValueError:  # not integers
            data.MessHand.warning("Could not interpret input list")
    elif chan != 'all':  # one integer, or coma-separated list of int
        num = string.split(chan,',')
        toList = []
        try:
            for n in num:
                toList.append(int(n))
            chan = toList
        except ValueError:  # not integers
            data.MessHand.warning("Could not interpret input list")
        
    # message handler will write error message if not
    # a valid list
    data.BolometerArray.setCurrChanList(chan)
    chanList = answer

def selectRef(data):
    """ open a box to input ref. channel number
    INP: (object) data = data to process
    """

    global chanRef

    again = 1
    while(again):
        answer = fenetreInteractive(x0=785,y0=270,x1=915,y1=300,\
                                    prevText=str(chanRef),bgCol=16)
        try:
            newRef = int(answer)
            if newRef in data.BolometerArray.UsedChannels:
                chanRef = newRef
                again = 0
            else:
                data.MessHand.warning("Channel not in use")
        except ValueError:
            # not an integer
            data.MessHand.warning("Could not interpret input value")

def setInteg(data):
    """ open a box to input number of integrations for computing rms
    INP: (object) data = data to process
    """

    global integRms

    again = 1
    while(again):
        answer = fenetreInteractive(x0=620,y0=215,x1=700,y1=245,\
                                    prevText=str(integRms),bgCol=16)
        try:
            newInteg = int(answer)
            if newInteg > 1: # need at least 2 integ.
                integRms = newInteg
                again = 0
            else:
                # should use message handler
                data.MessHand.warning("Need at least two integrations to compute rms")
        except ValueError:
            # not an integer
            data.MessHand.warning("Could not interpret input value")
