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
NAM: BoaDir (module)
DES: contains all the method to handle list of MBFits file
"""
__version__=  '$Revision: 2794 $'
__date__=     '$Date: 2015-03-02 15:30:03 +0100 (Mon, 02 Mar 2015) $'

import BoaConfig
import os, sys, re, string
import BoaMessageHandler
import BoaMBFits, BoaMBFitsReader

MessHand=BoaMessageHandler.MessHand('BoaDir')

ListDir = []
CurrentList = []

# -------------------------------------------------------------------
def setInDir(inputDirectory=''): 
    """
    DES: set the input directory
    INP: (string) inputDirectory = name of input directory
    """

    global ListDir, CurrentList

    if not isinstance(inputDirectory,str):
        MessHand.error("invalid input directory"+str(inputDirectory))
        return

    if inputDirectory == '':
        MessHand.info("input directory : "+BoaConfig.inDir)
        return

    inputDirectory = inputDirectory.strip()

    if inputDirectory[-1] != '/':
        inputDirectory = inputDirectory+'/'

    if os.path.isdir(inputDirectory)==1:
        BoaConfig.inDir=inputDirectory
        MessHand.info("input directory : "+BoaConfig.inDir)
        ListDir = []
        CurrentList = []
    else: 
        MessHand.error("no such directory: "+inputDirectory)

# -------------------------------------------------------------------
def setProjectID(projectID=None):
    """
    DES: set the current project ID
    """

    if not projectID :
        MessHand.info("project ID : "+BoaConfig.projID)
        return

    # the format of the project ID changed several times
    # so, the following regular expression keeps getting more
    # complicated
    checkID = re.compile('^\w-\d{2,3}\.\w-((\d{4,4})|(\d{4,4}\w))-\d{4,4}$')

    if not checkID.match(projectID):
        MessHand.error("not a valid APEX project ID")
        return

    BoaConfig.projID = projectID

# -------------------------------------------------------------------
def setDate(date=None):
    """
    DES: set the current observing date
    """

    if not date :
        MessHand.info("obs. date : "+BoaConfig.date)
        return

    checkDate = re.compile('\d{4,4}-\d{2,2}-\d{2,2}')

    if not checkDate.match(date):
        MessHand.error("not a valid observing date")
        return

    BoaConfig.date = date


# -------------------------------------------------------------------
def setOutDir(outputDirectory=''): 
    """
    DES: set the output directory
    INP: (string) outputDirectory = name of output directory
    """
    if not isinstance(outputDirectory,str):
        MessHand.error("invalid output directory"+str(outputDirectory))
        return

    if outputDirectory == '':
        MessHand.info("output directory : "+BoaConfig.outDir)
        return

    outputDirectory = outputDirectory.strip()

    if outputDirectory[-1] != '/':
        outputDirectory = outputDirectory+'/'

    if os.path.isdir(outputDirectory)==1:
        BoaConfig.outDir=outputDirectory
        MessHand.info("output directory : "+BoaConfig.outDir)
    else: 
        MessHand.error("no such directory: "+outputDirectory)


# -------------------------------------------------------------------
def setInFile(inputFile='?'): 
    """
    DES: set the input file name
    INP: (string) inputFile = name of input file 
    """

    if not isinstance(inputFile,str):
        MessHand.error("invalid input file name"+str(inputFile))
        return

    if inputFile == '?':
        MessHand.info("input file : "+BoaConfig.inFile)
        return

    inputFile = inputFile.strip()

    if inputFile[-5::] == '.fits':
        inputFile = inputFile[:-5]
    
    if inputFile[-8::] == '.fits.gz':
        inputFile = inputFile[:-8]
    
    BoaConfig.inFile=''

    testfiles=[inputFile,\
               inputFile+'.fits',\
               inputFile+'.fits.gz',\
               'APEX-'+inputFile+'-'+BoaConfig.projID,\
               'APEX-'+inputFile+'-'+BoaConfig.projID+'.fits',\
               'APEX-'+inputFile+'-'+BoaConfig.projID+'.fits.gz',
               'APEX-'+inputFile+'-'+BoaConfig.date+'-'+BoaConfig.projID]
    
    for testname in testfiles:
        if os.path.isfile(BoaConfig.inDir+testname)==1: 
            BoaConfig.inFile=BoaConfig.inDir+testname
            MessHand.debug("input file = "+ testname)
            return
        elif os.path.isdir(BoaConfig.inDir+testname)==1: 
            BoaConfig.inFile=BoaConfig.inDir+testname
            MessHand.debug("input directory = "+ testname)
            return

    
    if BoaConfig.inFile == '':
        listing = os.listdir(BoaConfig.inDir)
        for filename in listing:
            if string.find(filename,inputFile) >= 0:
                BoaConfig.inFile = BoaConfig.inDir + filename
                MessHand.debug("input file = "+ testname)
                return

        MessHand.error("no such file or directory: "+ BoaConfig.inDir+str(inputFile))

# -------------------------------------------------------------------
def setOutFile(outputFile='?'): 
    """
    DES: set the output file name
    INP: (string) outputFile = name of output file 
    """

    if not isinstance(outputFile,str):
        MessHand.error("invalid output file name"+str(outputFile))
        return

    if outputFile == '?':
        MessHand.info("output file : "+BoaConfig.outFile)
        return

    outputFile = outputFile.strip()

    if outputFile[-5::] != '.fits':
        outputFile = outputFile+'.fits'

    if os.path.isfile(BoaConfig.inDir+outputFile)==1: 
        BoaConfig.inFile=BoaConfig.inDir+outputFile
        MessHand.info("output file = "+ BoaConfig.outFile)
    else:
        MessHand.error("no such file or directory: "+ BoaConfig.inDir+str(outputFile))


# -------------------------------------------------------------------
def listInDir(separator=' '):
    """
    DES: list the input directory
    INP: (str) field separator
         e.g. separator = '||' for moinmoin-style table
    """

    global CurrentList
    if CurrentList in [[]]:
        MessHand.warning('list empty, populating...')
        findInDir()

    toPrint        = ['Object','ScanType','NSubscans','ExpTime','FEBEList',\
                      'refChan','MJD' ,'Size']
    toPrintFormat  = ['%15s'  ,'%10s'    ,'%3i'      ,'%7.0f'  ,'%17s'    ,\
                      '%5s'    ,'%9.3f','(%5.1f MB)']
    
    for fitsfile in CurrentList:
        if fitsfile['status']:
            
            # Reformat the filename for display
            displayName = fitsfile['filename']
            if displayName[-5::] == '.fits':
                displayName = displayName[:-5]
            if displayName[-8::] == '.fits.gz':
                displayName = displayName[:-8]
            if displayName[:5] == 'APEX-' and displayName[-len(BoaConfig.projID)::] == BoaConfig.projID:
                displayName = displayName[5:-len(BoaConfig.projID)-1]

            outString = separator + "%7s " % displayName

            for i in range(len(toPrint)):
                outString = outString + separator + toPrintFormat[i] % fitsfile[toPrint[i]]

            print outString + separator

# -------------------------------------------------------------------
def listInDirNsamp(separator=' ', nsamp_febe='', outfile=''):
    """
    DES: list the input directory. Print number of samples. 
    INP: (str) field separator
         e.g. separator = '||' for moinmoin-style table
         (str) nsamp_febe : febe for counting number of samples, if wanted
         (str) outfile    : name of output file, if output should be written to file
    """

    global CurrentList
    if CurrentList in [[]]:
        MessHand.warning('list empty, populating...')
        findInDir(nsamp_febe=nsamp_febe)

    #toPrint        = ['Object','ScanType','NSubscans','ExpTime','FEBEList',\
    #                  'refChan','MJD' ,'Size']
    #toPrintFormat  = ['%15s'  ,'%10s'    ,'%3i'      ,'%7.0f'  ,'%17s'    ,\
    #                  '%5s'    ,'%9.3f','(%5.1f MB)']

    toPrint        = ['Scannum','Object','ScanType','ScanMode','NSubscans',\
                      'NSamp','MJD','FEBEList','refChan' ,'Size']
    toPrintFormat  = ['%6i'  ,'%31s'    ,'%12s'      ,'%12s'  ,'%5i'   ,\
                      '%10i','%12.4f', '%35s', '%12s', '(%5.1f MB)']

    if outfile:
        f=file(outfile,'w')

    for fitsfile in CurrentList:
        if fitsfile['status']:
            
            # Reformat the filename for display
            displayName = fitsfile['filename']
            if displayName[-5:] == '.fits':
                displayName = displayName[:-5]
            if displayName[-8:] == '.fits.gz':
                displayName = displayName[:-8]
            if displayName[:5] == 'APEX-' and displayName[-len(BoaConfig.projID):] == BoaConfig.projID:
                displayName = displayName[5:-len(BoaConfig.projID)-1]

            #outString = separator + "%7s " % displayName
            outString = ''

            for i in range(len(toPrint)):
                outString = outString + separator + toPrintFormat[i] % fitsfile[toPrint[i]]

            if outfile:
                f.write(outString + separator + '\n')
            else:
                print outString + separator

    if outfile:
        f.close()




# -------------------------------------------------------------------
def resetCurrentList():
    """
    DES: reset the CurrentList to the complete List
    """
    findInDir(update=1)


# -------------------------------------------------------------------
def selectInDir(type, value, test='eq'):
    """
    DES: Make a selection in the current List
    INP: (s)  type : on what type should we make the selection
         (s) value : the value of that type
         (s)  test : the test to do ( default: eq -- lt, gt, le, ge)
    """
    global CurrentList
    List = []
    for fitsfile in CurrentList:
        if fitsfile['status']:
            fitsfileType = fitsfile[type]
            if test == 'eq' and fitsfileType == value:
                List.append(fitsfile)
            if test == 'lt' and fitsfileType < value:
                List.append(fitsfile)
            if test == 'gt' and fitsfileType > value:
                List.append(fitsfile)
            if test == 'le' and fitsfileType <= value:
                List.append(fitsfile)
            if test == 'ge' and fitsfileType >= value:
                List.append(fitsfile)
            if test == 'ne' and fitsfileType != value:
                List.append(fitsfile)

    CurrentList = List
    MessHand.info("selected %i scan(s) with %s %s %s"%(len(List),type,test,value)) 

# -------------------------------------------------------------------
def removeScans(scanList):
    """
    DES: Remove scans from the current list
    INP: (l)  scanList : list of scans to remove
    """
    global CurrentList

    for i in range(len(CurrentList)):
        if CurrentList[i]['Scannum'] in scanList:
            CurrentList[i]['status'] = 0

# -------------------------------------------------------------------
def findInDir(update=0, nsamp_febe=''):
    """
    DES: Build the list of readable FITS files located in the input directory
    INP: (l) update : should we look only for new scans? (def.: NO)
    (str) nsamp_febe : specify a febe combination to derive the number of samples from
                            the corresponding arraydata tables
    """

    global ListDir, CurrentList

    if not update:
        ListDir=[]
    
    # Keep track of what is in the ListDir right now
    filenameInListDir = []
    for fitsfile in ListDir:
        filenameInListDir.append(fitsfile['filename'])

    listing = os.listdir(BoaConfig.inDir)
    fitsfile = [filename for filename in listing \
                if BoaMBFits.isDataset(BoaConfig.inDir+filename)]

    for filename in fitsfile:
        MessHand.debug("Processing "+filename)

        if filename not in filenameInListDir:
            fitsdesc = {'filename': filename}
            fitsdesc['status'] = 1
            try:
                checkFits(filename)

                dataset = BoaMBFits.importDataset(BoaConfig.inDir+filename)
                reader = BoaMBFitsReader.createReader(dataset)
                
                reader.openSubscan(subsnum=None)

                if nsamp_febe:
                    tablesARRAYDATA = dataset.getTables(EXTNAME='ARRAYDATA-MBFITS')
                    nsamp=0
                    for tab in tablesARRAYDATA:
                        tab.open()
                        k=tab.getKeyword("FEBE")
                        febe=k.getValue()
                        if (febe == nsamp_febe):
                            nsamp +=  tab.getNumRows()
                    fitsdesc['NSamp']      = nsamp
                else:
                    fitsdesc['NSamp']      = -999


                fitsdesc['MBFitsVer']  = reader.read('MBFitsVer')
                fitsdesc['Instrument'] = reader.read('Instrument')
                fitsdesc['ExpTime']    = reader.read('ExpTime')
                
                fitsdesc['Telescope']  = reader.read('Telescope')
                fitsdesc['Project']    = reader.read('Project')
                fitsdesc['Observer']   = reader.read('Observer')
                fitsdesc['Scannum']    = reader.read('ScanNum')
                fitsdesc['Date']       = reader.read('DateObs')
                fitsdesc['MJD']        = reader.read('ScanMJD')
                fitsdesc['LST']        = reader.read('ScanLST')
                fitsdesc['NSubscans']  = reader.read('NObs')
                fitsdesc['Object']     = reader.read('Object')
                fitsdesc['ScanType']   = reader.read('ScanType')
                fitsdesc['ScanMode']   = reader.read('ScanMode')
                #fitsdesc['ScanGeom']   = reader.read('ScanGeom')
                    
                fitsdesc['FEBEList']   = reader.read('Febes')

                fitsdesc['FEBE'] = []
                fitsdesc['refChan'] = []
                for febe in fitsdesc['FEBEList']:
                    febedesc = {}

                    febedesc['FEBE']      = reader.read('Febe', febe=febe)
                    febedesc['nFeed']     = reader.read('FebeFeed', febe=febe)
                    #febedesc['type']      = reader.read('FdTypCod', febe=febe)
                    febedesc['type']      = reader.read('FeedCode', febe=febe)
                    febedesc['nUsedFeed'] = reader.read('NUseFeeds', febe=febe)
                    febedesc['refChan']   = reader.read("RefFeed",  febe=febe)
                    fitsdesc['FEBE'].append(febedesc)
                    fitsdesc['refChan'].append(febedesc['refChan'])

                fitsdesc['Size'] = dataset.getSize()/1024.**2  # in mega bytes

                dataset.close()

                ListDir.append(fitsdesc)
            except Exception, explain:
                fitsdesc['status'] = 0
                MessHand.warning(filename+
                                 " is a bad MBFits file ("+str(explain)+
                                 "), removed from list")
        
    # Reset the current list and sort by MJD:
    tempDir = {} # key: MJD; value: list of indeces in ListDir
    for i in xrange(len(ListDir)):
        mjd = ListDir[i]['MJD']
        if not mjd in tempDir.keys(): 
            tempDir[mjd] = []
        tempDir[mjd].append(i)
    keys = tempDir.keys()
    keys.sort()
   
    CurrentList = []
    for key in keys:
        indices = tempDir[key]
        for index in indices:
            CurrentList.append(ListDir[index])


# -------------------------------------------------------------------

def checkFits(filename):
    """
    DES: check for MBFits name structure
    INP: (string) filename : the complete name of the fitsfile
    """
    dataset = None
    try:
        dataset = BoaMBFits.importDataset(BoaConfig.inDir+filename)
    
        MessHand.debug("Check fits basic structure : "+str(len(dataset.getTables()))+" tables in dataset")
        
        # Get the HduNumber for primary and other tables
        tablesSCAN      = dataset.getTables(EXTNAME='SCAN-MBFITS')
        tablesFEBEPAR   = dataset.getTables(EXTNAME='FEBEPAR-MBFITS')
        tablesARRAYDATA = dataset.getTables(EXTNAME='ARRAYDATA-MBFITS')
        tablesMONITOR   = dataset.getTables(EXTNAME='MONITOR-MBFITS')
        tablesDATAPAR   = dataset.getTables(EXTNAME='DATAPAR-MBFITS')
        
        # Check the basic structure of the fitsfile
        if not tablesSCAN :
            dataset.close()
            raise BoaMBFits.MBFitsError("no SCAN-MBFITS tables")
        if len(tablesSCAN) > 1:
            dataset.close()
            raise BoaMBFits.MBFitsError("more than one SCAN-MBFITS tables")
        if not tablesFEBEPAR :
            dataset.close()        
            raise BoaMBFits.MBFitsError("no FEBEPAR-MBFITS tables")
        if not tablesARRAYDATA :
            dataset.close()
            raise BoaMBFits.MBFitsError("no ARRAYDATA-MBFITS tables")
        if not tablesMONITOR :
            dataset.close()
            raise BoaMBFits.MBFitsError("no MONITOR-MBFITS tables")
        if not tablesDATAPAR :
            dataset.close()
            raise BoaMBFits.MBFitsError("no DATAPAR-MBFITS tables")
    ##     if len(tablesARRAYDATA)*len(tablesFEBEPAR) \
    ##            != len(tablesMONITOR) \
    ##            != len(tablesDATAPAR)*len(tablesFEBEPAR):
    ##         raise BoaMBFits.MBFitsError("number of ARRAYDATA, MONITOR and DATAPAR does not match")
    	
    ##     hdu = mbfits.Table(Parent=fitsFile)
        
        # Primary Header
        kwMBFTSVER = dataset.getKeyword('MBFTSVER')
        if not kwMBFTSVER:
            dataset.close()
            raise BoaMBFits.MBFitsError("Keyword MBFITSVER is missing")
        if kwMBFTSVER.getValue() < 1.56:
            dataset.close()
            raise BoaMBFits.MBFitsError("Old MBFits version :"+str(kwMBFTSVER.getValue()))
    	
        # SCAN-MBFITS
        MessHand.debug("Check SCAN-MBIFTS :")
        tableSCAN = tablesSCAN[0]
        tableSCAN.open()
        
        kwOBJECT = tableSCAN.getKeyword('OBJECT')
        if (not kwOBJECT) or (kwOBJECT.getValue() == ''):
            #MessHand.info(str(filename)+" has no object name")
            raise BoaMBFits.MBFitsError(str(filename)+" has no object name")

        kwNOBS     = tableSCAN.getKeyword('NOBS')
        if not kwNOBS:
            dataset.close()
            raise BoaMBFits.MBFitsError("No NOBS keyword")
    
        kwNSUBS    = tableSCAN.getKeyword('NSUBS')
        if not kwNSUBS:
            dataset.close()
            raise BoaMBFits.MBFitsError("No NSUBS keyword")
    
        # From 1.57 one should only test for NSUBS, but since Fred is
        # a silly person and changed only the NOBS keyword..
        if kwNOBS.getValue() == 0 and kwNSUBS.getValue() == 0:
            dataset.close()
            raise BoaMBFits.MBFitsError("NOBS keyword set to 0")
    
        # check the silly bug of cfitsio
    #    if kwNSUBS.getValue() > 98:
    #        dataset.close()
    #        raise BoaMBFits.MBFitsError("More than 98 subscans, cfitsio can not (yet) handle that")
            
    
        kwDIAMETER = tableSCAN.getKeyword('DIAMETER')
        if (not kwDIAMETER) or (kwDIAMETER.getValue() == '') or (kwDIAMETER.getValue() == 0):
            dataset.close()
            raise BoaMBFits.MBFitsError("DIAMETER keyword missing")
    
        nRows = tableSCAN.getNumRows()
        if nRows == 0:
            dataset.close()
            raise BoaMBFits.MBFitsError("No FEBEPAR in SCAN-MBFITS")
        if nRows != len(tablesFEBEPAR):
            dataset.close()
            raise BoaMBFits.MBFitsError("SCAN-MBFITS mismatch the number of FEBEPAR")
    
        allFEBE = tableSCAN.getColumn('FEBE').read()
        SCANNUM = tableSCAN.getKeyword('SCANNUM').getValue()
    	    
        tableSCAN.close()
    
        # FEBEPAR-MBFITS
        MessHand.debug("Check FEBEPAR tables ("+str(len(tablesFEBEPAR))+"table(s))")
    
        for tableFEBEPAR in tablesFEBEPAR:
    ##         MessHand.debug("Check FEBEPAR ( hdu ="+str(numHdu)+")")
            tableFEBEPAR.open()
    
            if tableFEBEPAR.getNumRows() == 0:
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBEPAR table empty")
            kwFEBE   = tableFEBEPAR.getKeyword('FEBE')
            if (not kwFEBE) or (kwFEBE.getValue() in ['','None']):
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBE name not defined")
            if not kwFEBE.getValue() in allFEBE:
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBE not listed in SCAN-MBFITS")
            lSCANNUM = tableFEBEPAR.getKeyword('SCANNUM').getValue()
            if lSCANNUM != SCANNUM:
                dataset.close()
                raise BoaMBFits.MBFitsError("SCANNUM does not match SCAN-MBFITS")
    	    
            tableFEBEPAR.close()
            
        # ARRAYDATA-MBFITS
        MessHand.debug("Check ARRAYDATA tables ("+str(len(tablesARRAYDATA))+"table(s))")
    
        arraydata_empty = 1
    
        for tableARRAYDATA in tablesARRAYDATA:
    ##         MessHand.debug("Check ARRAYDATA ( hdu ="+str(numHdu)+")")
            tableARRAYDATA.open()
            
            if tableARRAYDATA.getNumRows() != 0:
                arraydata_empty = 0
            RESTFREQ = tableARRAYDATA.getKeyword('RESTFREQ').getValue()
            if RESTFREQ == 0:
                dataset.close()
                raise BoaMBFits.MBFitsError("RESTFREQ null")
            kwFEBE   = tableARRAYDATA.getKeyword('FEBE')
            if (not kwFEBE) or (kwFEBE.getValue() in ['','None']):
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBE name not defined")
            if not kwFEBE.getValue() in allFEBE:
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBE not listed in SCAN-MBFITS")
            lSCANNUM = tableARRAYDATA.getKeyword('SCANNUM').getValue()
            if lSCANNUM != SCANNUM:
                dataset.close()
                raise BoaMBFits.MBFitsError("SCANNUM does not match SCAN-MBFITS")
        
            tableARRAYDATA.close()
        
        if arraydata_empty:
            dataset.close()
            raise BoaMBFits.MBFitsError("All ARRAYDATA tables empty")
    
    	
        # DATAPAR-MBFITS
        MessHand.debug("Check DATAPAR tables ("+str(len(tablesDATAPAR))+"table(s))")
    
        datapar_empty = 1
        for tableDATAPAR in tablesDATAPAR:
    ##         MessHand.debug(" Check DATAPAR ( hdu ="+str(numHdu)+")")
            tableDATAPAR.open()
            
            if tableDATAPAR.getNumRows() != 0:
                datapar_empty = 0
            kwFEBE   = tableDATAPAR.getKeyword('FEBE')
            if (not kwFEBE) or (kwFEBE.getValue() in ['','None']):
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBE name not defined")
            if not kwFEBE.getValue() in allFEBE:
                dataset.close()
                raise BoaMBFits.MBFitsError("FEBE not listed in SCAN-MBFITS")
            lSCANNUM = tableDATAPAR.getKeyword('SCANNUM').getValue()
            if lSCANNUM != SCANNUM:
                dataset.close()
                raise BoaMBFits.MBFitsError("SCANNUM does not match SCAN-MBFITS")
    
            tableDATAPAR.close()
            
        if datapar_empty:
            dataset.close()
            raise BoaMBFits.MBFitsError("All DATAPAR tables empty")
    
        # Check that corresponding ARRAYDATA and DATAPAR tables have same numbers of rows
        nSubs = kwNSUBS.getValue()
        for iSubs in xrange(1, nSubs+1):
            for febe in allFEBE:
                tableARRAYDATA = dataset.getTables(EXTNAME="ARRAYDATA-MBFITS", SUBSNUM=iSubs, FEBE=febe)
                tableDATAPAR   = dataset.getTables(EXTNAME="DATAPAR-MBFITS", SUBSNUM=iSubs, FEBE=febe)
                tableARRAYDATA = tableARRAYDATA[0]
                tableDATAPAR   = tableDATAPAR[0]
    
    	        tableARRAYDATA.open()
                tableDATAPAR.open()
    
                nRows1 = tableARRAYDATA.getNumRows()
                nRows2 = tableDATAPAR.getNumRows()
    
    	        tableARRAYDATA.close()
                tableDATAPAR.close()
    
                if nRows1 != nRows2:
                    MessHand.warning(filename+" different number of rows in ARRAYDATA and DATAPAR"+\
                                     " for subscan %i - still readable"%(iSubs))
                                    
                kwRESTFREQ = tableARRAYDATA.getKeyword('RESTFREQ')
                if not kwRESTFREQ:
                    dataset.close()
                    raise BoaMBFits.MBFitsError("No RESTFREQ keyword")
                    
                if kwRESTFREQ.getValue() == 0:
                    dataset.close()
                    raise BoaMBFits.MBFitsError("RESTFREQ set to 0")
        # Close open files
        dataset.close()
    except:
        # Close open files
        if dataset:
            dataset.close()
        raise
    
