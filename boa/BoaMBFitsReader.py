"""
NAM: BoaMBFitsReader.py (file)
DES: Provides classes and methodes for the high level access to MBFits datasets.
     The Module Interface consists of 
     - Function createReader
     - Class MBFitsReader
     - Class ApexMBFitsReader
     - Class IramMBFitsReader
     - Class MBFitsReaderError

     MBFitsReader is the parent class of ApexMBFitsReader and IramMBFitsReader and
     contains the public interface for the subclasses. To read the contents of a MBFits
     dataset, use the concrete subclasses of MBFitsReader.
"""                    
# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------

import BoaMBFits
import Numeric
import arrayfns


# ---------------------------------------------------------------------
# ---- Function createReader ------------------------------------------
# ---------------------------------------------------------------------

def createReader(dataset):
    """
    DES: Creates a new MBFitsReader for the specified BoaMBFits.Dataset object dataset.
         Depending on the contents of the dataset, either an ApexMBFitsReader or an
         IramMBFitsReader will be created.
    INP: (BoaMBFits.Dataset) dataset: The dataset to be read by the created MBFitsReader.
    OUT: (MBFitsReader)             : The reader object of the proper concrete subtype of
                                      MBFitsReader.
    """
    try:
        kwMBFTSVER = dataset.getKeyword("MBFTSVER")
        if kwMBFTSVER:
            return ApexMBFitsReader(dataset)
        
        kwIMBFTSVE = dataset.getKeyword("IMBFTSVE")
        if kwIMBFTSVE:
            return IramMBFitsReader(dataset)
    except Exception, data:
        raise MBFitsReaderError(str(data))


# ---------------------------------------------------------------------
# ---- Class MBFitsReader ---------------------------------------------
# ---------------------------------------------------------------------

class MBFitsReader:
    """
    DES: Parent reader class for ApexMBFitsReader and IramMBFitsReader.
         Contains the public interface for the subclasses and common private methods.
         To read the contents of a MBFits dataset, use the concrete subclasses of this
         class.
    """
    
    def __init__(self, dataset):
        self._dataset = dataset
        
        # The dictionary that maps item keys to functions and arguments are
        # filled in the constructors of the subclasses
        self._itemKeyToFuncArgsMap = {}

    # -----------------------------------------------------------------

    def getType(self):
        return self.__class__.__name__
    
    # -----------------------------------------------------------------

    def getBlankInt(self):
        """
        DES: Returns blanking value for integers as used in MBFITS
        """
        return -999

    # -----------------------------------------------------------------

    def getBlankFloat(self):
        """
        DES: Returns blanking value for floats as used in MBFITS
        """
        return -999.

    # -----------------------------------------------------------------

    def openSubscan(self, subsnum=None):
        """
        DES: Opens all tables related with the specified subscan.
        INP: (int) subsnum: The number of the subscan to be opened.
                            If subsnum = None (the default) all tables
                            that are related with the scan itself instead
                            of a subscan are opened.
        OUT: (int)        : The number of tables opened
        """
        try:
            if subsnum==None:
                tables = self._dataset.getTables(OBSNUM=None)
            else:
                tables = self._dataset.getTables(OBSNUM=subsnum)
            
            nOpenedTables = 0
            for table in tables:
                if not table.isOpen():
                    nOpenedTables += 1
                    table.open()
            return nOpenedTables
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def closeSubscan(self, subsnum=None):
        """
        DES: Closes all (open) tables related with the specified subscan.
        INP: (int) subsnum: The number of the subscan to be closed.
                            If subsnum = None (the default) all tables
                            that are related with the scan itself instead
                            of a subscan are closed.
        """
        try:
            if subsnum==None:
                tables = self._dataset.getTables(OBSNUM=None)
            else:
                tables = self._dataset.getTables(OBSNUM=subsnum)
            for table in tables:
                if table.isOpen():
                    table.close()
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def read(self, itemKey, **kargs):
        """
        DES: Reads item itemKey from the dataset, using the additional
             arguments in **kargs.
        INP: (string) itemKey: The key of the item to be read.
                               Inspect the init method of the concrete
                               subclasses for valid keys.
                      **kargs: Additional arguments necessary to read
                               the specified item. 
                               Inspect the init method of the concrete
                               subclasses for necessary additional arguments.
        OUT:                 : The read item.
                               May be of virtually every data type depending
                               on the datatype and location of the item in the
                               MBFitsFile.
        """
        try:
            # Get the information from _itemKeyToFuncArgsMap:
            fa =  self._itemKeyToFuncArgsMap[itemKey]
            func = fa["function"]
            setArgs = fa["setArgs"]
            addArgs = fa["addArgs"]

            # Check that all necessary argumets are passed:
            for addArg in addArgs:
                if not addArg in kargs.keys():
                    raise "Argument '%s' missing" % addArg
                setArgs[addArg] = kargs[addArg]

            # Retrieve the data using the function specified in
            # _itemKeyToFuncArgsMap:
            result = apply(func, [], setArgs)
            return result

        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readKeyword(self, key, extname=None, **kargs):
        # Read a keyword:
        try:
            if not extname:
                kw = self._dataset.getKeyword(key)
            else:
                tbl = self._getTable(extname=extname, **kargs)
                kw = tbl.getKeyword(key)
            if not kw:
                raise "Could not read keyword"
            result = kw.getValue()
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))
        
    # -----------------------------------------------------------------

    def _readKeywordTuple(self, keys, extname=None, **kargs):
        # Read a tuple of keywords:
        try:
            result = []
            for key in keys:
                result.append(self._readKeyword(key=key, \
                                                extname=extname, \
                                                **kargs))
            return tuple(result)
        except Exception, data:
            raise MBFitsReaderError(str(data))
        
    # -----------------------------------------------------------------

    def _readColumn(self, colname, extname, **kargs):
        # Read a column:
        try:
            tbl = self._getTable(extname=extname, **kargs)
            col = tbl.getColumn(colname)
            if not col:
                raise "Could not read column"
            result = col.read()
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _getTable(self, **kargs):
        try:
            kargsT = {}
            if "extname" in kargs.keys():
                kargsT["EXTNAME"] = kargs["extname"]
            if "subsnum" in kargs.keys():
                kargsT["OBSNUM"] = kargs["subsnum"]
            if "febe" in kargs.keys():
                kargsT["FEBE"] = kargs["febe"]
            if "baseband" in kargs.keys():
                kargsT["BASEBAND"] = kargs["baseband"]
            
            tbls = self._dataset.getTables(**kargsT)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            return tbls[0]
        except Exception, data:
            raise MBFitsReaderError(str(data))
        
# ---------------------------------------------------------------------
# ---- Class MBFitsReaderError ----------------------------------------
# ---------------------------------------------------------------------

class MBFitsReaderError(Exception):

    def __init__(self, msg):
        self.msg = msg
        
    # -----------------------------------------------------------------

    def __str__(self):
        return `self.msg`




# ---------------------------------------------------------------------
# ---- Class ApexMBFitsReader -----------------------------------------
# ---------------------------------------------------------------------

class ApexMBFitsReader(MBFitsReader):
    """
    DES: Reader class for (APEX) MBFITS 1.60 and earlier.
         Consult the documentation of the superclass MBFitsReader
         and the source code of the init method to find out what the
         class does.
    """
    
    def __init__(self, dataset):
        MBFitsReader.__init__(self, dataset)

        # The following dictionary maps each itemKey to
        # a method of MBFitsReader together with information
        # about arguments of the call to this method.

        # Each entry into the dictionary must be of the form:
        # itemKey: {"function": function,
        #           "setArgs": setArgsDict,
        #           "addArgs": addArgsList}

        # itemKey:  The key to be used in the read method
        # function: The method of MBFitsReader which is called
        #           by the read method
        # setArgs:  A dictionary that specifies names and values of
        #           argumets to be passed to function
        # addArgs:  A list with the names of arguments to be passed
        #           to function. The requested data can only be
        #           retrieved from the dataset, if the values of
        #           these arguments are known at runtime.
        #           In contrast to setArgs, the values
        #           of these arguments cannot be specified during
        #           implementation time. Hence, the values must be passed to
        #           the read method.
        
        self._itemKeyToFuncArgsMap = {

            # Primary header keywords:

            "ExpTime":    {"function": self._readKeyword,
                           "setArgs": {"key": "EXPTIME"},
                           "addArgs": []},

            "Instrument": {"function": self._readKeyword,
                           "setArgs": {"key": "INSTRUME"},
                           "addArgs": []},

            "MBFitsVer":  {"function": self._readKeyword,
                           "setArgs": {"key": "MBFTSVER"},
                           "addArgs": []},

            "StartAz":    {"function": self._readKeyword,
                           "setArgs": {"key": "HIERARCH ESO TEL AZ"},
                           "addArgs": []},

            "StartAlt":   {"function": self._readKeyword,
                           "setArgs": {"key": "HIERARCH ESO TEL ALT"},
                           "addArgs": []},

            "StartRa":    {"function": self._readKeyword,
                           "setArgs": {"key": "RA"},
                           "addArgs": []},

            "StartDec":   {"function": self._readKeyword,
                           "setArgs": {"key": "DEC"},
                           "addArgs": []},



            # SCAN-MBFITS keywords:

            "Coord":   {"function": self._readKeywordTuple,
                          "setArgs": {"keys": ("CRVAL1", "CRVAL2"),
                          #"setArgs": {"keys": ("BLONGOBJ", "BLATOBJ"),
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "Basis":   {"function": self._readKeywordTuple,
                          "setArgs": {"keys": ("CTYPE1", "CTYPE2"),
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "Equinox":   {"function": self._readKeyword,
                          "setArgs": {"key": "EQUINOX",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "DateObs":   {"function": self._readKeyword,
                          "setArgs": {"key": "DATE-OBS",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "Diameter":  {"function": self._readKeyword,
                          "setArgs": {"key": "DIAMETER",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanLST":   {"function": self._readKeyword,
                          "setArgs": {"key": "LST",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanMJD":   {"function": self._readKeyword,
                          "setArgs": {"key": "MJD",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "NObs":      {"function": self._readKeyword,
                          "setArgs": {"key": "NOBS",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "NSubs":     {"function": self._readKeyword,
                          "setArgs": {"key": "NSUBS",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "Object":    {"function": self._readKeyword,
                          "setArgs": {"key": "OBJECT",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "Observer":  {"function": self._readKeyword,
                          "setArgs": {"key": "OBSID",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "DeltaCA":   {"function": self._readKeyword,
                           "setArgs": {"key": "PDELTACA",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "DeltaIE":   {"function": self._readKeyword,
                           "setArgs": {"key": "PDELTAIE",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "Project":   {"function": self._readKeyword,
                           "setArgs": {"key": "PROJID",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "ScanGeom":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANGEOM",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanLen":   {"function": self._readKeyword,
                          "setArgs": {"key": "SCANLEN",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanLine":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANLINE",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanMode":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANMODE",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanNum":   {"function": self._readKeyword,
                          "setArgs": {"key": "SCANNUM",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanPar1":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANPAR1",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanPar2":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANPAR2",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanType":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANTYPE",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanXsp":   {"function": self._readKeyword,
                          "setArgs": {"key": "SCANXSPC",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanXVel":  {"function": self._readKeyword,
                          "setArgs": {"key": "SCANXVEL",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "ScanYsp":   {"function": self._readKeyword,
                          "setArgs": {"key": "SCANYSPC",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "SiteElev":  {"function": self._readKeyword,
                          "setArgs": {"key": "SITEELEV",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "SiteLat":   {"function": self._readKeyword,
                          "setArgs": {"key": "SITELAT",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "SiteLong":  {"function": self._readKeyword,
                          "setArgs": {"key": "SITELONG",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "Telescope": {"function": self._readKeyword,
                          "setArgs": {"key": "TELESCOP",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},

            "WobCycle":   {"function": self._readKeyword,
                           "setArgs": {"key": "WOBCYCLE",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "WobMode":    {"function": self._readKeyword,
                           "setArgs": {"key": "WOBMODE",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "WobThrow":   {"function": self._readKeyword,
                           "setArgs": {"key": "WOBTHROW",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "WobUsed":    {"function": self._readKeyword,
                           "setArgs": {"key": "WOBUSED",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "TAIUTC":     {"function": self._readKeyword,
                           "setArgs": {"key": "TAI2UTC",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},

            "UTCUT1":     {"function": self._readKeyword,
                           "setArgs": {"key": "UTC2UT1",
                                       "extname": "SCAN-MBFITS"},
                           "addArgs": []},


            # The following are used for writing .dat files used for
            # computing pointing model
            "PDeltaCA":    {"function": self._readKeyword,
                            "setArgs": {"key": "PDELTACA",
                                        "extname": "SCAN-MBFITS"},
                            "addArgs": []},
            "PDeltaIE":    {"function": self._readKeyword,
                            "setArgs": {"key": "PDELTAIE",
                                        "extname": "SCAN-MBFITS"},
                            "addArgs": []},
            "FDeltaCA":    {"function": self._readKeyword,
                            "setArgs": {"key": "FDELTACA",
                                        "extname": "SCAN-MBFITS"},
                            "addArgs": []},
            "FDeltaIE":    {"function": self._readKeyword,
                            "setArgs": {"key": "FDELTAIE",
                                        "extname": "SCAN-MBFITS"},
                            "addArgs": []},

            # SCAN-MBFITS columns:

            "Febes":     {"function": self._readColumn,
                          "setArgs": {"colname": "FEBE",
                                      "extname": "SCAN-MBFITS"},
                          "addArgs": []},




            # FEBEPAR-MBFITS keywords:

            "Colstart":  {"function": self._readKeyword,
                          "setArgs": {"key": "CARX",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "FdTypCod":  {"function": self._readKeyword,
                          "setArgs": {"key": "FDTYPCOD",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "Febe":      {"function": self._readKeyword,
                          "setArgs": {"key": "FEBE",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "FebeFeed":  {"function": self._readKeyword,
                          "setArgs": {"key": "FEBEFEED",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "Nula":      {"function": self._readKeyword,
                          "setArgs": {"key": "IARX",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "Nule":      {"function": self._readKeyword,
                          "setArgs": {"key": "IERX",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "FeedCode":  {"function": self._readKeyword,
                          "setArgs": {"key": "FDTYPCOD",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "DewUser":  {"function": self._readKeyword,
                          "setArgs": {"key": "DEWUSER",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "RefOffX":  {"function": self._readKeyword,
                          "setArgs": {"key": "REFOFFX",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},

            "RefOffY":  {"function": self._readKeyword,
                          "setArgs": {"key": "REFOFFY",
                                      "extname": "FEBEPAR-MBFITS"},
                          "addArgs": ["febe"]},


            # FEBEPAR-MBFITS columns:

            "BolCalFc":  {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "BOLCALFC"},
                          "addArgs": ["febe"]},

            "FeedOffX":  {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "FEEDOFFX"},
                          "addArgs": ["febe"]},

            "FeedOffY":  {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "FEEDOFFY"},
                          "addArgs": ["febe"]},

            "FlatField": {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "FLATFIEL"},
                          "addArgs": ["febe"]},

            "NUseFeeds": {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "NUSEFEED"},
                          "addArgs": ["febe"]},

            "RefFeed":   {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "REFFEED"},
                          "addArgs": ["febe"]},

            "UseBand":   {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "USEBAND"},
                          "addArgs": ["febe"]},

            "UseFeed":   {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "USEFEED"},
                          "addArgs": ["febe"]},

            "FeedType":  {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "FEEDTYPE"},
                          "addArgs": ["febe"]},

            "DCoffset":  {"function": self._readColumnFebepar,
                          "setArgs": {"colname": "BOLDCOFF"},
                          "addArgs": ["febe"]},

            # DATAPAR-MBFITS keywords:

            "SubsStart": {"function": self._readKeyword,
                          "setArgs": {"key": "DATE-OBS",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "ObsNum":    {"function": self._readKeyword,
                          "setArgs": {"key": "OBSNUM",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "ObsType":   {"function": self._readKeyword,
                          "setArgs": {"key": "OBSTYPE",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "ObsStatus": {"function": self._readKeyword,
                          "setArgs": {"key": "OBSTATUS",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "DewExtra": {"function": self._readKeyword,
                          "setArgs": {"key": "DEWEXTRA",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            # DATAPAR-MBFITS columns:

            "Az":       {"function": self._readColumn,
                          "setArgs": {"colname": "AZIMUTH",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "BasLat":       {"function": self._readColumn,
                          "setArgs": {"colname": "BASLAT",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "BasLon":       {"function": self._readColumn,
                          "setArgs": {"colname": "BASLONG",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "El":       {"function": self._readColumn,
                          "setArgs": {"colname": "ELEVATIO",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "LatOff":   {"function": self._readColumn,
                          "setArgs": {"colname": "LATOFF",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "LonOff":    {"function": self._readColumn,
                          "setArgs": {"colname": "LONGOFF",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "LST":       {"function": self._readColumn,
                          "setArgs": {"colname": "LST",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "MJD":       {"function": self._readColumn,
                          "setArgs": {"colname": "MJD",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "LatPole":       {"function": self._readColumn,
                          "setArgs": {"colname": "MLATPOLE",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "LonPole":       {"function": self._readColumn,
                          "setArgs": {"colname": "MLONPOLE",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "RotAngle":  {"function": self._readColumn,
                          "setArgs": {"colname": "ROTANGLE",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "MVAL1":  {"function": self._readColumn,
                          "setArgs": {"colname": "MCRVAL1",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "MVAL2":  {"function": self._readColumn,
                          "setArgs": {"colname": "MCRVAL2",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "Ra":     {"function": self._readColumn,
                          "setArgs": {"colname": "RA",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

            "Dec":    {"function": self._readColumn,
                          "setArgs": {"colname": "DEC",
                                      "extname": "DATAPAR-MBFITS"},
                          "addArgs": ["subsnum", "febe"]},

#             "WobblerSta": {"function": self._readColumn,
#                           "setArgs": {"colname": "PHASE",
#                                       "extname": "DATAPAR-MBFITS"},
#                           "addArgs": ["subsnum", "febe"]},



            # ARRAYDATA-MBFITS keywords:
            
            "RestFreq":  {"function": self._readKeyword,
                          "setArgs": {"key": "RESTFREQ",
                                      "extname": "ARRAYDATA-MBFITS"},
                          "addArgs": ["subsnum", "febe", "baseband"]},

            "NUseFeed":  {"function": self._readKeyword,
                          "setArgs": {"key": "NUSEFEED",
                                      "extname": "ARRAYDATA-MBFITS"},
                          "addArgs": ["subsnum", "febe", "baseband"]},

            # ARRAYDATA-MBFITS columns:

            "Data":      {"function": self._readData,
                          "setArgs": {},
                          "addArgs": ["subsnum", "febe", "baseband"]},


            # MONITOR-MBFITS keywords:
            
            "UsrFrame":  {"function": self._readKeyword,
                          "setArgs": {"key": "USRFRAME",
                                      "extname": "MONITOR-MBFITS"},
                          "addArgs": ["subsnum"]},

            # Custom:

            "FocX":      {"function": self._readDFocus,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum", "febe"]},

            "FocY":      {"function": self._readDFocus,
                          "setArgs": {"mpIndex": 1},
                          "addArgs": ["subsnum", "febe"]},

            "FocZ":      {"function": self._readDFocus,
                          "setArgs": {"mpIndex": 2},
                          "addArgs": ["subsnum", "febe"]},

            "PhiX":      {"function": self._readDPhi,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum", "febe"]},

            "PhiY":      {"function": self._readDPhi,
                          "setArgs": {"mpIndex": 1},
                          "addArgs": ["subsnum", "febe"]},

            "PhiZ":      {"function": self._readDPhi,
                          "setArgs": {"mpIndex": 2},
                          "addArgs": ["subsnum", "febe"]},

            "Subscans":  {"function": self._readSubscans,
                          "setArgs": {},
                          "addArgs": []},

            "NInteg":    {"function": self._readNInteg,
                          "setArgs": {},
                          "addArgs": ["subsnum", "febe", "baseband"]},

            "Refract":   {"function": self._readRefract,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum"]},

            "TempHe":   {"function": self._readHeTemp,
                          "setArgs": {},
                          "addArgs": ["subsnum"]},
            
            "WindSpeedDir": {"function": self._readWindSpeedDir,
                             "setArgs": {},
                             "addArgs": ["subsnum"]},
            
            "AzEl0":   {"function": self._readZeroPos,
                        "setArgs": {},
                        "addArgs": ["subsnum"]},

            "T_amb":   {"function": self._readTamb,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum"]},

            "BiasAmplitude": {"function": self._readBiasAmpl,
                              "setArgs": {"mpIndex": 0},
                              "addArgs": ["subsnum"]},

            "QdacAmplitude": {"function": self._readQdacAmpl,
                              "setArgs": {"mpIndex": 0},
                              "addArgs": ["subsnum"]},

            "BiasPotsetting": {"function": self._readBiasPotset,
                               "setArgs": {"mpIndex": 0},
                               "addArgs": ["subsnum"]},

            "GainSetting": {"function": self._readGainSetting,
                               "setArgs": {"mpIndex": 0},
                               "addArgs": ["subsnum"]},

            "PWV":            {"function": self._readPWV,
                               "setArgs": {},
                               "addArgs": ["subsnum"]},
            
            "WobblerSta": {"function": self._readWobblerSta,
                          "setArgs": {},
                          "addArgs": ["subsnum", "febe"]},

            "ScanDir":   {"function": self._readScanDir,
                          "setArgs": {},
                          "addArgs": ["subsnum", "febe"]},

            }

    # -----------------------------------------------------------------

    def _readColumnFebepar(self, colname, febe):
        # Read a column from a FEBEPAR table:

        # FEBEPAR tables are special in the following respects:
        # - They contain only one line; hence, it can be argued that the data from
        #   these columns should be arrays (or lists) of dimension 1 less than
        #   ordinary columns. This is implemented here.
        # - They contain a variaty of scalar, 1-dim, and 2-dim data, not allways
        #   with the correct TDIM keywords.

        # This make a seperate routine to read these columns necessary

        colsScalar = ["REFFEED"]
        
        try:
            tbl = self._getTable(extname="FEBEPAR-MBFITS", febe=febe)
            col = tbl.getColumn(colname)
            if not col:
                raise "Could not read column"

            rows = col.read()
            shape = col.getDim()

            row = rows[0]

            if colname in colsScalar:
                try:
                    result = row[0]
                except:
                    result = row
            else:
                if type(row) == type(Numeric.array([])):
                    result = row
                else:
                    try:
                        result = Numeric.array(row)
                    except:
                        result = row

            if shape:
                try:
                    result.shape = shape
                except:
                    pass

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readData(self, subsnum, febe, baseband):
        # Special routine to read DATA column from ARRAYDATA table 
        # Necessary since the ordering of data changed in May 2006
        try:

            nChannels = self._readKeyword("CHANNELS", \
                                          extname="ARRAYDATA-MBFITS", \
                                          subsnum=subsnum, \
                                          febe=febe, \
                                          baseband=baseband)

            data = self._readColumn(colname="DATA", \
                                    extname="ARRAYDATA-MBFITS", \
                                    subsnum=subsnum, \
                                    febe=febe, \
                                    baseband=baseband)

            if nChannels == data.shape[2]:
                # Since 2006/5/05 the spectral channels are the 3rd dim of data: 
                nInteg, nUsfd, nCh = data.shape
                if nCh == 1:
                    result = data
                    result.shape = (nInteg, nUsfd)
                else:
                    # More than 1 spectral channel:
                    result = data[:,:,0]
            else:
                # Before 2006/5/05 the spectral channels are the 2nd dim of data: 
                nInteg, nCh, nUsfd = data.shape
                if nCh == 1:
                    result = data
                    result.shape = (nInteg, nUsfd)
                else:
                    # More than 1 spectral channel:
                    result = data[:,0,::]

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDFocus(self, subsnum, febe, mpIndex):
        # Special routine to read Focus data necessary since
        # - data was originally stored in column DFOCUS in
        #   DATAPAR, but in columns DFOCUS_X, DFOCUS_Y, DFOCUS_Z
        #   later
        # - special case for FOCUS Scans: Here, the values from
        #   the MONITOR table are returned.
        try:
            # DFocusZ is used to decide whether the values
            # from DATAPAR or MONITOR is used:
            focusZ = self._readDFocusDatapar(subsnum=subsnum, \
                                             febe=febe, \
                                             mpIndex=2)

            # Read from DATAPAR:
            if mpIndex==2:
                result = focusZ
            else:
                result = self._readDFocusDatapar(subsnum=subsnum, \
                                                 febe=febe, \
                                                 mpIndex=mpIndex)

            scanType = self.read("ScanType")
            if scanType.find('FOCUS') > -1:
                if focusZ[0] == self.getBlankFloat():
                    # Read first value from MONITOR and resize it to correct shape:
                    monFocusXYZ = self._readMonValue(monpoint="DFOCUS_X_Y_Z",
                                                     subsnum=subsnum)
                    monFocus = monFocusXYZ[0][mpIndex]
                    result = Numeric.array([monFocus]*len(result),
                                           typecode=result.typecode())

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDFocusDatapar(self, subsnum, febe, mpIndex):
        try:
            tbls = self._dataset.getTables(EXTNAME="DATAPAR-MBFITS", \
                                           OBSNUM=subsnum, \
                                           FEBE=febe)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            if "DFOCUS" in tbls[0].getColumnNames():
                # Old format: All three values in single column:
                dfocus = self._readColumn(colname="DFOCUS", \
                                          extname="DATAPAR-MBFITS", \
                                          subsnum=subsnum, \
                                          febe=febe)
                try:
                    result = dfocus[:,mpIndex]
                except:
                    # If blanking value, the column contains only 1-dim data:
                    result = dfocus
            else:
                if mpIndex==0:
                    result = self._readColumn(colname="DFOCUS_X", \
                                              extname="DATAPAR-MBFITS", \
                                              subsnum=subsnum, \
                                              febe=febe)
                elif mpIndex==1:
                    result = self._readColumn(colname="DFOCUS_Y", \
                                              extname="DATAPAR-MBFITS", \
                                              subsnum=subsnum, \
                                              febe=febe)
                else:
                    result = self._readColumn(colname="DFOCUS_Z", \
                                              extname="DATAPAR-MBFITS", \
                                              subsnum=subsnum, \
                                              febe=febe)
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDPhi(self, subsnum, febe, mpIndex):
        # Special routine to read Phi data necessary since
        # - data was originally stored in column DPHI in
        #   DATAPAR, but in columns DPHI_X, DPHI_Y, DPHI_Z
        #   later
        # - special case for FOCUS Scans: Here, the values from
        #   the MONITOR table are returned.
        try:
            # DFocusZ is used to decide whether the values
            # from DATAPAR or MONITOR is used:
            focusZ = self._readDFocusDatapar(subsnum=subsnum, \
                                             febe=febe, \
                                             mpIndex=2)

            # Read from DATAPAR:
            result = self._readDPhiDatapar(subsnum=subsnum, \
                                           febe=febe, \
                                           mpIndex=mpIndex)

            scanType = self.read("ScanType")
            if scanType.find('FOCUS') > -1:
                if focusZ == self.getBlankFloat():
                    # Read frist value from MONITOR and resize it to correct shape:
                    monPhiXYZ = self._readMonValue(monpoint="DPHI_X_Y_Z",
                                                   subsnum=subsnum)
                    monPhi = monPhiXYZ[0][mpIndex]
                    result = Numeric.array([monPhi]*len(result),
                                           typecode=result.typecode())

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDPhiDatapar(self, subsnum, febe, mpIndex):
        try:
            tbls = self._dataset.getTables(EXTNAME="DATAPAR-MBFITS", \
                                           OBSNUM=subsnum, \
                                           FEBE=febe)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            if "DPHI" in tbls[0].getColumnNames():
                # Old format: All three values in single column:
                dphi = self._readColumn(colname="DPHI", \
                                        extname="DATAPAR-MBFITS", \
                                        subsnum=subsnum, \
                                        febe=febe)
                try:
                    result = dphi[:,mpIndex]
                except:
                    # If blanking value, the column contains only 1-dim data:
                    result = dphi
            else:
                if mpIndex==0:
                    result = self._readColumn(colname="DPHI_X", \
                                              extname="DATAPAR-MBFITS", \
                                              subsnum=subsnum, \
                                              febe=febe)
                elif mpIndex==1:
                    result = self._readColumn(colname="DPHI_Y", \
                                              extname="DATAPAR-MBFITS", \
                                              subsnum=subsnum, \
                                              febe=febe)
                else:
                    result = self._readColumn(colname="DPHI_Z", \
                                              extname="DATAPAR-MBFITS", \
                                              subsnum=subsnum, \
                                              febe=febe)
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readWobblerSta(self, subsnum, febe):
        # Special routine necessary since data was originally
        # stored in column ISWITCH in DATAPAR, while later
        # column PHASE was used.
        # In the latter case, PHASE must be translated to the corresponding
        # string.
        try:
            tbls = self._dataset.getTables(EXTNAME="DATAPAR-MBFITS", \
                                           OBSNUM=subsnum, \
                                           FEBE=febe)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            if "ISWITCH" in tbls[0].getColumnNames():
                # Old format: String in column:
                result = self._readColumn(colname="ISWITCH", \
                                          extname="DATAPAR-MBFITS", \
                                          subsnum=subsnum, \
                                          febe=febe)
                for i in xrange(len(result)):
                    if result[i] == "NONE":
                        result[i] = ""
            else:
                # New format: Phase number in column:
                phases = self._readColumn(colname="PHASE", \
                                          extname="DATAPAR-MBFITS", \
                                          subsnum=subsnum, \
                                          febe=febe)
                # Need to read the phase description (which one is ON,
                # which one is OFF). Defined in SCAN header, or re-defined
                # in DATAPAR header if different for that specific subscan
                # Try first to get it from DATAPAR
                phas1 = tbls[0].getKeyword('PHASE1')
                phas2 = tbls[0].getKeyword('PHASE2')
                if (not phas1) or (not phas2):
                    # at least one is missing, search in SCAN header
                    scantbl = self._dataset.getTables(EXTNAME="SCAN-MBFITS")
                    phas1 = scantbl[0].getKeyword('PHASE1')
                    phas2 = scantbl[0].getKeyword('PHASE2')

                phas1 = phas1.getValue()
                phas2 = phas2.getValue()
                    
                result = []
                thePhases = []
                if phas1 == 'WON':
                    thePhases.append('ON')
                elif phas1 == 'WOFF':
                    thePhases.append('OFF')
                else:
                    thePhases.append('')
                if phas2 == 'WON':
                    thePhases.append('ON')
                elif phas2 == 'WOFF':
                    thePhases.append('OFF')
                else:
                    thePhases.append('')
                
                # Translate:
                for phase in phases:
                    if phase == self.getBlankInt():
                        result.append("")
                    elif phase == 1 or phase == 2:
                        result.append(thePhases[phase-1])
                    else:
                        result.append('')
                                
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------
    def _readScanDir(self, subsnum, febe):
        # Another special case where the information can be either
        # in the SCAN header, or in DATAPAR headers if it changes
        # from subscan to subscan
        try:
            tbls = self._dataset.getTables(EXTNAME="DATAPAR-MBFITS", \
                                           OBSNUM=subsnum, \
                                           FEBE=febe)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            # Try first to get values from DATAPAR
            onedir = tbls[0].getKeyword('SCANDIR')
            if not onedir:
                # search in SCAN header
                scantbl = self._dataset.getTables(EXTNAME="SCAN-MBFITS")
                onedir = scantbl[0].getKeyword('SCANDIR')
                
            return onedir.getValue()
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readMonValue(self, subsnum, monpoint,getTime = 0):
        try:
            tbls = self._dataset.getTables(EXTNAME="MONITOR-MBFITS", \
                                           OBSNUM=subsnum)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"
            col = tbls[0].getColumn("MONVALUE")
            if not col:
                raise "Could not read column"
            if getTime:
                colTime = tbls[0].getColumn("MJD")
            tbls[0].setSelection("MONPOINT=='%s'" % monpoint)
            result = col.read()
            if getTime:
                resTime = colTime.read()
                result = (result,resTime)
            tbls[0].clearSelection()

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readSubscans(self):
        # Subscan numbers from MONITOR tables
        try:
            subscans = []
            tbls = self._dataset.getTables(EXTNAME="MONITOR-MBFITS")
            for tbl in tbls:
                # Use OBSNUM instead of SUBSNUM, since some files have
                # SUBSNUM=0 for all subscans
                kw = tbl.getKeyword("OBSNUM")
                if not kw:
                    raise "Could not read keyword 'OBSNUM'"
                subscans.append(kw.getValue())
            return subscans
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readNInteg(self, subsnum, febe, baseband):
        # Number of integrations from ARRAYDATA table
        try:
            tbls = self._dataset.getTables(EXTNAME="ARRAYDATA-MBFITS", \
                                           OBSNUM=subsnum,
                                           FEBE=febe,
                                           BASEBAND=baseband)
            if len(tbls) == 0:
                msg  = "No table selected for "
                msg += "EXTNAME=%s " % str("ARRAYDATA-MBFITS")
                msg += "OBSNUM=%s " % str(subsnum)
                msg += "FEBE=%s " % str(febe)
                msg += "BASEBAND=%s " % str(baseband)
                raise msg
            elif len(tbls) > 1:
                msg  = "No unique table selected for "
                msg += "EXTNAME=%s " % str("ARRAYDATA-MBFITS")
                msg += "OBSNUM=%s " % str(subsnum)
                msg += "FEBE=%s " % str(febe)
                msg += "BASEBAND=%s " % str(baseband)
                raise msg

            nbRowsArraydata = tbls[0].getNumRows()

            tbls = self._dataset.getTables(EXTNAME="DATAPAR-MBFITS", \
                                           OBSNUM=subsnum,
                                           FEBE=febe)
            if len(tbls) == 0:
                msg  = "No table selected for "
                msg += "EXTNAME=%s " % str("DATAPAR-MBFITS")
                msg += "OBSNUM=%s " % str(subsnum)
                msg += "FEBE=%s " % str(febe)
                raise msg
            elif len(tbls) > 1:
                msg  = "No unique table selected for "
                msg += "EXTNAME=%s " % str("DATAPAR-MBFITS")
                msg += "OBSNUM=%s " % str(subsnum)
                msg += "FEBE=%s " % str(febe)
                raise msg

            nbRowsDatapar = tbls[0].getNumRows()

            result = min((nbRowsArraydata,nbRowsDatapar))
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readRefract(self, subsnum, mpIndex):
        try:
            monRefrac = self._readMonValue(monpoint="REFRACTION",
                                           subsnum=subsnum)
            monRefrac = monRefrac[0][mpIndex]
            result = Numeric.array([monRefrac],
                                   typecode=Numeric.Float32)
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))
        
    # -----------------------------------------------------------------

    def _readTamb(self, subsnum, mpIndex):
        try:
            monTambPHum = self._readMonValue(monpoint="TAMB_P_HUMID",
                                             subsnum=subsnum)
            monTamb = monTambPHum[0][mpIndex]
            result = Numeric.array([monTamb],'f')
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))
        
    # -----------------------------------------------------------------

    def _readHeTemp(self, subsnum):
        try:
            monHeT,timeHe = self._readMonValue(monpoint="LABOCA_HE3TEMP",\
                                               subsnum=subsnum,\
                                               getTime=1)
            result = []
            times = []
            nb = len(timeHe)
            for i in range(nb):
                if monHeT[i][0] != self.getBlankFloat() and \
                       timeHe[i] != self.getBlankFloat():
                    result.append(monHeT[i][0])
                    # 32s shift in the timestamps!
                    times.append(timeHe[i]-32./86400.)
                    # TODO: remove this -32./86400. when fixed by Giorgio
                    
            times  = Numeric.array((times),'d')
            result = Numeric.array((result),'f')
            return result, times
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readWindSpeedDir(self, subsnum):
        try:
            windSpeedDir,timeWind = self._readMonValue(monpoint="WIND_SPEED_DIR",\
                                                       subsnum=subsnum,\
                                                       getTime=1)

            
            aWindSpeedDir = Numeric.array(windSpeedDir)
            return aWindSpeedDir[:,0], aWindSpeedDir[:,1], Numeric.array(timeWind)
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readGains(self, febe):

        tbl = self._getTable(extname="FEBEPAR-MBFITS", febe=febe)

        # Frontend Gain
        kwFeGain = tbl.getKeyword('FEGAIN')
        if kwFeGain:
            kwFeGain = kwFeGain.getValue()
        else:
            kwFeGain = 0.

        # Backend Gain
        try:
            colBeGain = self._readColumnFebepar(colname="BEGAIN", \
                                                febe=febe)
        except MBFitsReaderError:
            colBeGain = self._readColumnFebepar(colname="BOLREFGN", \
                                                febe=febe)
        if not colBeGain:
            raise "Could not read keyword BE Gain"

        if kwFeGain == self.getBlankFloat():
            kwFeGain = 0.
        if colBeGain == self.getBlankFloat():
            colBeGain = 1.

        return kwFeGain, colBeGain

    # -----------------------------------------------------------------
    def _readZeroPos(self,subsnum):

        """
        Read and return ANTENNA and ENCODE Az, El at scan start
        """
        #try:
        monEncoder = self._readMonValue(monpoint="ENCODER_AZ_EL",\
                                        subsnum=subsnum)
        monAntenna = self._readMonValue(monpoint="ANTENNA_AZ_EL",\
                                        subsnum=subsnum)
        return monAntenna[0],monEncoder[0]
        #except Exception, data:
        #    raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readBiasAmpl(self, subsnum, mpIndex):
        try:
            monBiasAmpl = self._readMonValue(monpoint="BOLOSZ_BIASAMPLITUDE",
                                             subsnum=subsnum)
            if monBiasAmpl:
                result = monBiasAmpl[mpIndex]
            else:
                result = None
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readQdacAmpl(self, subsnum, mpIndex):
        try:
            monQdacAmpl = self._readMonValue(monpoint="BOLOSZ_QDACAMPLITUDE",
                                             subsnum=subsnum)
            if monQdacAmpl:
                result = monQdacAmpl[mpIndex]
            else:
                result = None
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))


    # -----------------------------------------------------------------

    def _readBiasPotset(self, subsnum, mpIndex):
        try:
            monBiasPotset = self._readMonValue(monpoint="BOLOSZ_BIASPOTSETTING",
                                               subsnum=subsnum)
            if monBiasPotset:
                result = monBiasPotset[mpIndex]
            else:
                result = None
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readGainSetting(self, subsnum, mpIndex):
        try:
            monGainSetting = self._readMonValue(monpoint="BOLOSZ_DEMODATTEN",
                                               subsnum=subsnum)
            if monGainSetting:
                result = monGainSetting[mpIndex]
            else:
                result = None
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))     
        
    # -----------------------------------------------------------------

    def _readPWV(self, subsnum):
        try:
            monPWV = self._readMonValue(monpoint="PWV",
                                        subsnum=subsnum)
            return monPWV
        except Exception, data:
            raise MBFitsReaderError(str(data))

# ---------------------------------------------------------------------
# ---- Class IramMBFitsReader -----------------------------------------
# ---------------------------------------------------------------------

class IramMBFitsReader(MBFitsReader):
    """
    DES: Reader class for IRAM-MBFITS
         Consult the documentation of the superclass MBFitsReader
         and the source code of the init method to find out what the
         class does.
    """
    
    def __init__(self, dataset):
        MBFitsReader.__init__(self, dataset)
        
        # Two dictionary used to cache data for the interpolations:
        self._MJD = {} # key: subscan number; item: MJDs from IMBF-backend
        self._MJD_dap = {} # key: (subscan, febe); item: MJDs from IMBF-antenna

        # The following dictionary maps each itemKey to
        # a method of MBFitsReader together with information
        # about arguments of the call to this method.

        # Each entry into the dictionary must be of the form:
        # itemKey: {"function": function,
        #           "setArgs": setArgsDict,
        #           "addArgs": addArgsList}

        # itemKey:  The key to be used in the read method
        # function: The method of MBFitsReader which is called
        #           by the read method
        # setArgs:  A dictionary that specifies names and values of
        #           argumets to be passed to function
        # addArgs:  A list with the names of arguments to be passed
        #           to function. The requested data can only be
        #           retrieved from the dataset, if the values of
        #           these arguments are known at runtime.
        #           In contrast to setArgs, the values
        #           of these arguments cannot be specified during
        #           implementation time. Hence, the values must be passed to
        #           the read method.
        
        self._itemKeyToFuncArgsMap = {

            # Primary header keywords:

            "ExpTime":    {"function": self._readKeyword,
                           "setArgs": {"key": "EXPTIME"},
                           "addArgs": []},

            "Instrument": {"function": self._readKeyword,
                           "setArgs": {"key": "INSTRUME"},
                           "addArgs": []},

            "MBFitsVer":  {"function": self._readKeyword,
                           "setArgs": {"key": "IMBFTSVE"},
                           "addArgs": []},

            "DateObs":   {"function": self._readKeyword,
                          "setArgs": {"key": "DATE-OBS"},
                          "addArgs": []},

            "ScanLST":   {"function": self._readKeyword,
                          "setArgs": {"key": "LST"},
                          "addArgs": []},

            "ScanMJD":   {"function": self._readKeyword,
                          "setArgs": {"key": "MJD-OBS"},
                          "addArgs": []},

            "NObs":      {"function": self._readKeyword,
                          "setArgs": {"key": "N_OBSP"},
                          "addArgs": []},

            "NSubs":     {"function": self._readKeyword,
                          "setArgs": {"key": "N_OBSP"},
                          "addArgs": []},

            "Object":    {"function": self._readKeyword,
                          "setArgs": {"key": "OBJECT"},
                          "addArgs": []},

            "Telescope": {"function": self._readKeyword,
                          "setArgs": {"key": "TELESCOP"},
                          "addArgs": []},

            "ScanType":  {"function": self._readKeyword,
                          "setArgs": {"key": "OBSTYPE"},
                          "addArgs": []},

            # IMBF-scan keywords:

            "Coord":   {"function": self._readKeywordTuple,
                          "setArgs": {"keys": ("LONGOBJ", "LATOBJ"),
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "Basis":   {"function": self._readKeywordTuple,
                          "setArgs": {"keys": ("CTYPE1", "CTYPE2"),
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "Equinox":   {"function": self._readKeyword,
                          "setArgs": {"key": "EQUINOX",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "Diameter":  {"function": self._readKeyword,
                          "setArgs": {"key": "TELSIZE",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "Observer":  {"function": self._readKeyword,
                          "setArgs": {"key": "OBSID",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "DeltaCA":   {"function": self._readKeyword,
                           "setArgs": {"key": "P2COR",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "DeltaIE":   {"function": self._readKeyword,
                           "setArgs": {"key": "P7COR",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "Project":   {"function": self._readKeyword,
                           "setArgs": {"key": "PROJID",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "ScanNum":   {"function": self._readKeyword,
                          "setArgs": {"key": "SCANNUM",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "SiteElev":  {"function": self._readKeyword,
                          "setArgs": {"key": "SITEELEV",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "SiteLat":   {"function": self._readKeyword,
                          "setArgs": {"key": "SITELAT",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "SiteLong":  {"function": self._readKeyword,
                          "setArgs": {"key": "SITELONG",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},

            "WobCycle":   {"function": self._readKeyword,
                           "setArgs": {"key": "WOBCYCLE",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "WobMode":    {"function": self._readKeyword,
                           "setArgs": {"key": "WOBMODE",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "WobThrow":   {"function": self._readKeyword,
                           "setArgs": {"key": "WOBTHROW",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "WobUsed":    {"function": self._readKeyword,
                           "setArgs": {"key": "WOBUSED",
                                       "extname": "IMBF-scan"},
                           "addArgs": []},

            "RefFeed":   {"function": self._readKeyword,
                          "setArgs": {"key": "RECON",
                                      "extname": "IMBF-scan"},
                          "addArgs": []},


            # IMBF-scan columns:
            # nothing in the IMBF-scan table... take instrument from Primary

            "Febes":     {"function": self._readFebes,
                           "setArgs": {},
                           "addArgs": []},



            # IMBF-frontend keywords:

            "Febe":      {"function": self._readKeyword,
                          "setArgs": {"key": "FEBE",
                                      "extname": "IMBF-frontend"},
                          "addArgs": ["febe"]},

            "FebeFeed":  {"function": self._readKeyword,
                          "setArgs": {"key": "FEBEFEED",
                                      "extname": "IMBF-frontend"},
                          "addArgs": ["febe"]},

            "NUseFeeds": {"function": self._readKeyword,
                          "setArgs": {"key": "NUSEFEED",
                                      "extname": "IMBF-frontend"},
                          "addArgs": ["febe"]},

            # IMBF-frontend columns:

            "UseFeed":   {"function": self._readColumnFrontend,
                          "setArgs": {"colname": "USEFEED"},
                          "addArgs": ["febe"]},

            "FeedType":  {"function": self._readColumnFrontend,
                          "setArgs": {"colname": "FEEDTYPE"},
                          "addArgs": ["febe"]},

            "RestFreq":  {"function": self._readColumnFrontend,
                          "setArgs": {"colname": "RESTFREQ"},
                          "addArgs": ["febe"]},

            # IMBF-antenna keywords:

            "ScanMode":  {"function": self._readKeyword,
                          "setArgs": {"key": "SUBSTYPE",
                                      "extname": "IMBF-antenna",
                                      "subsnum": 1},
                          "addArgs": []},

            "SubsStart": {"function": self._readKeyword,
                          "setArgs": {"key": "DATE-OBS",
                                      "extname": "IMBF-antenna"},
                          "addArgs": ["subsnum"]},

            "ObsNum":    {"function": self._readKeyword,
                          "setArgs": {"key": "OBSNUM",
                                      "extname": "IMBF-antenna"},
                          "addArgs": ["subsnum"]},

            "ObsType":   {"function": self._readKeyword,
                          "setArgs": {"key": "OBSTYPE",
                                      "extname": "IMBF-antenna"},
                          "addArgs": ["subsnum"]},

            # IMBF-antenna columns:

            "Az":       {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "AZIMUTH"},
                          "addArgs": ["subsnum", "febe"]},

            "BasLat":       {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "BASLAT"},
                          "addArgs": ["subsnum", "febe"]},

            "BasLon":       {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "BASLONG"},
                          "addArgs": ["subsnum", "febe"]},

            "El":       {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "ELEVATION"},
                          "addArgs": ["subsnum", "febe"]},

            "LatOff":   {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "LATOFF"},
                          "addArgs": ["subsnum", "febe"]},

            "LonOff":    {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "LONGOFF"},
                          "addArgs": ["subsnum", "febe"]},

            "LST":       {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "LST"},
                          "addArgs": ["subsnum", "febe"]},

            "MJD-dap":   {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "MJD"},
                          "addArgs": ["subsnum", "febe"]},

            "RotAngle":  {"function": self._readColumnAntenna,
                          "setArgs": {"colname": "PARANGLE"},
                          "addArgs": ["subsnum", "febe"]},


            # ARRAYDATA-MBFITS keywords:
            
            "NUseFeed":  {"function": self._readKeywordBackend,
                          "setArgs": {"key": "CHANNELS"},
                          "addArgs": ["subsnum", "febe"]},

            # IMBF-backend columns:

            "Data":      {"function": self._readColumnBackend,
                          "setArgs": {"colname": "DATA"},
                          "addArgs": ["subsnum", "febe"]},

            "MJD":       {"function": self._readColumnBackend,
                          "setArgs": {"colname": "MJD"},
                          "addArgs": ["subsnum", "febe"]},
            
            "WobblerSta": {"function": self._readColumnBackend,
                          "setArgs": {"colname": "ISWITCH"}, 
                          "addArgs": ["subsnum", "febe"]},

            # Custom:

            "FocX":      {"function": self._readDFocus,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum", "febe"]},

            "FocY":      {"function": self._readDFocus,
                          "setArgs": {"mpIndex": 1},
                          "addArgs": ["subsnum", "febe"]},

            "FocZ":      {"function": self._readDFocus,
                          "setArgs": {"mpIndex": 2},
                          "addArgs": ["subsnum", "febe"]},

            "PhiX":      {"function": self._readDPhi,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum", "febe"]},

            "PhiY":      {"function": self._readDPhi,
                          "setArgs": {"mpIndex": 1},
                          "addArgs": ["subsnum", "febe"]},

            "PhiZ":      {"function": self._readDPhi,
                          "setArgs": {"mpIndex": 2},
                          "addArgs": ["subsnum", "febe"]},

            "Subscans":  {"function": self._readSubscans,
                          "setArgs": {},
                          "addArgs": []},

            "NInteg":    {"function": self._readNInteg,
                          "setArgs": {},
                          "addArgs": ["subsnum", "febe", "baseband"]},

            "Refract":   {"function": self._readRefract,
                          "setArgs": {"mpIndex": 0},
                          "addArgs": ["subsnum"]},
            
            # Some items not relevant for IRAM data
            "UseBand":   {"function": self._readUseBand,
                          "setArgs": {},
                          "addArgs": []},

            "FeedCode":  {"function": self._readFeedcode,
                          "setArgs": {},
                          "addArgs": []},

            "BolCalFc":  {"function": self._readBolCal,
                          "setArgs": {},
                          "addArgs": []},

            "UsrFrame":  {"function": self._readUsrFrame,
                          "setArgs": {},
                          "addArgs": []},

            }


    # -----------------------------------------------------------------

    def _readFebes(self):
        # Take instrument from Primary:
        febe = self._readKeyword("INSTRUME")
        return [febe]

    # -----------------------------------------------------------------

    def _readColumnAntenna(self, colname, subsnum, febe):
        # Fill the MJDs for the interpolation if not already present:
        if not (subsnum in self._MJD_dap.keys()):
            self._MJD_dap[subsnum] = self._readColumn(colname="MJD", \
                                                      extname="IMBF-antenna", \
                                                      subsnum=subsnum)   
        if not ((subsnum, febe) in self._MJD.keys()):
            self._MJD[(subsnum, febe)] = self.read("MJD", \
                                                   subsnum=subsnum, \
                                                   febe=febe).astype(Numeric.Float)   

        tmpArray = self._readColumn(colname=colname, \
                                    extname="IMBF-antenna", \
                                    subsnum=subsnum)
        if colname in ("AZIMUTH", "ELEVATION"):
            tmpArray = tmpArray[:,0,0]
        
        tmpInterp = arrayfns.interp(tmpArray, \
                                    self._MJD_dap[subsnum], \
                                    self._MJD[(subsnum, febe)])
        
        if colname in ("AZIMUTH", "ELEVATION", "LONGOFF", "LATOFF", "BASLONG", "BASLAT", "PARANGLE"):
            tmpInterp *= 180./Numeric.pi
        
        return tmpInterp
    
    # -----------------------------------------------------------------


    def _readColumnFrontend(self, colname, febe):
        # Read a column from the IMBF-frontend table:

        # The IMBF-frontend table is special in the following respects:
        # - It contains only one line; hence, it can be argued that the data from
        #   these columns should be arrays (or lists) of dimension 1 less than
        #   ordinary columns. This is implemented here.
        # - It contains a variaty of scalar, 1-dim, and 2-dim data, not allways
        #   with the correct TDIM keywords.

        # This make a seperate routine to read these columns necessary

        colsScalar = ["RESTFREQ"]
        
        try:
            tbl = self._getTable(extname="IMBF-frontend", febe=febe)
            col = tbl.getColumn(colname)
            if not col:
                raise "Could not read column"

            rows = col.read()
            row = rows[0]

            if colname in colsScalar:
                result = row[0][0]
            else:
                if type(row) == type(Numeric.array([])):
                    result = [row[:,0]]
                else:
                    result = [Numeric.array(row)]

            if colname=="RESTFREQ":
                # Scale from GHz to Hz
                result *= 1.e9

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------
    # IMBF-backendXXX table :
    # need special routines because table name depends on backend
    def _readKeywordBackend(self, key, subsnum, febe):
        tmp = febe.split('-')
        be = tmp[1]
        tabName = "IMBF-backend"+be
        result = self._readKeyword(key, extname=tabName, subsnum=subsnum)
        return result
    
    def _readColumnBackend(self, colname, subsnum, febe):
        tmp = febe.split('-')
        be = tmp[1]
        tabName = "IMBF-backend"+be
        try:
            tbl = self._getTable(extname=tabName, subsnum=subsnum)
            col = tbl.getColumn(colname)
            if not col:
                raise "Could not read column"

            rows = col.read()
            shape = col.getDim()

            if type(rows) == type(Numeric.array([])):
                result = rows
            else:
                try:
                    result = Numeric.array(rows)
                except:
                    result = rows

            if colname=="DATA":
                result = result[:,:,0]
                
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDFocus(self, subsnum, febe, mpIndex):
        # Special routine to read Focus data necessary since
        # - data was originally stored in column DFOCUS in
        #   DATAPAR, but in columns DFOCUS_X, DFOCUS_Y, DFOCUS_Z
        #   later
        # - special case for FOCUS Scans: Here, the values from
        #   the MONITOR table are returned.
        try:
            # DFocusZ is used to decide whether the values
            # from DATAPAR or MONITOR is used:
            focusZ = self._readDFocusDatapar(subsnum=subsnum, \
                                             febe=febe, \
                                             mpIndex=2)

            # Read from DATAPAR:
            if mpIndex==2:
                result = focusZ
            else:
                result = self._readDFocusDatapar(subsnum=subsnum, \
                                                 febe=febe, \
                                                 mpIndex=mpIndex)

            scanType = self.read("ScanType")
            if scanType.find('FOCUS') > -1:
                if focusZ[0] == self.getBlankFloat():
                    # Read frist value from MONITOR and resize it to correct shape:
                    monFocusXYZ = self._readMonValue(monpoint="DFOCUS_X_Y_Z",
                                                     subsnum=subsnum)
                    monFocus = monFocusXYZ[0][mpIndex]
                    result = Numeric.array([monFocus]*len(result),
                                           typecode=result.typecode())

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDFocusDatapar(self, subsnum, febe, mpIndex):
        try:
            tbls = self._dataset.getTables(EXTNAME="IMBF-antenna", \
                                           OBSNUM=subsnum, \
                                           FEBE=febe)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            if "DFOCUS" in tbls[0].getColumnNames():
                # Old format: All three values in single column:
                dfocus = self._readColumn(colname="DFOCUS", \
                                          extname="IMBF-antenna", \
                                          subsnum=subsnum, \
                                          febe=febe)
                try:
                    result = dfocus[:,mpIndex]
                except:
                    # If blanking value, the column contains only 1-dim data:
                    result = dfocus
            else:
                if mpIndex==0:
                    result = self._readColumn(colname="DFOCUS_X", \
                                              extname="IMBF-antenna", \
                                              subsnum=subsnum, \
                                              febe=febe)
                elif mpIndex==1:
                    result = self._readColumn(colname="DFOCUS_Y", \
                                              extname="IMBF-antenna", \
                                              subsnum=subsnum, \
                                              febe=febe)
                else:
                    result = self._readColumn(colname="DFOCUS_Z", \
                                              extname="IMBF-antenna", \
                                              subsnum=subsnum, \
                                              febe=febe)
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDPhi(self, subsnum, febe, mpIndex):
        # Special routine to read Phi data necessary since
        # - data was originally stored in column DPHI in
        #   DATAPAR, but in columns DPHI_X, DPHI_Y, DPHI_Z
        #   later
        # - special case for FOCUS Scans: Here, the values from
        #   the MONITOR table are returned.
        try:
            # DFocusZ is used to decide whether the values
            # from DATAPAR or MONITOR is used:
            focusZ = self._readDFocusDatapar(subsnum=subsnum, \
                                             febe=febe, \
                                             mpIndex=2)

            # Read from DATAPAR:
            result = self._readDPhiDatapar(subsnum=subsnum, \
                                           febe=febe, \
                                           mpIndex=mpIndex)

            scanType = self.read("ScanType")
            if scanType.find('FOCUS') > -1:
                if focusZ == self.getBlankFloat():
                    # Read frist value from MONITOR and resize it to correct shape:
                    monPhiXYZ = self._readMonValue(monpoint="DPHI_X_Y_Z",
                                                   subsnum=subsnum)
                    monPhi = monPhiXYZ[0][mpIndex]
                    result = Numeric.array([monPhi]*len(result),
                                           typecode=result.typecode())

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readDPhiDatapar(self, subsnum, febe, mpIndex):
        try:
            tbls = self._dataset.getTables(EXTNAME="IMBF-antenna", \
                                           OBSNUM=subsnum, \
                                           FEBE=febe)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"

            if "DPHI" in tbls[0].getColumnNames():
                # Old format: All three values in single column:
                dphi = self._readColumn(colname="DPHI", \
                                        extname="IMBF-antenna", \
                                        subsnum=subsnum, \
                                        febe=febe)
                try:
                    result = dphi[:,mpIndex]
                except:
                    # If blanking value, the column contains only 1-dim data:
                    result = dphi
            else:
                if mpIndex==0:
                    result = self._readColumn(colname="DPHI_X", \
                                              extname="IMBF-antenna", \
                                              subsnum=subsnum, \
                                              febe=febe)
                elif mpIndex==1:
                    result = self._readColumn(colname="DPHI_Y", \
                                              extname="IMBF-antenna", \
                                              subsnum=subsnum, \
                                              febe=febe)
                else:
                    result = self._readColumn(colname="DPHI_Z", \
                                              extname="IMBF-antenna", \
                                              subsnum=subsnum, \
                                              febe=febe)
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------
    def _readMonValue(self, subsnum, monpoint):
        try:
            tbls = self._dataset.getTables(EXTNAME="MONITOR-MBFITS", \
                                           OBSNUM=subsnum)
            if len(tbls) == 0:
                raise "No table selected"
            elif len(tbls) > 1:
                raise "No unique table selected"
            col = tbls[0].getColumn("MONVALUE")
            if not col:
                raise "Could not read column"
            tbls[0].setSelection("MONPOINT=='%s'" % monpoint)
            result = col.read()
            tbls[0].clearSelection()

            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readSubscans(self):
        # Subscan numbers from IMBF-antenna tables
        try:
            subscans = []
            tbls = self._dataset.getTables(EXTNAME="IMBF-antenna")
            for tbl in tbls:
                kw = tbl.getKeyword("OBSNUM")
                if not kw:
                    raise "Could not read keyword 'OBSNUM'"
                subscans.append(kw.getValue())
            return subscans
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readNInteg(self, subsnum, febe, baseband):
        # Number of integrations from backend table
        # split FEBE name to retrieve backend
        tmp = febe.split('-')
        be = tmp[1]
        tabName = "IMBF-backend"+be
        try:
            tbls = self._dataset.getTables(EXTNAME=tabName, \
                                           OBSNUM=subsnum)
            if len(tbls) == 0:
                msg  = "No table selected for "
                msg += "EXTNAME=%s " % str("ARRAYDATA-MBFITS")
                msg += "OBSNUM=%s " % str(subsnum)
                raise msg
            elif len(tbls) > 1:
                msg  = "No unique table selected for "
                msg += "EXTNAME=%s " % str("ARRAYDATA-MBFITS")
                msg += "OBSNUM=%s " % str(subsnum)
                raise msg
            result = tbls[0].getNumRows()
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))

    # -----------------------------------------------------------------

    def _readRefract(self, subsnum, mpIndex):
        try:
            monRefrac = self._readMonValue(monpoint="REFRACTION",
                                           subsnum=subsnum)
            monRefrac = monRefrac[0][mpIndex]
            result = Numeric.array([monRefrac],
                                   typecode=Numeric.Float32)
            return result
        except Exception, data:
            raise MBFitsReaderError(str(data))
        
    # -----------------------------------------------------------------

    def _readGains(self, febe):

        tbls = self._dataset.getTables(EXTNAME="IMBF-backendABBA")
        if len(tbls) == 0:
            raise "No table selected"
        elif len(tbls) > 1:
            tbls = tbls[0]

        # No Frontend Gain for MAMBO
        kwFeGain = 0.

        # Backend Gain
        kwBeGain = tbls.getKeyword('BOLGAIN')
        if kwBeGain:
            kwBeGain = kwBeGain.getValue()
        else:
            kwBeGain = 0.
        if not kwBeGain:
            raise "Could not read keyword BE Gain"

        if kwFeGain == self.getBlankFloat():
            kwFeGain = 0.
        if kwBeGain == self.getBlankFloat():
            kwBeGain = 1.

        return kwFeGain, kwBeGain

    def _readUseBand(self):
        # Only one band in MAMBO data
        return [0]

    def _readFeedcode(self):
        return '1:AC,2:DC,3:N'  # use fixed pseudo-code

    def _readBolCal(self):
        return 1.
    
    def _readUsrFrame(self):
        return "EQEQHO"
    
