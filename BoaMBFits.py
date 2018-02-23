"""
NAM: BoaMBFits.py (file)
DES: Provides classes and methodes for the low level access to MBFits datasets.
     The Module Interface consists of 
     - Function isDataset
     - Function importDataset
     - Function openDataset
     - Class Dataset
     - Class Table
     - Class Keyword
     - Class Column
     - Class ColumnInfo
     - Class MBFitsError
     Only these functions and classes should be used by clients!     

     The rest of the module contains classes that implement the module's
     functionality. These classes should not be used directly from outside the module,
     but only through the module's interface!
"""                    
# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------

import os
import Numeric
import cfitsio

# ---------------------------------------------------------------------
# ---- Module Interface -----------------------------------------------
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# ---- Function isDataset ---------------------------------------------
# ---------------------------------------------------------------------

def isDataset(datasetName):
    """
    NAM: isDataset (Function)
    DES: Checks if datasetName corresponds to a dataset
    INP: datasetName (str) : The full path to the directory or file to be checked
    OUT: (long)            : 1, if datasetName is a dataset, 0 else
    """
    if os.path.isdir(datasetName):
        return os.path.isfile(datasetName + "/GROUPING.fits") or \
               os.path.isfile(datasetName + "/GROUPING.fits.gz")
    else:
        return (datasetName[-5::] == ".fits" or datasetName[-8::] == ".fits.gz")

# ---------------------------------------------------------------------
# ---- Function importDataset -----------------------------------------
# ---------------------------------------------------------------------

def importDataset(datasetName, iomode=0):
    """
    NAM: importDataset (Function)
    DES: Creates a new Dataset object and imports it from disk files specified by
         datasetName. The return Dataset object is opened for reading (iomode=0) or
         reading and writing (iomode=1); however, the Dataset's tables are not opened.
    INP: datasetName (str): The full path of the directory or file of the dataset
         iomode     (long): 0: the returned Dataset is open for reading
                            1: the returned Dataset is open for reading and writing
    OUT: (Dataset)        : The imported Dataset object    
    """
    datasetImpl = _Dataset(datasetName)
    if os.path.isdir(datasetName):
        if os.path.isfile(datasetName + "/GROUPING.fits"):
            datasetImpl.imprt(datasetName + "/GROUPING.fits", iomode)
        elif os.path.isfile(datasetName + "/GROUPING.fits.gz"):
            datasetImpl.imprt(datasetName + "/GROUPING.fits.gz", iomode)
        else:
            raise MBFitsError("Direcotry %s no Dataset" % datasetName)
    else:
        datasetImpl.imprt(datasetName, iomode)

    return Dataset(datasetImpl)

# ---------------------------------------------------------------------
# ---- Function createDataset -----------------------------------------
# ---------------------------------------------------------------------

def createDataset(datasetName, filename, keywords, groupName):
    """
    NAM: createDataset (Function)
    DES: Creates a new Dataset object without Tables. The created dataset
         is opened for reading and writing.
    INP: datasetName       (str): The full path of the directory or file of the dataset
         filename          (str): The full name of the file containig the primary header;
                                  may be identical to datasetName
         keywords (Keyword list): List of Keyword objects to be added to the dataset
         groupName         (str): Name of the hierarchical group defined in the dataset's
                                  Grouping Table. If None, no Grouping Table is created.
    OUT: (Dataset)              : The created Dataset object
    """
    datasetImpl = _Dataset(datasetName)
    datasetImpl.construct(filename, keywords, groupName)
    return Dataset(datasetImpl)

# ---------------------------------------------------------------------
# ---- Class Dataset --------------------------------------------------
# ---------------------------------------------------------------------

class Dataset:
    """
    NAM: Dataset (Class)
    DES: Represents a MBFits dataset.
    """
    def __init__(self, implementation):
        """
        NAM: Dataset.__init__ (Method)
        DES: Constructor of Class Dataset
             Do not use this constructor from outside the module; instead, create
             Dataset objects via the functions importDataset or createDataset
        """
        self._implementation = implementation

    # -----------------------------------------------------------------

    def isOpen(self):
        """
        NAM: Dataset.isOpen (Method)
        DES: Check if Dataset is open for reading
        OUT: (long): 1, if Dataset is open for reading, 0 otherwise
        """
        if not self._implementation:
            return 0
        return self._implementation.isOpen()
    
    # -----------------------------------------------------------------

    def isWriteOpen(self):
        """
        NAM: Dataset.isWriteOpen (Method)
        DES: Check if Dataset is open for reading and writing
        OUT: (long): 1, if Dataset is open for reading and writing, 0 otherwise
        """
        if not self._implementation:
            return 0
        return self._implementation.isWriteOpen()
    
    # -----------------------------------------------------------------

    def open(self, iomode=0):
        """
        NAM: Dataset.open (Method)
        DES: Open the Dataset for reading or reading plus writing.
             A Dataset must be open to perform Dataset.getTables.
             A Dataset must be open for reading and writing in order to
             perform Dataste.addTable.
             Tables can only be opened if the corresponding Dataset is open.
             Note that this method does not open the Dataset's Tables.
        INP: iomode (long): 0: open Dataset for reading
                            1: open Dataset for reading and writing
                            Note that the Dataset's Tables can only be opened
                            for reading and writing if the Dataset was opened
                            with iomode=1
        """
        if self._implementation:
            self._implementation.open(iomode=iomode)

    # -----------------------------------------------------------------

    def close(self):
        """
        NAM: Dataset.close (Method)
        DES: Close the Dataset and all its Tables.
        """
        if self._implementation:
            self._implementation.close()

    # -----------------------------------------------------------------

    def exprt(self):
        """
        NAM: Dataset.exprt (Method)
        DES: Export the Dataset to disk.
             Export writes the keywords of the Primary Header and of the
             Grouping Table to disk; in addition, it creates the necessary
             columns of the Grouping Table if necessary.
             Note that this method does not export the Dataset's Tables.
        """
        if self._implementation:
            self._implementation.exprt()

    # -----------------------------------------------------------------

    def getName(self):
        """
        NAM: Dataset.getName (Method)
        DES: Returns the Dataset's name 
        OUT: (str): The full path of the directory or file of the Dataset
        """
        if not self._implementation:
            return None
        return self._implementation.getName()

    # -----------------------------------------------------------------

    def getSize(self):
        """
        NAM: Dataset.getSize (Method)
        DES: Returns the sum of the sizes of all files that make up the Dataset.
        OUT: (long): The size (in byte) of the Dataset.
        """
        if not self._implementation:
            return 0
        return self._implementation.getSize()

    # -----------------------------------------------------------------

    def getKeywordNames(self):
        """
        NAM: Dataset.getKeywordNames (Method)
        DES: Returns the names of the Keywords in the Dataset's Primary Header
             in the correct order
        OUT: (str list): The names of Keywords in the Dataset's Primary Header
        """
        if not self._implementation:
            return None
        return self._implementation.getKeywordNames()
    
    # -----------------------------------------------------------------

    def getKeyword(self, keyname):
        """
        NAM: Dataset.getKeyword (Method)
        DES: Returns a Keyword from the Dataset's Primary Header
        INP: keyname (str): The keyname of the requested Keyword
        OUT: (Keyword)    : The specified Keyword object from the Primary Header
        """
        if not self._implementation:
            return None
        return self._implementation.getKeyword(keyname)
    
    # -----------------------------------------------------------------

    def getTables(self, **keywords):
        """
        NAM: Dataset.getTables (Method)
        DES: Returns Tables from the Dataset.
             getTables(keyname1=value1, keyname2=value2,...) returns all Tables with
             matching keywords, so that for each returned Table t
             t.getKeyword("keyword1")==value1 and t.getKeyword("keyword2")==value2 and ...
             is true.
             Without argument, getTables() returns all Tables.
        INP: **keywords  : Arbitrary number of keyword arguments of type
                           keyname=keyvalue.
        OUT: (Table list): The specified Table objects
        """
        if not self._implementation:
            return None
        tableImpls = self._implementation.getTables(**keywords)
        tables = []
        for tableImpl in tableImpls:
            if tableImpl:
                tables.append(Table(tableImpl))
        return tables

    # -----------------------------------------------------------------

    def addTable(self, filename=None, keywords=[], colinfos=[]):
        """
        NAM: Dataset.addTable (Method)
        DES: Adds a Table to the Dataset.
             The added Table is created in the file specified by filename. If necessary,
             the file is created (including dirctories along its path). If no filename is 
             given, the new table is created in the file that holds the Dataset's
             Primary Header.
             The new Table will have Keywords and Columns as specified by the arguments keywords
             and colinfos. Note that in the course of creation of the Table, some additional
             Keywords may be created automatically; these Keywords can be accessed via
             Table.getKeyword exactly as the Keywords specified explicitely in
             Dataset.addTable.
             During execution of Dataset.addTable, the new Table is added to the Dataset's
             Grouping Table, if this exists. The new Table will be open for reading and
             writing, and is exported during execution of Dataset.addTable.
             The Dataset must be open for reading and writing.
        INP: filename (str)            : The full path of the file in which the new Table is created.
                                         If None, the new Table is created in the file with the
                                         Dataset's Primary Header.
             keywords (Keyword list)   : Keywords to be added to the new Table.
             colinfos (ColumnInfo list): Colinfo objects that describe the columns of the 
                                         new Table.
        OUT: (Table)                   : The newly created Table
        """
        if not self._implementation:
            return None
        tableImpl = self._implementation.addTable(filename, keywords, colinfos)
        if not tableImpl:
            return None
        table = Table(tableImpl)
        return table

# ---------------------------------------------------------------------
# ---- Class Table ----------------------------------------------------
# ---------------------------------------------------------------------

class Table:
    """
    NAM: Table (Class)
    DES: Represents a table of a MBFits dataset. 
    """
    def __init__(self, implementation):
        """
        NAM: Table.__init__ (Method)
        DES: Constructor of Class Table
             Do not use this constructor from outside the module; instead, create
             Table objects via the methods Dataset.getTables and .addTable
        """
        self._implementation = implementation

    # -----------------------------------------------------------------

    def isOpen(self):
        """
        NAM: Table.isOpen (Method)
        DES: Check if Table is open for reading
        OUT: (long): 1, if Table is open for reading, 0 otherwise
        """
        if not self._implementation:
            return 0
        return self._implementation.isOpen()

    # -----------------------------------------------------------------

    def isWriteOpen(self):
        """
        NAM: Table.isWriteOpen (Method)
        DES: Check if Table is open for reading and writing
        OUT: (long): 1, if Table is open for reading and writing, 0 otherwise
        """
        if not self._implementation:
            return 0
        return self._implementation.isWriteOpen()
    
    # -----------------------------------------------------------------

    def open(self, iomode=0):
        """
        NAM: Table.open (Method)
        DES: Open the Table for reading or reading plus writing.
             A Table must be open to perform Table.close, .getColumn,
             .getNumRows, and .setSelection.
             A Table must be open for reading and writing in order to
             perform Table.exprt.
             A Table can only be opened if the corresponding Dataset is open.
             A Table can only be opened for reading and writing if the corresponding
             Dataset is open for reading and writing.
        INP: iomode (long): 0: open Table for reading
                            1: open Table for reading and writing
                            Note that a Table can only be opened
                            for reading and writing if the Dataset was opened
                            with iomode=1
        """
        if self._implementation:
            self._implementation.open(iomode=iomode)

    # -----------------------------------------------------------------

    def close(self):
        """
        NAM: Table.close (Method)
        DES: Close the Table.
        """
        if self._implementation:
            self._implementation.close()

    # -----------------------------------------------------------------

    def exprt(self):
        """
        NAM: Table.exprt (Method)
        DES: Export the Table to disk.
             Export writes the Table's keywords to disk; in addition, it creates 
             the Table's columns if necessary.
             The Table must be open for reading and writing to call this method.
        """
        if self._implementation:
            self._implementation.exprt()

    # -----------------------------------------------------------------

    def getKeywordNames(self):
        """
        NAM: Table.getKeywordNames (Method)
        DES: Returns the names of the Keywords in the correct order
        OUT: (str list): The names of Keywords in the Table
        """
        if not self._implementation:
            return None
        return self._implementation.getKeywordNames()
    
    # -----------------------------------------------------------------

    def getKeyword(self, keyname):
        """
        NAM: Table.getKeyword (Method)
        DES: Returns a Keyword from the Table
        INP: keyname (str): The keyname of the requested Keyword
        OUT: (Keyword)    : The specified Keyword object
        """
        if not self._implementation:
            return None
        return self._implementation.getKeyword(keyname)
    
    # -----------------------------------------------------------------

    def getNumColumns(self):
        """
        NAM: Table.getNumColumns (Method)
        DES: Returns the number of Columns
        OUT: (long): The number of Columns in the Table
        """
        if not self._implementation:
            return None
        return self._implementation.getNumColumns()
 
    # -----------------------------------------------------------------

    def getColumnNames(self):
        """
        NAM: Table.getColumnNames (Method)
        DES: Returns the names of the Columns in the correct order.
        OUT: (str list): The names of Columns in the Table
        """
        if not self._implementation:
            return None
        return self._implementation.getColumnNames()

    # -----------------------------------------------------------------

    def getColumn(self, colName):
        """
        NAM: Table.getColumn (Method)
        DES: Returns the Column object with name colName.
             The Table must be open to call this method.
        INP: colName (str): The name of the requesetd Column
        OUT: (Column):      The Column object
        """
        if not self._implementation:
            return None
        columnImpl = self._implementation.getColumn(colName)
        if not columnImpl:
            return None
        column = Column(columnImpl)
        return column
   
    # -----------------------------------------------------------------

    def getNumRows(self):
        """
        NAM: Table.getNumRows (Method)
        DES: Returns the number of rows in the Table's columns.
             The Table must be open to call this method.
        OUT: (long): The number of rows in the Table's columns.
        """
        if not self._implementation:
            return None
        return self._implementation.getNumRows()
   
    # -----------------------------------------------------------------

    def getOptimalRowsize(self):
        """
        NAM: Table.getOptimalRowsize (Method)
        DES: Returns the optimal number of rows to be read or written in a single call
             of Column.read or .write for this Table.
             This number depends on the number of files that are open at the time of
             reading or writing. Hence for optimal performance, call this function after
             opening or closing a Table.
        OUT: (long): The optimal number of rows for reading and writing.
        """
        if not self._implementation:
            return None
        return self._implementation.getOptimalRowsize()

    # -----------------------------------------------------------------

    def setSelection(self, expression, firstRow=1, numRows=0):
        """
        NAM: Table.setSelection (Method)
        DES: Sets a selection to the Table.
             If a selection is set to a Table, subsequent calls of
             Column.read will return only those rows, where the specified boolean
             expression is true.
             For a full description of expression, firstRow and numRows, see
             the documentation of the CFITSIO routine fits_find rows.
             The Table must be open to call this method.
        INP: expression (str): Boolean expression which defines the selection
             firstRow (long):  First Row to which the selection applies.
             numRows(long):    Number of rows, to which the selection applies;
                               if None or 0, the selection applies to all rows
                               after firstRow
        """
        if self._implementation:
            self._implementation.setSelection(expression, firstRow, numRows)

    # -----------------------------------------------------------------

    def clearSelection(self):
        """
        NAM: Table.clearSelection (Method)
        DES: Clears a selection
        """
        if self._implementation:
            self._implementation.clearSelection()

    # -----------------------------------------------------------------

    def hasSelection(self):
        """
        NAM: Table.hasSelection (Method)
        DES: Check if a Table has a selection
        OUT: 1 if true; 0 else
        """
        if not self._implementation:
            return None
        return self._implementation.hasSelection()


# ---------------------------------------------------------------------
# ---- Class Keyword --------------------------------------------------
# ---------------------------------------------------------------------

class Keyword:
    """
    NAM: Keyword (Class)
    DES: Represents a keyword of a MBFits dataset
    """
    def __init__(self, name='', datatype='', value=None, format='', comment='', unit=''):
        """
        NAM: Keyword.__init__
        DES: Constructor of class Keyword
        INP: name (str):     The Keyword's name
             datatype (str): The Keywords datatype coded in the spirit of CFITSIO:
                             'A' - String
                             'L' - Logical
                             'I' - Integer
                             'J' - Long
                             'E' - Float
                             'D' - Double
             value:          The Keyword's value of the appropriate datatype.
             format (long):  The Keywords format used when writing the Keyword to disk.
                             For datatype
                             'A': The maximum length of the string
                             'E': The number of decimal digits
                             'D': The number of decimal digits
                             For all other datatypes unused.
             comment (str):  The comment of the Keyword; Note that it may or may not contain
                             the string describing the unit.
             unit (str):     The unit of the keyword.
        """
        self._name = name
        self._datatype = datatype
        self._value = value
        self._format = format
        self._comment = comment
        self._unit = unit
    
    # -----------------------------------------------------------------

    def __str__(self):
        return str((self._name, \
                    self._value, \
                    self._unit, \
                    self._format, \
                    self._comment, \
                    self._datatype))

    # -----------------------------------------------------------------

    def getName(self):
        """
        NAM: Keyword.getName
        DES: Retuns the keyword's name
        """
        return self._name
    
    # -----------------------------------------------------------------

    def getDatatype(self):
        """
        NAM: Keyword.getDatatype
        DES: Retuns the keyword's datatype
        """
        return self._datatype
    
    # -----------------------------------------------------------------

    def getValue(self):
        """
        NAM: Keyword.getValue
        DES: Retuns the keyword's value
        """
        return self._value
    
    # -----------------------------------------------------------------

    def getFormat(self):
        """
        NAM: Keyword.getFormat
        DES: Retuns the keyword's format
        """
        return self._format
    
    # -----------------------------------------------------------------

    def getComment(self):
        """
        NAM: Keyword.getComment
        DES: Retuns the keyword's comment
        """
        return self._comment
    
    # -----------------------------------------------------------------

    def getUnit(self):
        """
        NAM: Keyword.getUnit
        DES: Retuns the keyword's unit
        """
        return self._unit
    
    # -----------------------------------------------------------------

    def setValue(self, value):
        """
        NAM: Keyword.setValue
        DES: Sets the keyword's value
        """
        self._value = value
    
    # -----------------------------------------------------------------

    def setFormat(self, format):
        """
        NAM: Keyword.setFormat
        DES: Sets the keyword's format
        """
        self._format = format
    
    # -----------------------------------------------------------------

    def setComment(self, comment):
        """
        NAM: Keyword.setComment
        DES: Sets the keyword's comment
        """
        self._comment = comment
    
    # -----------------------------------------------------------------

    def setUnit(self, unit):
        """
        NAM: Keyword.setUnit
        DES: Sets the keyword's unit
        """
        self._unit = unit
    

# ---------------------------------------------------------------------
# ---- Class Column ---------------------------------------------------
# ---------------------------------------------------------------------

class Column:
    """
    NAM: Column (Class)
    DES: Represents a column of a MBFits dataset
    """
    def __init__(self, implementation):
        """
        NAM: Column.__init__ (Method)
        DES: Constructor of Class Column
             Do not use this constructor from outside the module; instead, create
             Column objects via the methods Table.getColumn
        """
        self._implementation = implementation
        
    # -----------------------------------------------------------------

    def getColnum(self):
        """
        NAM: Column.getColnum (Method)
        DES: Returns the column number. The first number of a Table has column number 1
        OUT: (long): The column number
        """
        if not self._implementation:
            return None
        return self._implementation.getColnum()
    
    # -----------------------------------------------------------------

    def getName(self):
        """
        NAM: Column.getName (Method)
        DES: Returns the column's name.
        OUT: (str): The column's nane
        """
        if not self._implementation:
            return None
        return self._implementation.getName()
    
    # -----------------------------------------------------------------

    def getDatatype(self):
        """
        NAM: Column.getDatatype (Method)
        DES: Returns the column's datatype.
             See Class ColumnInfo for further documentation.
        OUT: (str): The column's datatype
        """
        if not self._implementation:
            return None
        return self._implementation.getDatatype()
    
    # -----------------------------------------------------------------

    def getRepeat(self):
        """
        NAM: Column.getRepeat (Method)
        DES: Returns the column's repeat count. 
             See Class ColumnInfo for further documentation.
        OUT: (str): The column's repeat count
        """
        if not self._implementation:
            return None
        return self._implementation.getRepeat()

    # -----------------------------------------------------------------

    def getDim(self):
        """
        NAM: Column.getDim (Method)
        DES: Returns the column's dimension.
             See Class ColumnInfo for further documentation.
        OUT: (str): The column's dimension
        """
        if not self._implementation:
            return None
        return self._implementation.getDim()

    # -----------------------------------------------------------------

    def getDescription(self):
        """
        NAM: Column.getDescription (Method)
        DES: Returns the column's description.
             See Class ColumnInfo for further documentation.
        OUT: (str): The column's description.
        """
        if not self._implementation:
            return None
        return self._implementation.getDescription()

    # -----------------------------------------------------------------

    def getUnit(self):
        """
        NAM: Column.getUnit (Method)
        DES: Returns the column's unit.
             See Class ColumnInfo for further documentation.
        OUT: (str): The column's unit.
        """
        if not self._implementation:
            return None
        return self._implementation.getUnit()

    # -----------------------------------------------------------------

    def read(self, firstRow=1, numRows=0):
        """
        NAM: Column.read (Method)
        DES: Reads data from a Column.
             Reads numRows rows starting at row firstRow. If numRows=0, all rows
             starting from firstRow are read. Note that row numbers start at 1.
             The datatype of the output depends on the Column's datatype, dimension,
             and repeat count. Whenever possible, a Numeric.array is returned with the
             relation of the column's datatype and the Numeric.array's typecode as follows:
                 'L'                        -> Numeric.Int0
                 'B'                        -> Numeric.UnsignedInteger8
                 'I', 'I', 'UI'             -> Numeric.Int16
                 'K', 'UK', 'J', 'JJ', 'UJ' -> Numeric.Int32
                 'E', 'F'                   -> Numeric.Float32
                 'D', 'G'                   -> Numeric.Float64
             The shape of the returned Numeric array depends on the column's dimension; in general,
             the rank of the array is 1 higher than the rank of the column's dimension. This
             is also true, if only one row was read.
             For the column's datatype = 'A', a list of strings is returned.
             For variable length arrays, a list of 1-dim Numeric.arrays of the appropriate size
             is returned.
        INP: firstRow (long): The number of the first Row to be read
             numRows (long):  The number of rows to be read. If 0, all rows starting from firstRow
                              are read.
        OUT: See DES
        """
        if not self._implementation:
            return None
        return self._implementation.read(firstRow, numRows)

    # -----------------------------------------------------------------

    def write(self, firstRow, data):
        """
        NAM: Column.write (Method)
        DES: Writes data to a Column.
             Writing starts at row firstRow. Note that row numbers start at 1.
             For the datatype of data the corresponding rules to Column.read apply:
             For columns that hold numeric data of not-variable length, data must be a
             Numeric array. If the Numeric array's typecode does not match the typecode
             as specified in Column write, data is cast to the correct typecode. (Note that
             this introduces a performance penalty!). The shape of the Numeric array need
             not match the dimension of the Column; however, the number of elements in data
             must be correct.
             For columns with varable length numeric data, data must be a list of Numeric.arrays.
             Each element of the list will be written into a separate row. Concerning
             typecode the above said applies.
             For string data (both of fixed length and variable length), data must be a 
             list of strings.
        INP: firstRow (long): The number of the first Row to which data is written
             data:            See DES
        """
        if self._implementation:
            self._implementation.write(firstRow, data)
            
# ---------------------------------------------------------------------
# ---- Class ColumnInfo -----------------------------------------------
# ---------------------------------------------------------------------

class ColumnInfo:
    """
    NAM: ColumnInfo (Class)
    DES: Contains the information that is needed for the creation of a single column
         during the execution of Dataset.addTable.
         The meaning of the data elements is:
         - name:        The name of the column as in Keyword TTYPEn
         - datatype:    The datatype code for the column as specified in Keyword TFORMn.
                        Must not contain the repeat count nor the length for variable
                        length arrays. The repeat count is specccified in repeat, wheras
                        the length is added automatically when the Table is closed.
                        See Column.read and the CFITSIO documentation for valid datatypes.
         - repeat:      The repeat count as specified in Keyword TFORMn. In case of string
                        data: The maximal length of the string.
         - dim:         The dimension as specified in Keyword TDIMn.
         - description: A descriptive text for the column. Stored as comment of the Keyword
                        TTYPEn.
         - unit:        The unit for the column. Stored as part of the comment of the Keyword
                        TTYPEn.
         If both dim and repeat are specified, the repeat count in Keyword TFORMn is evaluated
         from dim. repeat is ignored in this case.
    """
    def __init__(self, name='', datatype='', repeat=0, dim=[], description='', unit=''):
        self.name = name
        self.datatype = datatype
        self.repeat = repeat
        self.dim = dim
        self.description = description
        self.unit = unit

# ---------------------------------------------------------------------
# ---- Class MBFitsError ----------------------------------------------
# ---------------------------------------------------------------------

class MBFitsError(Exception):
    """
    NAM: MBFitsError (Class)
    DES: Exception class for exceptions related with module BoaMBFits
    """
    def __init__(self, msg, code=None):
        self.code = code
        self.msg = msg
        
    # -----------------------------------------------------------------

    def __str__(self):
        return `self.msg`


# ---------------------------------------------------------------------
# ---- Implementation Classes -----------------------------------------
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# ---- Class _Dataset ----------------------------------------------
# ---------------------------------------------------------------------

class _Dataset:
    # NAM: _Dataset (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def __init__(self, name):
        self._name = name
        
        self._isopen = 0
        self._iomode = 0
        
        self._primaryHeader = None
        self._groupingTable = None
        self._tableHdus = []

    # -----------------------------------------------------------------

    def construct(self, filename, keywords=[], groupName=None):
        # Construct a new Dataset with Primary Header and Grouping Table
        # in the specified file.
        # The Primary Header will contain the specified Keywords.
        # The Grouping Table will contain a group with name groupName.
        # If groupName=None, no Grouping Table will be created.

        if self._primaryHeader:
            raise MBFitsError("Dataset %s exists already" % self._name)

        if os.path.exists(filename):
            raise MBFitsError("File %s exists already" % filename)
        
        self._isopen = 1
        self._iomode = 1
        
        fptr = self._createFile(filename, isDatasetFile=1)
        self._primaryHeader = _Hdu(fptr, self, iomode=1)
        for keyword in keywords:
            self._primaryHeader.addKeyword(keyword)
        self._primaryHeader.exprt()
        
        if groupName:
            self._groupingTable = self._createGroupingTable(filename, groupName)
            self._groupingTable.exprt()

        # Flush the file and reimport to synchronize keywords:
        status = cfitsio.fits_flush_file(fptr)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        self._primaryHeader.imprt()
        if self._groupingTable:
            self._groupingTable.imprt()
            
    # -----------------------------------------------------------------

    def isOpen(self):
        return self._isopen
        
    # -----------------------------------------------------------------

    def isWriteOpen(self):
        return self.isOpen() and self._iomode
        
    # -----------------------------------------------------------------

    def open(self, iomode=0):
        if self.isOpen():
            raise MBFitsError("Dataset %s already open" % self._name)

        self._isopen = 1
        self._iomode = iomode

        if not self._primaryHeader:
            raise MBFitsError("Dataset %s has no Primary Header" % self._name)
        self._primaryHeader.open(iomode)

        if self._groupingTable:
            self._groupingTable.open(iomode)

    # -----------------------------------------------------------------

    def close(self):
        if not self.isOpen():
            raise MBFitsError("Dataset %s not open" % self._name)

        for table in self._tableHdus:
            if table.isOpen():
                table.close()

        if self._groupingTable:
            if self._groupingTable.isOpen():
                self._groupingTable.close()
        
        if self._primaryHeader.isOpen():
            self._primaryHeader.close()
            
        self._isopen = 0
        
    # -----------------------------------------------------------------

    def imprt(self, filename, iomode):
        # Import Dataset from disk. filename specifies the file that contains the Primary Header.
        # On exit, the Dataset will be open, while all Tables will be closed.

        if self._primaryHeader:
            raise MBFitsError("Dataset %s exists already" % self._name)

        # Open the specified file:
        status, fptr = cfitsio.fits_open_file(filename, iomode)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        if not fptr:
            raise MBFitsError("Failed to open file %s" % filename)

        self._isopen = 1
        self._iomode = iomode
        
        status, nHdus = cfitsio.fits_get_num_hdus(fptr)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        # Walk through the specified file and create _TableHdu-objects:
        for iHdu in xrange(1, nHdus+1):
            # Reopen the file to get a new file pointer to be used for this HDU:
            status, fptrHdu = cfitsio.fits_reopen_file(fptr)
            if status:    
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
            if not fptrHdu:
                raise MBFitsError("Failed to reopen file %s for HDU %d" % (filename, iHdu))

            # Move the new file pointer to the HDU: 
            status, cpHduNum = cfitsio.fits_movabs_hdu(fptrHdu, iHdu)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            if iHdu == 1:
                # Primary Header:
                self._primaryHeader = _Hdu(fptrHdu, self, iomode)
                self._primaryHeader.imprt()
                # The primary header is left open:
                #self._primaryHeader.close()
            else:
                # Tables:
                if not _GroupingHdu.isGroupingHdu(fptrHdu):
                    # Create _TableHdu for this Hdu:
                    tableHdu = _TableHdu(fptrHdu, self, iomode)
                    tableHdu.imprt()
                    tableHdu.close()
                    self._tableHdus.append(tableHdu)
                else:
                    if self._groupingTable:
                        raise MBFitsError("More than 1 GROUPING HDU in file %s" % filename)
                    # Create _GroupingHdu for this HDU:
                    self._groupingTable = _GroupingHdu(fptrHdu, self, iomode)
                    # Create _TableHdus for all HDUs referenced in this Grouping HDU:
                    memberTableHdus = self._groupingTable.importMemberTableHdus()
                    self._tableHdus.extend(memberTableHdus)
                    # The Grouping HDU is left open:
                    #self._groupingTable.close()

        # Close the dataset file:
        # Primary Header and Grouping Table remain open!
        status = cfitsio.fits_close_file(fptr)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

    # -----------------------------------------------------------------

    def exprt(self):
        # Export the Dataset, ie. the Primary Header and the Grouping Table.
        # The Dataset's Tables will nopt be exported.

        if not self.isWriteOpen():
            raise MBFitsError("Dataset %s not open to write" % self._name)

        self._primaryHeader.exprt()
        
        if self._groupingTable:
            self._groupingTable.exprt()

    # -----------------------------------------------------------------

    def getName(self):
        return self._name

    # -----------------------------------------------------------------

    def getFilename(self):
        if self._primaryHeader:
            return self._primaryHeader.getFilename()
        else:
            return None

    # -----------------------------------------------------------------

    def getSize(self):
        size = 0
        if self._primaryHeader:
            filenames = [self._primaryHeader.getFilename()]
            if  self._groupingTable:
                filename = self._groupingTable.getFilename()
                if not filename in filenames:
                    filenames.append(filename)
            for tableHdu in self._tableHdus:
                filename = tableHdu.getFilename()
                if not filename in filenames:
                    filenames.append(filename)

            for filename in filenames:
                size += os.path.getsize(filename)

        return size
            
    # -----------------------------------------------------------------

    def getKeywordNames(self):
        if not self._primaryHeader:
            return None
        return self._primaryHeader.getKeywordNames()

    # -----------------------------------------------------------------

    def getKeyword(self, keyname):
        if not self._primaryHeader:
            return None
        return self._primaryHeader.getKeyword(keyname)

    # -----------------------------------------------------------------

    def getTables(self, **kargs):
        if not self.isOpen():
            raise MBFitsError("Dataset %s not open" % self._name)

        if not kargs:
            # Return all tableHdus if no arguments are passed:
            result = self._tableHdus
        else:
            # Return only matching tableHdus if arguments are passed:
            result = []
            for tableHdu in self._tableHdus:
                keywordsMatch = 1
                for karg in kargs.keys():
                    tableKeyword = tableHdu.getKeyword(karg)
                    if kargs[karg] != None:
                        #  keyname=value in kargs
                        if tableKeyword == None:
                            # Keyword keyname not present in tableHDU:
                            keywordsMatch = 0
                            break
                        if tableKeyword.getValue() != kargs[karg]:
                            # Keyword keyname present in tableHDU,
                            # but value does not match:
                            keywordsMatch = 0
                            break
                    else:
                        #  keyname=None in kargs
                        if tableKeyword != None:
                            # Keyword keyname present in tableHDU:
                            keywordsMatch = 0
                            break
                if keywordsMatch:
                    result.append(tableHdu)
        return result

    # -----------------------------------------------------------------

    def addTable(self, filename=None, keywords=[], colinfos=[]):
        # Creates a new Table in the specified file. If filename=None, the Table
        # will be created in the file that contains the Primary Header.
        # The Table will contain the specified Keywords and Columns.
        # The Table will be added to the Grouping Table, if a Grouping Table exists.
        # The Table will be open for reading and writing on exit from this method.

        if not self.isWriteOpen():
            raise MBFitsError("Dataset %s not open to write" % self._name)

        if not filename:
            filename = self._primaryHeader.getFilename()
            fileExists = 1
        else:
            # See if the specified file already exists:
            hdus = [self._primaryHeader, self._groupingTable]
            hdus.extend(self._tableHdus)
            
            fileExists = 0
            for hdu in hdus:
                if hdu:
                    if hdu.getFilename() == filename:
                        fileExists = 1
                        break

        # Create or open the file:
        if not fileExists:
            fptr = self._createFile(filename)
        else:
            iomode = 1 # read-write
            status, fptr = cfitsio.fits_open_file(filename, iomode)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

        # Create the new _TableHdu:
        tableHdu = self._createTable(fptr, keywords, colinfos)
        
        # Add the new _TableHdu to the Grouping Table:
        if self._groupingTable:
            self._groupingTable.addMemberTableHdu(fptr, keywords)
            
        # Add the new _TableHdu to the 
        self._tableHdus.append(tableHdu)
        
        # Export the new _TableHdu:
        tableHdu.exprt()

        # Flush the file and reimport to synchronize keywords:
        status = cfitsio.fits_flush_file(fptr)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        tableHdu.imprt()

        return tableHdu

    # -----------------------------------------------------------------

    def _createFile(self, filename, isDatasetFile=0):
        # Create Directory if necessary:
        (dir, file) = os.path.split(filename)
        if not os.path.isdir(dir):
            os.makedirs(dir, mode=0755)
 
        # Create Fits file:
        status, fptr = cfitsio.fits_create_file(filename)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        # Create primary image in file
        status=cfitsio.fits_create_img(fptr,32,0,[])
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        status, cpHduNum = cfitsio.fits_movabs_hdu(fptr, 1)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
    
        if not isDatasetFile:
            # Add keyword to primary image:
            keyname = "MBFITS"
            keyval = self.getName().split('/')[-1]
            keycom = "Name of MBFitsFile"
            status = cfitsio.fits_update_key_str(fptr, \
                                                 keyname, \
                                                 keyval, \
                                                 keycom)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
          
        return fptr
    
    # -----------------------------------------------------------------

    def _createTable(self, fptr, keywords, colinfos):
        # Get the extension name from the keywords:
        extname = ""
        for keyword in keywords:
            if keyword.getName() == "EXTNAME":
                extname = keyword.getValue()
                break
        
        # From fitsio.h:
        #define IMAGE_HDU  0  /* Primary Array or IMAGE HDU */
        #define ASCII_TBL  1  /* ASCII table HDU  */
        #define BINARY_TBL 2  /* Binary table HDU */
        #define ANY_HDU   -1  /* matches any HDU type */
        binary_tbl = 2
        status = cfitsio.fits_create_tbl(fptr, \
                                         binary_tbl, \
                                         0, \
                                         0, \
                                         [], \
                                         [], \
                                         [], \
                                         extname)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        
        tableHdu = _TableHdu(fptr, self, self._iomode)
        
        for keyword in keywords:
            tableHdu.addKeyword(keyword)
            
        for colinfo in colinfos:
            tableHdu.addColumn(colinfo)
            
        return tableHdu

    # -----------------------------------------------------------------

    def _createGroupingTable(self, filename, groupName):
        status, fptr = cfitsio.fits_open_file(filename, self._iomode)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        
        groupType =_GroupingHdu.getGroupType()
        status = cfitsio.fits_create_group(fptr, groupName, groupType)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        groupingTable = _GroupingHdu(fptr, self, self._iomode)
        return groupingTable

# ---------------------------------------------------------------------
# ---- Class _Hdu -----------------------------------------------------
# ---------------------------------------------------------------------

class _Hdu:
    # NAM: _Hdu (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def __init__(self, fptr, dataset, iomode=0):
        self._dataset = dataset
        self._fptr = fptr
        self._header = _Header(fptr)

        status, filename = cfitsio.fits_file_name(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        if os.path.isfile(filename):
            self._filename = filename
        else:
            gzFilename = filename + ".gz"
            if os.path.isfile(gzFilename):
                self._filename = gzFilename
            else:
                raise MBFitsError("Hdu %s, %d: File does not exist" % (filename, self._hdunum))

        dummy, self._hdunum = cfitsio.fits_get_hdu_num(self._fptr) # Strange output from cfitsio routine

        if iomode:
            if not dataset.isWriteOpen():
                raise MBFitsError("Hdu %s, %d: Dataset not open for write" % (self._filename, self._hdunum))
        self._iomode = iomode

    # -----------------------------------------------------------------

    def isOpen(self):
        if self._fptr:
            return 1
        else:
            return 0
        
    # -----------------------------------------------------------------

    def isWriteOpen(self):
        return self.isOpen() and self._iomode
    
    # -----------------------------------------------------------------

    def open(self, iomode=0):
        if self.isOpen():
            raise MBFitsError("Hdu %s, %d already open" % (self._filename, self._hdunum))

        if not self._dataset.isOpen():
            raise MBFitsError("Hdu %s, %d: Dataset not open" % (self._filename, self._hdunum))

        if iomode:
            if not self._dataset.isWriteOpen():
                raise MBFitsError("Hdu %s, %d: Dataset not open for write" % (self._filename, self._hdunum))
        self._iomode = iomode

        status, fptr = cfitsio.fits_open_file(self._filename, self._iomode)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        if not fptr:
            raise MBFitsError("Failed to open file %s" % self._filename)
       
        status, cpHduNum = cfitsio.fits_movabs_hdu(fptr, self._hdunum)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        self._header.open(fptr)

        self._fptr = fptr
        
    # -----------------------------------------------------------------

    def close(self):
        if not self.isOpen():
            raise MBFitsError("Hdu %s, %d not open" % (self._filename, self._hdunum))

        self._header.close()

        status = cfitsio.fits_close_file(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        self._fptr = None
        
    # -----------------------------------------------------------------

    def imprt(self):
        if not self.isOpen():
            raise MBFitsError("Hdu %s, %d not open" % (self._filename, self._hdunum))

        self._header.imprt()

    # -----------------------------------------------------------------

    def exprt(self):
        if not self.isWriteOpen():
            raise MBFitsError("Hdu %s, %d not open to write" % (self._filename, self._hdunum))

        self._header.exprt()

    # -----------------------------------------------------------------

    def getFilename(self):
        return self._filename

    # -----------------------------------------------------------------

    def getKeywordNames(self):
        return self._header.getKeywordNames()

    # -----------------------------------------------------------------

    def getKeyword(self, key):
        return self._header.getKeyword(key)

    # -----------------------------------------------------------------

    def addKeyword(self, keyword):
        if not self.isWriteOpen():
            raise MBFitsError("Hdu %s, %d not open to write" % (self._filename, self._hdunum))

        self._header.addKeyword(keyword)

# ---------------------------------------------------------------------
# ---- Class _TableHdu ------------------------------------------------
# ---------------------------------------------------------------------

class _TableHdu(_Hdu):
    # NAM: _TableHdu (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def __init__(self, fptr, dataset, iomode=0):
        _Hdu.__init__(self, fptr, dataset, iomode)

        self._columns  = {} # key: colname;   value: _Column-object
        self._colnames = {} # key: colnumber; value: colname

        self._wasExported = 0

        self._selection = None

    # -----------------------------------------------------------------

    def open(self, iomode=0):
        _Hdu.open(self, iomode)

        for column in self._columns.values():
            column.open(self._fptr, iomode)
        
    # -----------------------------------------------------------------

    def close(self):
        if not self.isOpen():
            raise MBFitsError("Table %s, %d not open" % (self._filename, self._hdunum))

        if self.hasSelection():
            self.clearSelection()
            
        for column in self._columns.values():
            column.close()

        _Hdu.close(self)
        
    # -----------------------------------------------------------------

    def imprt(self):
        _Hdu.imprt(self)

        status, nCols = cfitsio.fits_get_num_cols(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        for iCol in xrange(1, nCols+1):
            colinfo = self._importColinfo(iCol)
            col = _Column(self._fptr, iCol, colinfo, iomode=self._iomode)
            self._columns[colinfo.name] = col
            self._colnames[iCol] = colinfo.name

    # -----------------------------------------------------------------

    def exprt(self):
        _Hdu.exprt(self)
        
        if not self._wasExported:
            # Create the columns only during the first exprt:
            nCols = len(self._columns)
            for iCol in xrange(1, nCols+1):
                self._exportColinfo(iCol)

        self._wasExported += 1
        
    # -----------------------------------------------------------------

    def _importColinfo(self, colnumber):
        # Get some parameters here:
        status, name, dmyUnit, datatype, repeat, dmyScale, dmyZero, dmyNulval, dmyTdisp = \
            cfitsio.fits_get_bcolparms(self._fptr, colnumber)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        # Check whether TDIMn keyword exists:
        tdimExists = 1
        status, strTDim, comment = \
            cfitsio.fits_read_keyword(self._fptr, "TDIM"+str(colnumber))
        if status:
            if status == 202:
                # KEY_NO_EXIST
                tdimExists = 0
            else:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

        if tdimExists:
            status, tdim = cfitsio.fits_read_tdim(self._fptr, colnumber, 100)
            if status:
                # For variable length columns, status indicates an error,
                # while tdim is correct anyway. Hence, ignore status here:
                pass

            # Numeric and Fits differ in the ordering of indices:
            tdim.reverse() 
        else:
            tdim = []
            
        # Get unit and description from TTYPEn keyword:
        status, dmyName, description = \
            cfitsio.fits_read_keyword(self._fptr, "TTYPE"+str(colnumber))
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        status, unit = \
            cfitsio.fits_read_key_unit(self._fptr, "TTYPE"+str(colnumber))
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        
        colinfo = ColumnInfo()
        colinfo.name = name
        colinfo.datatype = datatype
        colinfo.repeat = repeat
        colinfo.dim = tdim
        colinfo.unit = unit
        colinfo.description = description
        
        return colinfo
   
    # -----------------------------------------------------------------

    def _exportColinfo(self, colnumber):
        col = self._columns[self._colnames[colnumber]]
        colinfo = col.getColinfo()

        # Add column to table:
        ttype = colinfo.name
        tform = colinfo.datatype
        # Modify tform if necessary:
        if colinfo.dim:
            repeat = 1
            for idim in colinfo.dim:
                repeat *= idim
        elif colinfo.repeat:
            repeat = colinfo.repeat
        else:
            repeat = ''
        tform = str(repeat) + tform
            
        status = cfitsio.fits_insert_col(self._fptr, \
                                         colnumber, \
                                         ttype, \
                                         tform)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        # Set TDIMn keyword if necessary:
        if colinfo.dim:
            tdim = list(colinfo.dim)
             # Numeric and Fits differ in the ordering of indices:
            tdim.reverse()
            if len(tdim)>0:
                strtdim = "("
                for i in tdim:
                    strtdim += str(i) + ","
                strtdim = strtdim[:-1]+")"
            else:
                strtdim = "()"
            status = cfitsio.fits_update_key_str(self._fptr, \
                                                 "TDIM" + str(colnumber), \
                                                 strtdim, \
                                                 "dimension of field")
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
        
        # Modify TTYPEn keyword to hold description and unit:
        if colinfo.description:
            keyname = "TTYPE" + str(colnumber)
            status = cfitsio.fits_update_key_str(self._fptr, \
                                                 keyname, \
                                                 ttype, \
                                                 colinfo.description)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
            
        if colinfo.unit:
            keyname = "TTYPE" + str(colnumber)
            status = cfitsio.fits_write_key_unit(self._fptr, keyname, colinfo.unit)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
    
    # -----------------------------------------------------------------

    def getNumColumns(self):
        return len(self._columns)

    # -----------------------------------------------------------------

    def getColumnNames(self):
        result = []
        # sort according to colnumber:
        for iCol in xrange(1, len(self._colnames)+1):
            result.append(self._colnames[iCol])
        return result
   
    # -----------------------------------------------------------------

    def getColumn(self, colName):
        if not self.isOpen():
            raise MBFitsError("Hdu %s, %d not open" % (self._filename, self._hdunum))

        return self._columns[colName]
   
    # -----------------------------------------------------------------

    def addColumn(self, colinfo):
        if not self.isWriteOpen():
            raise MBFitsError("Hdu %s, %d not open to write" % (self._filename, self._hdunum))

        if colinfo.name in self._columns.keys():
            raise MBFitsError("Keyword %s already present" % colinfo.name)
        
        iCol = len(self._colnames) + 1
        col = _Column(self._fptr, iCol, colinfo, iomode=self._iomode)
        self._columns[colinfo.name] = col
        self._colnames[iCol] = colinfo.name
    
    # -----------------------------------------------------------------

    def getNumRows(self):
        if not self.isOpen():
            raise MBFitsError("Table %s, %d not open" % (self._filename, self._hdunum))

        status, nRows = cfitsio.fits_get_num_rows(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        return nRows

    # -----------------------------------------------------------------

    def getOptimalRowsize(self):
        if not self.isOpen():
            raise MBFitsError("Table %s, %d not open" % (self._filename, self._hdunum))
        
        status, nrows = cfitsio.fits_get_rowsize(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        return nrows

    # -----------------------------------------------------------------

    def setSelection(self, expression, firstRow, numRows):
        if not self.isOpen():
            raise MBFitsError("Table %s, %d not open" % (self._filename, self._hdunum))

        if not numRows:
            numRows = self.getNumRows()

        self._selection = _Selection(self._fptr, expression, firstRow, numRows)
        for column in self._columns.values():
            column.setSelection(self._selection)

    # -----------------------------------------------------------------

    def clearSelection(self):
        for column in self._columns.values():
            column.clearSelection()

        self._selection = None

    # -----------------------------------------------------------------

    def hasSelection(self):
        if self._selection:
            return 1
        else:
            return 0

# ---------------------------------------------------------------------
# ---- Class _GroupingHdu ---------------------------------------------
# ---------------------------------------------------------------------

class _GroupingHdu(_Hdu):
    # NAM: _GroupingHdu (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def isGroupingHdu(cls, fptr):
        """
        Check the Keyword EXTNAME:
        """
        status, extname, comment = \
                cfitsio.fits_read_keyword(fptr, 'EXTNAME')
        if status:
            return 0
        return extname=="'GROUPING'"
    isGroupingHdu = classmethod(isGroupingHdu)

    def getGroupType(cls):
        # Grouptype parameters from fitsio.h:
        #define GT_ID_ALL_URI  0
        #define GT_ID_REF      1
        #define GT_ID_POS      2
        #define GT_ID_ALL      3
        #define GT_ID_REF_URI 11
        #define GT_ID_POS_URI 12
        return 12
    getGroupType = classmethod(getGroupType)

    # -----------------------------------------------------------------

    def __init__(self, fptr, dataset, iomode=0):
        _Hdu.__init__(self, fptr, dataset, iomode)

        # Aditional index columns:
        self._indexColTtypes = ["EXTNAME", \
                                "SUBSNUM", \
                                "FEBE", \
                                "BASEBAND"]
        self._indexColTforms = ["30A", \
                                "1J", \
                                "30A", \
                                "1J"]
        self._indexColWriteFunctions = {
            "EXTNAME": "cfitsio.fits_write_col_str",
            "SUBSNUM": "cfitsio.fits_write_col_lng",
            "FEBE": "cfitsio.fits_write_col_str", 
            "BASEBAND": "cfitsio.fits_write_col_lng"}
        self._indexColBlankValues = {
            "EXTNAME": "", \
            "SUBSNUM": -999,
            "FEBE": "", \
            "BASEBAND": -999}

        self._wasExported = 0

    # -----------------------------------------------------------------

    def exprt(self):
        _Hdu.exprt(self)
        
        if not self._wasExported:
            # Add the additionla index columns only during the first exprt:
            status, nCols = cfitsio.fits_get_num_cols(self._fptr)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
    
            status = cfitsio.fits_insert_cols(self._fptr, \
                                              nCols+1, \
                                              len(self._indexColTtypes), \
                                              self._indexColTtypes, \
                                              self._indexColTforms)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

        self._wasExported += 1

    # -----------------------------------------------------------------

    def importMemberTableHdus(self):
        """
        Impoorts all Tables of the Group
        """
        if not self.isOpen():
            raise MBFitsError("GroupingTable %s, %d not open" % (self._filename, self._hdunum))

        memberTableHdus = []

        status, nGroupMembers = cfitsio.fits_get_num_members(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)
        
        # Open the members and get their file pointers at the same time:
        for iMember in xrange(1, nGroupMembers+1):
            status, fptrMember = \
                    cfitsio.fits_open_member(self._fptr, iMember)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            # Create _TableHdu for this member:
            tableHdu = _TableHdu(fptrMember, self._dataset, self._iomode)
            tableHdu.imprt()
            tableHdu.close()
            memberTableHdus.append(tableHdu)

        return memberTableHdus
    
    # -----------------------------------------------------------------

    def addMemberTableHdu(self, fptrMember, keywords):
        if not self.isWriteOpen():
            raise MBFitsError("Hdu %s, %d not open to write" % (self._filename, self._hdunum))

        dummy, hdunumMember = cfitsio.fits_get_hdu_num(fptrMember) # Strange output from cfitsio routine

        status = cfitsio.fits_add_group_member(self._fptr, fptrMember, hdunumMember)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        status, numrows = cfitsio.fits_get_num_rows(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        # Fill indexing columns:
        keywordDict = {}
        for keyword in keywords:
            keywordDict[keyword.getName()] = keyword
            
        for colName in self._indexColTtypes:
            if colName in keywordDict.keys():
                colValue = keywordDict[colName].getValue()
            else:
                colValue = self._indexColBlankValues[colName]

            # String matching parameters from fitsio.h:
            #define CASESEN   1   /* do case-sensitive string match */
            #define CASEINSEN 0   /* do case-insensitive string match */
            casesen = 1
            status, colnum = cfitsio.fits_get_colnum(self._fptr, casesen, colName)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            args = "(self._fptr, colnum, numrows, 1, 1, [colValue])"
            status = eval(self._indexColWriteFunctions[colName]+args)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
   
# ---------------------------------------------------------------------
# ---- Class _Header --------------------------------------------------
# ---------------------------------------------------------------------

class _Header:
    # NAM: _Header (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def __init__(self, fptr):
        self._fptr = fptr
        
        self._keywords = {} # key: keyword name;   value: Keyword-object
        self._keywordNames = {} # key: keyword number; value: keyword name; 
                                # keyword name may be None. hence the length 
                                # of the two dictionaries may be different!
        
    # -----------------------------------------------------------------

    def isOpen(self):
        if self._fptr:
            return 1
        else:
            return 0
        
    # -----------------------------------------------------------------

    def open(self, fptr):
        self._fptr = fptr

    # -----------------------------------------------------------------

    def close(self):
        if not self.isOpen():
            raise MBFitsError("Hdu %s, %d not open" % (self._filename, self._hdunum))

        self._fptr = None

    # -----------------------------------------------------------------

    def imprt(self):
        if not self.isOpen():
            raise MBFitsError("Hdu %s, %d not open" % (self._filename, self._hdunum))

        status, nKeys, nMoreKeys = cfitsio.fits_get_hdrspace(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        for iKey in xrange(1, nKeys+1):
            # This does not work for "HIERARCH ESO" keywords:
            #status, keyname, value, comment = cfitsio.fits_read_keyn(self._fptr, iKey)
            
            status, record = cfitsio.fits_read_record(self._fptr, iKey)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            if ((record.find('COMMENT') == 0) or (record.find('CONTINUE') == 0)):
                # COMMENT or CONTINUE line
                self._keywordNames[iKey] = None
                continue

            keyname = record.split("=")[0].strip()

            status, value, comment = cfitsio.fits_read_keyword(self._fptr, keyname)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            if len(value) == 0:
                # special case
                datatype = 'A'
                value = ''
            else:
                # Get keytype, translate keytype to datatype, 
                # and convert value to corresponding type:
                status, keytype = cfitsio.fits_get_keytype(value)
                if status:
                    msg = cfitsio.fits_get_errstatus(status)
                    raise MBFitsError("CFITSIO: %s" % msg, status)

                if keytype == 'C':
                    datatype = 'A'
                    value = value[1:-1].strip()
                    if (len(value) > 0):
                        if (value[-1] == '&'):
                            status, value, comment = cfitsio.fits_read_key_longstr(self._fptr, keyname)
                            if status:
                                msg = cfitsio.fits_get_errstatus(status)
                                raise MBFitsError("CFITSIO: %s" % msg, status)
                elif keytype == 'L':
                    datatype = 'L'
                    if value == 'T':
                        value = 1
                    else:
                        value = 0
                elif keytype == 'I':
                    datatype = 'J'
                    value = int(value)
                elif keytype == 'F':
                    datatype = 'D'
                    value = float(value)
                else:
                    raise MBFitsError("Invalid keytype %s for keyword %s" % (keytype, keyname))

            # Get unit:
            status, unit = \
                cfitsio.fits_read_key_unit(self._fptr, keyname)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            # Create keyword:
            kw = Keyword(name=keyname, \
                         value=value, \
                         comment=comment, \
                         datatype=datatype, \
                         unit=unit)

            self._keywords[keyname] = kw
            self._keywordNames[iKey] = keyname
            
    # -----------------------------------------------------------------

    def exprt(self):
        if not self.isOpen():
            raise MBFitsError("Hdu %s, %d not open" % (self._filename, self._hdunum))

        nKeys = len(self._keywordNames)
        for iKey in xrange(1, nKeys+1):
            keywordName = self._keywordNames[iKey]
            if keywordName:
                kw = self._keywords[keywordName]
                if kw:
                    if 'A' in kw.getDatatype():
                        if kw.getFormat():
                            value = kw.getValue()[:kw.getFormat()]
                        else:
                            value = kw.getValue()
                        status = cfitsio.fits_update_key_str(self._fptr, \
                                                             kw.getName(), \
                                                             value, \
                                                             kw.getComment())
                    elif 'L' in kw.getDatatype():
                        status = cfitsio.fits_update_key_log(self._fptr, \
                                                             kw.getName(), \
                                                             kw.getValue(), \
                                                             kw.getComment())
                    elif 'I' in kw.getDatatype():
                        status = cfitsio.fits_update_key_lng(self._fptr, \
                                                             kw.getName(), \
                                                             kw.getValue(), \
                                                             kw.getComment())
                    elif 'J' in kw.getDatatype():
                        status = cfitsio.fits_update_key_lng(self._fptr, \
                                                             kw.getName(), \
                                                             kw.getValue(), \
                                                             kw.getComment())
                    elif 'E' in kw.getDatatype():
                        if kw.getFormat():
                            format = kw.getFormat()
                        else:
                            format = -7
                        status = cfitsio.fits_update_key_flt(self._fptr, \
                                                             kw.getName(), \
                                                             kw.getValue(), \
                                                             format, \
                                                             kw.getComment())
                    elif 'D' in kw.getDatatype():
                        if kw.getFormat():
                            format = kw.getFormat()
                        else:
                            format = -7
                        status = cfitsio.fits_update_key_dbl(self._fptr, \
                                                             kw.getName(), \
                                                             kw.getValue(), \
                                                             format, \
                                                             kw.getComment())
                    else:
                        msg = "Invalid datatype %s in keyword %s" % (kw.getDatatype(), kw.getName())
                        raise MBFitsError(msg)
                    if status:
                        msg = cfitsio.fits_get_errstatus(status)
                        raise MBFitsError("CFITSIO: %s" % msg, status)
               
                    if kw.getUnit():
                        status = cfitsio.fits_write_key_unit(self._fptr, kw.getName() , kw.getUnit())
                        if status:
                            msg = cfitsio.fits_get_errstatus(status)
                            raise MBFitsError("CFITSIO: %s" % msg, status)

    # -----------------------------------------------------------------

    def getKeywordNames(self):
        keynames = []
        
        iKeys = self._keywordNames.keys()
        iKeys.sort()
        for iKey in iKeys:
            keyname = self._keywordNames[iKey]
            if keyname!=None:
                keynames.append(keyname)
        return keynames

    # -----------------------------------------------------------------

    def getKeyword(self, keyname):
        if keyname in self._keywords.keys():
            return self._keywords[keyname]
        else:
            return None
    
    # -----------------------------------------------------------------

    def addKeyword(self, keyword):
        if keyword.getName() in self._keywords.keys():
            raise MBFitsError("Keyword %s already present" % keyword.getName())
        
        iKey = len(self._keywordNames) + 1
        self._keywordNames[iKey] = keyword.getName()
        self._keywords[keyword.getName()] = keyword
    
# ---------------------------------------------------------------------
# ---- Class _Column --------------------------------------------------
# ---------------------------------------------------------------------

class _Column:
    # NAM: _Column (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def __init__(self, fptr, colnum, colinfo, iomode=0):
        self._fptr = fptr
        self._colnum = colnum
        self._colinfo = colinfo
        self._iomode = iomode

        if 'A' in self._colinfo.datatype:
            self._readFunction = "cfitsio.fits_read_col_str"
            self._writeFunction = "cfitsio.fits_write_col_str"
            self._numericTypecode = None
        elif 'L' in self._colinfo.datatype:
            self._readFunction = "cfitsio.fits_read_col_log"
            self._writeFunction = "cfitsio.fits_write_col_log"
            self._numericTypecode = Numeric.Int0
        elif 'B' in self._colinfo.datatype:
            # 'B'
            self._readFunction = "cfitsio.fits_read_col_byt"
            self._writeFunction = "cfitsio.fits_write_col_byt"
            self._numericTypecode = Numeric.UnsignedInt8
        elif ('I' in self._colinfo.datatype): 
            # 'I', 'UI'
            self._readFunction = "cfitsio.fits_read_col_int"
            self._writeFunction = "cfitsio.fits_write_col_int"
            self._numericTypecode = Numeric.Int16
        elif ('K' in self._colinfo.datatype): 
            # 'K', 'UK'
            self._readFunction = "cfitsio.fits_read_col_lng"
            self._writeFunction = "cfitsio.fits_write_col_lng"
            self._numericTypecode = Numeric.Int32
        elif ('J' in self._colinfo.datatype):
            # 'J', 'JJ', 'UJ'
            self._readFunction = "cfitsio.fits_read_col_lng"
            self._writeFunction = "cfitsio.fits_write_col_lng"
            self._numericTypecode = Numeric.Int32
        elif ('E' in self._colinfo.datatype):
            self._readFunction = "cfitsio.fits_read_col_flt"
            self._writeFunction = "cfitsio.fits_write_col_flt"
            self._numericTypecode = Numeric.Float32
        elif ('F' in self._colinfo.datatype):
            if ('FC' in self._colinfo.datatype):
                self._readFunction = None
                self._writeFunction = None
                self._numericTypecode = None
            elif ('FM' in self._colinfo.datatype):
                self._readFunction = None
                self._writeFunction = None
                self._numericTypecode = None
            else:
                self._readFunction = "cfitsio.fits_read_col_flt"
                self._writeFunction = "cfitsio.fits_write_col_flt"
                self._numericTypecode = Numeric.Float32
        elif ('D' in self._colinfo.datatype):
            self._readFunction = "cfitsio.fits_read_col_dbl"
            self._writeFunction = "cfitsio.fits_write_col_dbl"
            self._numericTypecode = Numeric.Float64
        elif ('G' in self._colinfo.datatype):
            self._readFunction = "cfitsio.fits_read_col_dbl"
            self._writeFunction = "cfitsio.fits_write_col_dbl"
            self._numericTypecode = Numeric.Float64
        else:
            self._readFunction = None
            self._writeFunction = None
            self._numericTypecode = None

        self._selection = None

    # -----------------------------------------------------------------

    def isOpen(self):
        if self._fptr:
            return 1
        else:
            return 0
        
    # -----------------------------------------------------------------

    def isWriteOpen(self):
        return self.isOpen() and self._iomode
    
    # -----------------------------------------------------------------

    def open(self, fptr, iomode=0):
        self._fptr = fptr
        self._iomode = iomode

    # -----------------------------------------------------------------

    def close(self):
        if not self.isOpen():
            raise MBFitsError("Column %d not open" % self._colnum)

        self._fptr = None

    # -----------------------------------------------------------------

    def getColnum(self):
        return self._colnum
    
    # -----------------------------------------------------------------

    def getColinfo(self):
        return self._colinfo
    
    # -----------------------------------------------------------------

    def getName(self):
        return self._colinfo.name
    
    # -----------------------------------------------------------------

    def getDatatype(self):
        return self._colinfo.datatype
    
    # -----------------------------------------------------------------

    def getRepeat(self):
        return self._colinfo.repeat

    # -----------------------------------------------------------------

    def getDim(self):
         return self._colinfo.dim

    # -----------------------------------------------------------------

    def getDescription(self):
        return self._colinfo.description

    # -----------------------------------------------------------------

    def getUnit(self):
        return self._colinfo.unit
    
    # -----------------------------------------------------------------

    def setSelection(self, selection):
        if not self.isOpen():
            raise MBFitsError("Column %d not open" % self._colnum)

        self._selection = selection

    # -----------------------------------------------------------------

    def clearSelection(self):
        self._selection = None

    # -----------------------------------------------------------------

    def hasSelection(self):
        if self._selection:
            return 1
        else:
            return 0

    # -----------------------------------------------------------------

    def read(self, firstRow, numRows):
        if not self.isOpen():
            raise MBFitsError("Column %s not open" % self._colinfo.name)

        if not self._readFunction:
            return None

        status, numRowsPresent = cfitsio.fits_get_num_rows(self._fptr)
        
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        if not numRows:
            numRows = numRowsPresent

        if firstRow + numRows  > numRowsPresent + 1:
            numRows = numRowsPresent + 1 - firstRow

        if self._selection:
            selectedRows = self._selection.getSelectedRows()
            if 'A' in self._colinfo.datatype:
                return self._readSelectionA(firstRow, numRows, selectedRows)
            else:
                if 'P' in self._colinfo.datatype:
                    return self._readSelectionP(firstRow, numRows, selectedRows)
                else:
                    return self._readSelection(firstRow, numRows, selectedRows)
        else:
            if 'P' in self._colinfo.datatype:
                if 'A' in self._colinfo.datatype:
                    return self._readPA(firstRow, numRows)
                else:
                    return self._readP(firstRow, numRows)
            else:
                if 'A' in self._colinfo.datatype:
                    return self._readA(firstRow, numRows)
                else:
                    return self._read(firstRow, numRows)
            
    
    # -----------------------------------------------------------------

    def _read(self, firstRow, numRows):
        """
        Returns Numeric.array
        """
        if self._colinfo.dim:
            repeat = 1
            for idim in self._colinfo.dim:
                repeat *= idim
        elif self._colinfo.repeat:
            repeat = self._colinfo.repeat
        else:
            repeat = 1

        args = "(self._fptr, self._colnum, firstRow, 1, numRows*repeat, 0)"
        status, result = eval(self._readFunction+args)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        result = self._toNumericArray(result)


        shape = [numRows]
        if self._colinfo.dim:
            shape.extend(self._colinfo.dim)
        else:
            if self._colinfo.repeat > 1:
                shape.extend([self._colinfo.repeat])
        result.shape = shape

        return result

    # -----------------------------------------------------------------

    def _readA(self, firstRow, numRows):
        """
        Returns list of strings
        """
        args = "(self._fptr, self._colnum, firstRow, 1, numRows, 0)"
        status, result = eval(self._readFunction+args)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        if not result:
            result = []

        if type(result) != type([]):
            result = [result]

        return result

    # -----------------------------------------------------------------

    def _readP(self, firstRow, numRows):
        """
        Returns list of Numeric.arrays
        """
        result = []
        for iRow in xrange(firstRow, firstRow+numRows):
            status, rowRepeat, offset = cfitsio.fits_read_descript(self._fptr, \
                                                                   self._colnum, \
                                                                   iRow)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            args = "(self._fptr, self._colnum, iRow, 1, rowRepeat, 0)"
            status, rowResult = eval(self._readFunction+args)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            rowResult = self._toNumericArray(rowResult)
            result.append(rowResult)

        return result

    # -----------------------------------------------------------------

    def _readPA(self, firstRow, numRows):
        """
        Returns list of strings.
        """
        result = []
        for iRow in xrange(firstRow, firstRow+numRows):
            args = "(self._fptr, self._colnum, iRow, 1, 1, 0)"
            status, rowResult = eval(self._readFunction+args)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            result.append(rowResult)

        return result

    # -----------------------------------------------------------------

    def _readSelection(self, firstRow, numRows, selectedRows):
        """
        Returns list of Numeric.arrays
        """
        result = []

        if self._colinfo.dim:
            repeat = 1
            for idim in self._colinfo.dim:
                repeat *= idim
        elif self._colinfo.repeat:
            repeat = self._colinfo.repeat
        else:
            repeat = 1

        for iRow in xrange(firstRow, firstRow+numRows):
            if iRow in selectedRows:
                args = "(self._fptr, self._colnum, iRow, 1, repeat, 0)"
                status, rowResult = eval(self._readFunction+args)
                if status:
                    msg = cfitsio.fits_get_errstatus(status)
                    raise MBFitsError("CFITSIO: %s" % msg, status)

                rowResult = self._toNumericArray(rowResult)
                result.append(rowResult)

        shape = [len(result)]
        if self._colinfo.dim:
            shape.extend(self._colinfo.dim)
        else:
            if self._colinfo.repeat > 1:
                shape.extend([self._colinfo.repeat])
        result = Numeric.array(result)
        result.shape = shape

        return result

    # -----------------------------------------------------------------

    def _readSelectionA(self, firstRow, numRows, selectedRows):
        """
        Returns list of strings
        """
        result = []
        for iRow in xrange(firstRow, firstRow+numRows):
            if iRow in selectedRows:
                args = "(self._fptr, self._colnum, iRow, 1, 1, 0)"
                status, rowResult = eval(self._readFunction+args)
                if status:
                    msg = cfitsio.fits_get_errstatus(status)
                    raise MBFitsError("CFITSIO: %s" % msg, status)

                result.append(rowResult)

        return result

    # -----------------------------------------------------------------

    def _readSelectionP(self, firstRow, numRows, selectedRows):
        """
        Returns list of Numeric.arrays
        """
        result = []
        for iRow in xrange(firstRow, firstRow+numRows):
            if iRow in selectedRows:
                status, rowRepeat, offset = cfitsio.fits_read_descript(self._fptr, \
                                                                       self._colnum, \
                                                                       iRow)
                if status:
                    msg = cfitsio.fits_get_errstatus(status)
                    raise MBFitsError("CFITSIO: %s" % msg, status)

                args = "(self._fptr, self._colnum, iRow, 1, rowRepeat, 0)"
                status, rowResult = eval(self._readFunction+args)
                if status:
                    msg = cfitsio.fits_get_errstatus(status)
                    raise MBFitsError("CFITSIO: %s" % msg, status)

                rowResult = self._toNumericArray(rowResult)
                result.append(rowResult)

        return result

    # -----------------------------------------------------------------

    def write(self, firstRow, data):
        if not self.isWriteOpen():
            raise MBFitsError("Column %s not open to write" % self._colinfo.name)

        if self._writeFunction:
            if 'P' in self._colinfo.datatype:
                if 'A' in self._colinfo.datatype:
                    self._writePA(firstRow, data)
                else:
                    self._writeP(firstRow, self._toNumericArray(data))
            else:
                if 'A' in self._colinfo.datatype:
                    self._writeA(firstRow, data)
                else:
                    self._write(firstRow, self._toNumericArray(data))
            
    # -----------------------------------------------------------------

    def _write(self, firstRow, data):
        """
        data must be Numeric.array
        """
        # Cast to proper type if necessary:
        if data.typecode() != self._numericTypecode:
            coldata = data.astype(self._numericTypecode)
        else:
            coldata = data

        args = "(self._fptr, self._colnum, firstRow, 1, Numeric.size(coldata), coldata)"
        status = eval(self._writeFunction+args)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

    # -----------------------------------------------------------------

    def _writeA(self, firstRow, data):
        """
        data must be list of strings
        """
        args = "(self._fptr, self._colnum, firstRow, 1, len(data), data)"
        status = eval(self._writeFunction+args)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

    # -----------------------------------------------------------------

    def _writeP(self, firstRow, data):
        """
        data must be list of Numeric.arrays
        """
        rownum = firstRow
        for dataelement in data:
            # Cast to proper type if necessary:
            if dataelement.typecode() != self._numericTypecode:
                rowdata = dataelement.astype(self._numericTypecode)
            else:
                rowdata = dataelement
                
            args = "(self._fptr, self._colnum, rownum, 1, len(rowdata), rowdata)"
            status = eval(self._writeFunction+args)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)
            
            rownum += 1

    # -----------------------------------------------------------------

    def _writePA(self, firstRow, data):
        """
        data must be list of strings
        """
        rownum = firstRow
        for rowdata in data:
            args = "(self._fptr, self._colnum, rownum, 1, len(rowdata), rowdata)"
            status = eval(self._writeFunction+args)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise MBFitsError("CFITSIO: %s" % msg, status)

            rownum += 1

    # -----------------------------------------------------------------

    def _toNumericArray(self, x):
        """
        Convert x to a Numeric array, if possible.
        """
        if type(x) == type(Numeric.array([])):
            result = x
        else:
            try:
                # This works for lists and tuples, but raises an exception for scalars
                result = Numeric.array(list(x))
            except:
                # This works for scalars:
                result = Numeric.array([x])
        return result

# ---------------------------------------------------------------------
# ---- Class _Selection -----------------------------------------------
# ---------------------------------------------------------------------

class _Selection:
    # NAM: _Selection (Class)
    # DES: Implementation class; do not use from outside of module BoaMBFits!

    def __init__(self, fptr, expression, firstRow, numRows):
        # self._selectedRows is list of row numbers which conform to the 
        # boolean expression
        status, numRowsFound, rowBooleanResults = \
            cfitsio.fits_find_rows(fptr, expression, firstRow, numRows)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise MBFitsError("CFITSIO: %s" % msg, status)

        self._selectedRows = \
            [i+1 for i in xrange(len(rowBooleanResults)) \
             if rowBooleanResults[i] != '0']

    # -----------------------------------------------------------------

    def getSelectedRows(self):
        return self._selectedRows
    
