
# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------

import os
import Numeric
import cfitsio

# ---------------------------------------------------------------------
# ---- Module Interface -----------------------------------------------
# ---------------------------------------------------------------------
"""
The Module Interface contains:

- Function createDataset
- Class Dataset
- Class Image
- Class Keyword
- Class FitsError

Only these functions and classes should be used by clients!
"""

# ---------------------------------------------------------------------
# ---- Function createDataset -----------------------------------------
# ---------------------------------------------------------------------

def createDataset(datasetName):
    """
    Creates a new Dataset object and create new disk files specified by
    datasetName.
    """
    datasetImpl = _Dataset()
    datasetImpl.create(datasetName)
    return Dataset(datasetImpl)

# ---------------------------------------------------------------------
# ---- Class Dataset --------------------------------------------------
# ---------------------------------------------------------------------

class Dataset:
    
    def __init__(self, implementation):
        self._implementation = implementation

    # -----------------------------------------------------------------

    def close(self):
        if self._implementation:
            self._implementation.close()

    # -----------------------------------------------------------------

    def isOpen(self):
        if not self._implementation:
            return 0
        return self._implementation.isOpen()

    # -----------------------------------------------------------------

    def getName(self):
        if not self._implementation:
            return None
        return self._implementation.getName()

    # -----------------------------------------------------------------

    def createImage(self, imageData):
        if not self._implementation:
            return None
        
        imageImpl = self._implementation.createImage(imageData)
        return Image(imageImpl)    
        
    # -----------------------------------------------------------------

    def getImages(self, **keywords):
        if not self._implementation:
            return None
        imageImpls = self._implementation.getImages(**keywords)
        images = []
        for imageImpl in imageImpls:
            if imageImpl:
                images.append(Image(imageImpl))
        return images

# ---------------------------------------------------------------------
# ---- Class Image ----------------------------------------------------
# ---------------------------------------------------------------------

class Image:

    def __init__(self, implementation):
        self._implementation = implementation

    # -----------------------------------------------------------------

    def createKeyword(self, keyname, value, unit="", comment=""):
        if not self._implementation:
            return None
        keywordImpl = self._implementation.createKeyword(keyname, \
                                                         value, \
                                                         unit, \
                                                         comment)
        if not keywordImpl:
            return None
        keyword = Keyword(keywordImpl)
        return keyword

    # -----------------------------------------------------------------

    def createKeywordDate(self):
        if not self._implementation:
            return None
        keywordImpl = self._implementation.createKeywordDate()
        if not keywordImpl:
            return None
        keyword = Keyword(keywordImpl)
        return keyword

    # -----------------------------------------------------------------

    def getShape(self):
        if self._implementation:
            return self._implementation.getShape()

    # -----------------------------------------------------------------

    def writeImage(self):
        if self._implementation:
            self._implementation.writeImage()
    
# ---------------------------------------------------------------------
# ---- Class Keyword --------------------------------------------------
# ---------------------------------------------------------------------

class Keyword:
    
    def __init__(self, implementation):
        self._implementation = implementation

    # -----------------------------------------------------------------

    def __str__(self):
        return str(self._implementation)
        
    # -----------------------------------------------------------------

    def getKeyname(self):
        if not self._implementation:
            return None
        return self._implementation.getKeyname()

    # -----------------------------------------------------------------

    def getDatatype(self):
        if not self._implementation:
            return None
        return self._implementation.getDatatype()

    # -----------------------------------------------------------------

    def getValue(self):
        if not self._implementation:
            return None
        return self._implementation.getValue()

    # -----------------------------------------------------------------

    def getUnit(self):
        if not self._implementation:
            return None
        return self._implementation.getUnit()

    # -----------------------------------------------------------------

    def getComment(self):
        if not self._implementation:
            return None
        return self._implementation.getComment()

# ---------------------------------------------------------------------
# ---- Class FitsError ------------------------------------------------
# ---------------------------------------------------------------------

class FitsError(Exception):

    def __init__(self, msg, code=None):
        self.code = code
        self.msg = msg
        
    # -----------------------------------------------------------------

    def __str__(self):
        return `self.msg`


# ---------------------------------------------------------------------
# ---- Implementation Classes -----------------------------------------
# ---------------------------------------------------------------------
"""
The rest of the module contains classes that implement the module's
functionality.

These classes should not be used directly from outside the module,
but only through the module's interface!
"""

# ---------------------------------------------------------------------
# ---- Class _Dataset ----------------------------------------------
# ---------------------------------------------------------------------

class _Dataset:

    def __init__(self):
        self._name = None
        self._fptr = None
        self._hdus = []

    # -----------------------------------------------------------------

    def create(self, filename):
        if self._fptr:
            raise FitsError("File %s already open" % filename)

        self._name = filename

        status, self._fptr = cfitsio.fits_create_file(filename)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)
        if not self._fptr:
            raise FitsError("Failed to create file %s" % filename)

    # -----------------------------------------------------------------

    def close(self):
        if not self._fptr:
            raise FitsError("Dataset is not open")

        if self._hdus:
            # If no HDU exists, fits_close_file will complain
            # because a valid FITS file has to contain at least
            # the primary HDU.
            # Hence, we cheat here and simply skip fits_close_file.
            for hdu in self._hdus:
                hdu.close()
            self._hdus = []

            status = cfitsio.fits_close_file(self._fptr)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise FitsError("CFITSIO: %s" % msg, status)

        self._fptr = None
        self._name = None


    # -----------------------------------------------------------------

    def isOpen(self):
        if self._fptr:
            isOpen = 1
        else:
            isOpen = 0
        return isOpen

    # -----------------------------------------------------------------

    def getName(self):
        if not self._fptr:
            raise FitsError("Dataset is not open")

        return self._name

    # -----------------------------------------------------------------

    def createImage(self, imageData):
        # Reopen the file to get a new file pointer to be used for the new HDU:
        status, fptrHdu = cfitsio.fits_reopen_file(self._fptr)
        if status:    
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)
        if not fptrHdu:
            raise FitsError("Failed to reopen file %s" % self._name)

        imageHdu = _ImageHdu(fptrHdu)
        imageHdu.create(imageData)

        self._hdus.append(imageHdu)
        return imageHdu
        
    # -----------------------------------------------------------------

    def getImages(self, **kargs):
        result = []
        for hdu in self._hdus:
            keywordsMatch = 1
            for karg in kargs.keys():
                if kargs[karg]:
                    # Arguments that are None are ignored
                    hduKeyword = hdu.getKeyword(karg)
                    if not hduKeyword:
                        keywordsMatch = 0
                        break
                    if hduKeyword.getValue() != kargs[karg]:
                        keywordsMatch = 0
                        break
            if keywordsMatch:
                result.append(hdu)
        return result



# ---------------------------------------------------------------------
# ---- Class _Hdu -----------------------------------------------------
# ---------------------------------------------------------------------

class _Hdu:

    def __init__(self, fptr):
        self._fptr = fptr
        self._header = _Header(fptr)

    # -----------------------------------------------------------------

    def close(self):
        status = cfitsio.fits_close_file(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)

        self._fptr = None
        self._header = None
        
    # -----------------------------------------------------------------

    def createKeyword(self, keyname, value, unit="", comment=""):
        return self._header.createKeyword(keyname, value, unit, comment)

    # -----------------------------------------------------------------

    def createKeywordDate(self):
        return self._header.createKeywordDate()

# ---------------------------------------------------------------------
# ---- Class _ImageHdu ------------------------------------------------
# ---------------------------------------------------------------------

class _ImageHdu(_Hdu):

    def __init__(self, fptr):
        _Hdu.__init__(self, fptr)

        self._imageData = None
        self._writeFunction = None

    # -----------------------------------------------------------------

    def create(self, imageData):
        self._imageData = imageData

        typecode = imageData.typecode()
        if (typecode == Numeric.Float0  or \
            typecode == Numeric.Float8  or \
            typecode == Numeric.Float16 or \
            typecode == Numeric.Float32):
            bitpix = -32
            self._writeFunction = cfitsio.fits_write_img_flt
        elif (typecode == Numeric.Float or \
              typecode == Numeric.Float64):
            bitpix = -64
            self._writeFunction = cfitsio.fits_write_img_dbl
        elif (typecode == Numeric.Int0 or \
              typecode == Numeric.Int8):
            bitpix = 8
            self._writeFunction = cfitsio.fits_write_img_sht
        elif (typecode == Numeric.Int16):
            bitpix = 16
            self._writeFunction = cfitsio.fits_write_img_int
        elif (typecode == Numeric.Int or \
              typecode == Numeric.Int32):
            bitpix = 32
            self._writeFunction = cfitsio.fits_write_img_lng
        else:
            msg = "Invalid typecode '%s' for image data" % imageData.typecode()
            raise FitsError(msg)

        naxes = list(imageData.shape)
        naxis = len(naxes)

        status = cfitsio.fits_create_img(self._fptr, \
                                         bitpix, \
                                         naxis, \
                                         naxes)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)
        
    # -----------------------------------------------------------------

    def getShape(self):
        if self._imageData:
            result = self._imageData.shape
        else:
            result = tuple([])
        return result
    
    # -----------------------------------------------------------------

    def writeImage(self):
        if not self._writeFunction:
            msg = "No writeFunction defined"
            raise FitsError(msg)

        imageData1d = Numeric.ravel(Numeric.transpose(self._imageData))
        status = self._writeFunction(self._fptr, \
                                     0, \
                                     1, \
                                     len(imageData1d), \
                                     imageData1d)
        imageData1d = None
        
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)
        
    
        
# ---------------------------------------------------------------------
# ---- Class _Header --------------------------------------------------
# ---------------------------------------------------------------------

class _Header:

    def __init__(self, fptr):
        self._fptr = fptr
        self._keywords = {}
        
    # -----------------------------------------------------------------

    def createKeyword(self, keyname, value, unit="", comment=""):
        if keyname in self._keywords.keys():
            raise FitsError("Keyword %s already present" % keyname)

        # Determination of datatype incomplete!
        # Boolean (datatype = 14) is missing!
        if type(value) == type(""):
            status = cfitsio.fits_write_key_str(self._fptr, \
                                                keyname, \
                                                value, \
                                                comment)
            datatype = 'C'
        elif type(value) == type(0):
            status = cfitsio.fits_write_key_lng(self._fptr, \
                                                keyname, \
                                                value, \
                                                comment)
            datatype = 'I'
        elif type(value) == type(0.):
            decimals = -7
            status = cfitsio.fits_write_key_dbl(self._fptr, \
                                                keyname, \
                                                value, \
                                                decimals, \
                                                comment)
            datatype = 'F'
        else:
            msg = "Cannot determine datatype for value %s" % str(value)
            raise FitsError(msg)

        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)

        if unit:
            status = cfitsio.fits_write_key_unit(self._fptr, \
                                                 keyname, \
                                                 unit)
            if status:
                msg = cfitsio.fits_get_errstatus(status)
                raise FitsError("CFITSIO: %s" % msg, status)

        self._keywords[keyname] = _Keyword(keyname, value, unit, comment, datatype)
        return self._keywords[keyname]
    
    # -----------------------------------------------------------------

    def createKeywordDate(self):
        keyname = "DATE"

        if keyname in self._keywords.keys():
            raise FitsError("Keyword %s already present" % keyname)

        status = cfitsio.fits_write_date(self._fptr)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)

        self._keywords[keyname] = self._readKeyword(keyname)
        return self._keywords[keyname]
    
    # -----------------------------------------------------------------

    def _readKeyword(self, keyname):
        status, value, comment = cfitsio.fits_read_keyword(self._fptr, keyname)
        if status:
            if status == 202:
                # KEY_NO_EXIST
                return None
            else:
                msg = cfitsio.fits_get_errstatus(status)
                raise FitsError("CFITSIO: %s" % msg, status)

        status, unit = cfitsio.fits_read_key_unit(self._fptr, keyname)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)

        status, datatype = cfitsio.fits_get_keytype(value)
        if status:
            msg = cfitsio.fits_get_errstatus(status)
            raise FitsError("CFITSIO: %s" % msg, status)

        # Rudimentary datatype conversion:
        if datatype == 'C':
            value = value[1:-1].strip()
        elif datatype == 'I':
            value = int(value)
        elif datatype == 'F':
            value = float(value)

        return _Keyword(keyname, value, unit, comment, datatype)
    
# ---------------------------------------------------------------------
# ---- Class _Keyword -------------------------------------------------
# ---------------------------------------------------------------------

class _Keyword:
    
    def __init__(self, keyname, value, unit, comment, datatype):
        self._keyname = keyname
        self._value = value
        self._unit = unit
        self._comment = comment
        self._datatype = datatype

    # -----------------------------------------------------------------

    def __str__(self):
        return str((self._keyname, \
                    self._value, \
                    self._comment, \
                    self._datatype))
        
    # -----------------------------------------------------------------

    def getKeyname(self):
        return self._keyname

    # -----------------------------------------------------------------

    def getDatatype(self):
        return self._datatype

    # -----------------------------------------------------------------

    def getValue(self):
        return self._value

    # -----------------------------------------------------------------

    def getUnit(self):
        return self._unit

    # -----------------------------------------------------------------

    def getComment(self):
        return self._comment


