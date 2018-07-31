# ==========================================================================
# This module performs unit tests for the GammaLib FITS module.
#
# Copyright (C) 2012-2018 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import os
import math
import gammalib
import test_support


# =================================== #
# Test class for GammaLib FITS module #
# =================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib FITS module
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Setup GFits container
    def _setup_fits(self):
        """
        Setup GFits container

        Returns
        -------
        fits : `~gammalib.GFits`
            GFits container
        """
        # Setup GFits container
        fits  = gammalib.GFits()
        image = gammalib.GFitsImageDouble()
        for i in range(10):
            image.extname('%s' % i)
            fits.append(image)

        # Return GFits container
        return fits

    # Setup GFitsHeader container
    def _setup_header(self):
        """
        Setup GFitsHeader container

        Returns
        -------
        header : `~gammalib.GFitsHeader`
            GFitsHeader container
        """
        # Setup GFitsHeader container
        header = gammalib.GFitsHeader()
        card   = gammalib.GFitsHeaderCard()
        for i in range(10):
            card.keyname('%s' % i)
            header.append(card)

        # Return GFitsHeader container
        return header

    # Test GFits class access operators
    def _test_fits_access(self):
        """
        Test GFits class parameter access
        """
        # Setup GFits container and element
        fits  = self._setup_fits()
        image = gammalib.GFitsImageDouble()

        # Loop over all elements using the container iterator and count the
        # number of iterations
        niter = 0
        for hdu in fits:
            self.test_value(hdu.extname(), '%s' % niter)
            niter += 1

        # Check that looping was successful
        self.test_value(fits.size(), 10, 'Check container size')
        self.test_value(niter, 10, 'Check container iterator')

        # Test access from start
        self.test_value(fits[3].extname(), '3')

        # Check access from end
        self.test_value(fits[-2].extname(), '8')

        # Check parameter setting by index from start
        image.extname('98')
        fits[3] = image
        self.test_value(fits[3].extname(), '98')

        # Check parameter setting by index from end
        image.extname('99')
        fits[-2] = image
        self.test_value(fits[-2].extname(), '99')

        # Return
        return

    # Test GFits class slicing
    def _test_fits_slicing(self):
        """
        Test GFits class slicing
        """
        # Setup XML container
        fits = self._setup_fits()

        # Test fits[start:end]
        self.test_value(len(fits[3:5]), 2)
        self.test_value(fits[3:5][0].extname(), '3')
        self.test_value(fits[3:5][1].extname(), '4')

        # Test fits[start:]
        self.test_value(len(fits[7:]), 3)
        self.test_value(fits[7:][0].extname(), '7')
        self.test_value(fits[7:][1].extname(), '8')
        self.test_value(fits[7:][2].extname(), '9')

        # Test fits[:end]
        self.test_value(len(fits[:2]), 2)
        self.test_value(fits[:2][0].extname(), '0')
        self.test_value(fits[:2][1].extname(), '1')

        # Test fits[:]
        self.test_value(len(fits[:]), 10)
        for i in range(10):
            self.test_value(fits[:][i].extname(), '%s' % i)

        # Test fits[start:end:step]
        self.test_value(len(fits[3:7:2]), 2)
        self.test_value(fits[3:7:2][0].extname(), '3')
        self.test_value(fits[3:7:2][1].extname(), '5')

        # Test fits[start:end:step]
        self.test_value(len(fits[6:3:-2]), 2)
        self.test_value(fits[6:3:-2][0].extname(), '6')
        self.test_value(fits[6:3:-2][1].extname(), '4')

        # Test fits[-start:]
        self.test_value(len(fits[-2:]), 2)
        self.test_value(fits[-2:][0].extname(), '8')
        self.test_value(fits[-2:][1].extname(), '9')

        # Test fits[:-end]
        self.test_value(len(fits[:-7]), 3)
        self.test_value(fits[:-7][0].extname(), '0')
        self.test_value(fits[:-7][1].extname(), '1')
        self.test_value(fits[:-7][2].extname(), '2')

        # Return
        return

    # Test GFitsHeader class access operators
    def _test_header_access(self):
        """
        Test GFitsHeader class parameter access
        """
        # Setup FITS header
        header = self._setup_header()

        # Loop over all elements using the container iterator and count the
        # number of iterations
        niter = 0
        for card in header:
            self.test_value(card.keyname(), '%s' % niter)
            niter += 1

        # Check that looping was successful
        self.test_value(header.size(), 10, 'Check container size')
        self.test_value(niter, 10, 'Check container iterator')

        # Test access from start
        self.test_value(header[3].keyname(), '3')

        # Check access from end
        self.test_value(header[-2].keyname(), '8')

        # Setup header card
        card = gammalib.GFitsHeaderCard()

        # Check parameter setting by index from start
        card.keyname('98')
        header[3] = card
        self.test_value(header[3].keyname(), '98')

        # Check parameter setting by index from end
        card.keyname('99')
        header[-2] = card
        self.test_value(header[-2].keyname(), '99')

        # Return
        return

    # Test GFitsHeader class slicing
    def _test_header_slicing(self):
        """
        Test GFitsHeader class slicing
        """
        # Setup FITS header
        header = self._setup_header()

        # Test header[start:end]
        self.test_value(len(header[3:5]), 2)
        self.test_value(header[3:5][0].keyname(), '3')
        self.test_value(header[3:5][1].keyname(), '4')

        # Test header[start:]
        self.test_value(len(header[7:]), 3)
        self.test_value(header[7:][0].keyname(), '7')
        self.test_value(header[7:][1].keyname(), '8')
        self.test_value(header[7:][2].keyname(), '9')

        # Test header[:end]
        self.test_value(len(header[:2]), 2)
        self.test_value(header[:2][0].keyname(), '0')
        self.test_value(header[:2][1].keyname(), '1')

        # Test header[:]
        self.test_value(len(header[:]), 10)
        for i in range(10):
            self.test_value(header[:][i].keyname(), '%s' % i)

        # Test header[start:end:step]
        self.test_value(len(header[3:7:2]), 2)
        self.test_value(header[3:7:2][0].keyname(), '3')
        self.test_value(header[3:7:2][1].keyname(), '5')

        # Test header[start:end:step]
        self.test_value(len(header[6:3:-2]), 2)
        self.test_value(header[6:3:-2][0].keyname(), '6')
        self.test_value(header[6:3:-2][1].keyname(), '4')

        # Test header[-start:]
        self.test_value(len(header[-2:]), 2)
        self.test_value(header[-2:][0].keyname(), '8')
        self.test_value(header[-2:][1].keyname(), '9')

        # Test header[:-end]
        self.test_value(len(header[:-7]), 3)
        self.test_value(header[:-7][0].keyname(), '0')
        self.test_value(header[:-7][1].keyname(), '1')
        self.test_value(header[:-7][2].keyname(), '2')

        # Return
        return

    # Test GFile class interface
    def _test_fits(self):
        """
        Test GFile class interface
        """
        # Get test data directory
        datadir = os.environ['TEST_DATA'] + '/'
        
        # Set test file names
        file = gammalib.GFilename(datadir+'file.fits')

        # Test loading of FITS file
        self.test_try('Test GFits file load constructor')
        try:
            fits = gammalib.GFits(file)
            self.test_try_success()
        except:
            self.test_try_failure('Unable to load FITS file.')

        # Test creation of FITS file
        self.test_try('Test GFits file creation constructor')
        try:
            fits = gammalib.GFits('test_file.fits', True)
            self.test_try_success()
        except:
            self.test_try_failure('Unable to create FITS file.')

        # Open FITS file
        fits = gammalib.GFits(file)
        self.test_value(fits['EVENTS'].nrows(), 1231, 'Check number of rows')

        # Open FITS file with event selection
        fits = gammalib.GFits(file+'[EVENTS][ENERGY>1.0]')
        self.test_value(fits['EVENTS'].nrows(), 152, 'Check number of rows')

        # Open FITS file with non-existing extension name
        self.test_try('Test non-existing extension name')
        try:
            fits = gammalib.GFits(file+'[DUMMY][ENERGY>1.0]')
            self.test_try_failure()
        except:
            self.test_try_success()

        # Return
        return

    # Test GFitsImage class interface
    def _test_fits_image(self):
        """
        Test GFitsImage class interface
        """
        # Set test file names
        file1 = gammalib.GFilename('test_python_fits_image_v1.fits')
        file2 = gammalib.GFilename('test_python_fits_image_v2.fits')

        # Remove test files
        file1.remove()
        file2.remove()

        # Create FITS file
        fits = gammalib.GFits(file1, True)

        # Create images
        nx   = 10
        ny   = 10
        img1 = gammalib.GFitsImageByte(nx, ny)
        img2 = gammalib.GFitsImageDouble(nx, ny)
        img3 = gammalib.GFitsImageFloat(nx, ny)
        img4 = gammalib.GFitsImageLong(nx, ny)
        img5 = gammalib.GFitsImageLongLong(nx, ny)
        img6 = gammalib.GFitsImageSByte(nx, ny)
        img7 = gammalib.GFitsImageShort(nx, ny)
        img8 = gammalib.GFitsImageULong(nx, ny)
        img9 = gammalib.GFitsImageUShort(nx, ny)
        for x in range(nx):
            for y in range(ny):
                img1[x, y] = x + y * nx
                img2[x, y] = x + y * nx
                img3[x, y] = x + y * nx
                img4[x, y] = x + y * nx
                img5[x, y] = x + y * nx
                img6[x, y] = x + y * nx
                img7[x, y] = x + y * nx
                img8[x, y] = x + y * nx
                img9[x, y] = x + y * nx
        img1.extname('Byte')
        img2.extname('Double')
        img3.extname('Float')
        img4.extname('Long')
        img5.extname('LongLong')
        img6.extname('SByte')
        img7.extname('Short')
        img8.extname('ULong')
        img9.extname('UShort')

        # Test image types by checking the bitpix values
        self.test_value(img1.bitpix(), 8)
        self.test_value(img2.bitpix(), -64)
        self.test_value(img3.bitpix(), -32)
        self.test_value(img4.bitpix(), 32)
        self.test_value(img5.bitpix(), 64)
        self.test_value(img6.bitpix(), 10)
        self.test_value(img7.bitpix(), 16)
        self.test_value(img8.bitpix(), 40)
        self.test_value(img9.bitpix(), 20)

        # Append images to FITS file
        fits.append(img1)
        fits.append(img2)
        fits.append(img3)
        fits.append(img4)
        fits.append(img5)
        # fits.append(img6) # Not supported in older cfitsio
        fits.append(img7)
        fits.append(img8)
        fits.append(img9)

        # Set header keywords
        img_byte = fits.image(0)
        img_byte.card('test', 'test-value', 'this is for testing')
        img_byte.card('real', 3.1415, 'a real value')
        img_byte.card('int', 41, 'an integer value')

        # Save FITS file
        fits.save()

        # Close FITS file
        fits.close()

        # Re-open FITS file
        fits = gammalib.GFits(file1)

        # Test image types by checking the bitpix values
        self.test_value(fits[0].bitpix(), 8)
        self.test_value(fits[1].bitpix(), -64)
        self.test_value(fits[2].bitpix(), -32)
        self.test_value(fits[3].bitpix(), 32)
        self.test_value(fits[4].bitpix(), 64)
        #self.test_value(fits[5].bitpix(), 10) # Not supported in older cfitsio
        self.test_value(fits[5].bitpix(), 16)
        self.test_value(fits[6].bitpix(), 40)
        self.test_value(fits[7].bitpix(), 20)

        # Get double precision image, take square root of pixel and save in
        # another file
        img_double = fits.image('Double')
        for x in range(nx):
            for y in range(ny):
                img_double[x, y] = math.sqrt(img_double[x, y])

        # Save into another FITS file
        fits.saveto(file2)

        # Close FITS file
        fits.close()

        # Test type conversion constructors
        img = gammalib.GFitsImageByte(img6)
        self.test_value(img.bitpix(), 8)
        img = gammalib.GFitsImageDouble(img1)
        self.test_value(img.bitpix(), -64)
        img = gammalib.GFitsImageFloat(img1)
        self.test_value(img.bitpix(), -32)
        img = gammalib.GFitsImageLong(img1)
        self.test_value(img.bitpix(), 32)
        img = gammalib.GFitsImageLongLong(img1)
        self.test_value(img.bitpix(), 64)
        img = gammalib.GFitsImageSByte(img1)
        self.test_value(img.bitpix(), 10)
        img = gammalib.GFitsImageShort(img1)
        self.test_value(img.bitpix(), 16)
        img = gammalib.GFitsImageULong(img1)
        self.test_value(img.bitpix(), 40)
        img = gammalib.GFitsImageUShort(img1)
        self.test_value(img.bitpix(), 20)

        # Return
        return

    # Test GFitsTable class interface
    def _test_fits_table(self):
        """
        Test GFitsTable class interface
        """
        # Set test file names
        file1 = gammalib.GFilename('test_python_fits_table_v1.fits')
        file2 = gammalib.GFilename('test_python_fits_table_v2.fits')

        # Remove test files
        file1.remove()
        file2.remove()

        # Create FITS file
        fits = gammalib.GFits(file1, True)

        # Create table columns
        nrows = 10
        col1 = gammalib.GFitsTableBitCol('BIT', nrows)
        col2 = gammalib.GFitsTableBoolCol('BOOLEAN', nrows)
        col3 = gammalib.GFitsTableByteCol('BYTE', nrows)
        col4 = gammalib.GFitsTableDoubleCol('DOUBLE', nrows)
        col5 = gammalib.GFitsTableFloatCol('FLOAT', nrows)
        col6 = gammalib.GFitsTableLongCol('LONG', nrows)
        col7 = gammalib.GFitsTableLongLongCol('LONGLONG', nrows)
        col8 = gammalib.GFitsTableShortCol('SHORT', nrows)
        col9 = gammalib.GFitsTableStringCol('STRING', nrows, 20)
        col10 = gammalib.GFitsTableULongCol('ULONG', nrows)
        col11 = gammalib.GFitsTableUShortCol('USHORT', nrows)
        for i in range(nrows):
            #col1[i] = i % 2        # Old swig version
            #col2[i] = i % 2        # Old swig version
            col1[i] = bool(i % 2)  # New swig version (3.x.y)
            col2[i] = bool(i % 2)  # New swig version (3.x.y)
            col3[i] = i
            col4[i] = i * 0.01
            col5[i] = i * 0.01
            col6[i] = i * 100
            col7[i] = i * 10000
            col8[i] = i * 100
            col9[i] = str(i * 100)
            col10[i] = i * 100
            col11[i] = i * 100

        # Set ASCII table
        tbl_ascii = gammalib.GFitsAsciiTable()
        # tbl_ascii.append(col1) # Need to implement ?/!
        # tbl_ascii.append(col2) # Need to implement ?/!
        tbl_ascii.append(col3)
        tbl_ascii.append(col4)
        tbl_ascii.append(col5)
        tbl_ascii.append(col6)
        tbl_ascii.append(col7)
        tbl_ascii.append(col8)
        tbl_ascii.append(col9)
        tbl_ascii.append(col10)
        tbl_ascii.append(col11)
        tbl_ascii.extname('ASCII table')
        fits.append(tbl_ascii)

        # Set binary table
        tbl_bin = gammalib.GFitsBinTable()
        tbl_bin.append(col1)
        tbl_bin.append(col2)
        tbl_bin.append(col3)
        tbl_bin.append(col4)
        tbl_bin.append(col5)
        tbl_bin.append(col6)
        tbl_bin.append(col7)
        tbl_bin.append(col8)
        tbl_bin.append(col9)
        tbl_bin.append(col10)
        tbl_bin.append(col11)
        tbl_bin.extname('Binary table')
        fits.append(tbl_bin)

        # Save FITS file
        fits.save()

        # Close FITS file
        fits.close()

        # Return
        return

    # Test FITS table columns
    def _test_fits_table_columns(self):
        """
        Test FITS table columns
        """
        # Create table columns
        nrows = 10
        col1  = gammalib.GFitsTableBitCol('BIT', nrows)
        col2  = gammalib.GFitsTableBoolCol('BOOLEAN', nrows)
        col3  = gammalib.GFitsTableByteCol('BYTE', nrows)
        col4  = gammalib.GFitsTableDoubleCol('DOUBLE', nrows)
        col5  = gammalib.GFitsTableFloatCol('FLOAT', nrows)
        col6  = gammalib.GFitsTableLongCol('LONG', nrows)
        col7  = gammalib.GFitsTableLongLongCol('LONGLONG', nrows)
        col8  = gammalib.GFitsTableShortCol('SHORT', nrows)
        col9  = gammalib.GFitsTableStringCol('STRING', nrows, 20)
        col10 = gammalib.GFitsTableULongCol('ULONG', nrows)
        col11 = gammalib.GFitsTableUShortCol('USHORT', nrows)

        # Test number of rows
        self.test_value(col1.nrows(),  nrows, 'Check number of rows in Bit column')
        self.test_value(col2.nrows(),  nrows, 'Check number of rows in Boolean column')
        self.test_value(col3.nrows(),  nrows, 'Check number of rows in Byte column')
        self.test_value(col4.nrows(),  nrows, 'Check number of rows in Double column')
        self.test_value(col5.nrows(),  nrows, 'Check number of rows in Float column')
        self.test_value(col6.nrows(),  nrows, 'Check number of rows in Long column')
        self.test_value(col7.nrows(),  nrows, 'Check number of rows in LongLong column')
        self.test_value(col8.nrows(),  nrows, 'Check number of rows in Short column')
        self.test_value(col9.nrows(),  nrows, 'Check number of rows in String column')
        self.test_value(col10.nrows(), nrows, 'Check number of rows in ULong column')
        self.test_value(col11.nrows(), nrows, 'Check number of rows in UShort column')

        # Test iterators
        for row in col1:
            pass
        for row in col2:
            pass
        for row in col3:
            pass
        for row in col4:
            pass
        for row in col5:
            pass
        for row in col6:
            pass
        for row in col7:
            pass
        for row in col8:
            pass
        for row in col9:
            pass
        for row in col10:
            pass
        for row in col11:
            pass

        # Test setting and retrieving of values
        col1[5] = True
        self.test_assert(col1[5], 'Check Bit column setting and retrieving')
        col2[5] = True
        self.test_assert(col2[5], 'Check Boolean column setting and retrieving')
        col3[5] = 5
        self.test_value(col3[5], 5, 'Check Byte column setting and retrieving')
        col4[5] = 3.14
        self.test_value(col4[5], 3.14, 1.0e-6, 'Check Double column setting and retrieving')
        col5[5] = 3.14
        self.test_value(col5[5], 3.14, 1.0e-6, 'Check Float column setting and retrieving')
        col6[5] = 314
        self.test_value(col6[5], 314, 'Check Long column setting and retrieving')
        col7[5] = 314
        self.test_value(col7[5], 314, 'Check LongLong column setting and retrieving')
        col8[5] = 314
        self.test_value(col8[5], 314, 'Check Short column setting and retrieving')
        col9[5] = 'Hallo'
        self.test_assert(col9[5] == 'Hallo', 'Check String column setting and retrieving')
        col10[5] = 314
        self.test_value(col10[5], 314, 'Check ULong column setting and retrieving')
        col11[5] = 314
        self.test_value(col11[5], 314, 'Check UShort column setting and retrieving')

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GFits())
        test_support.pickeling(self, gammalib.GFitsAsciiTable())
        test_support.pickeling(self, gammalib.GFitsBinTable())
        test_support.pickeling(self, gammalib.GFitsHeader())
        test_support.pickeling(self, gammalib.GFitsHeaderCard())
        test_support.pickeling(self, gammalib.GFitsImageByte())
        test_support.pickeling(self, gammalib.GFitsImageDouble())
        test_support.pickeling(self, gammalib.GFitsImageFloat())
        test_support.pickeling(self, gammalib.GFitsImageLong())
        test_support.pickeling(self, gammalib.GFitsImageLongLong())
        test_support.pickeling(self, gammalib.GFitsImageSByte())
        test_support.pickeling(self, gammalib.GFitsImageShort())
        test_support.pickeling(self, gammalib.GFitsImageULong())
        test_support.pickeling(self, gammalib.GFitsImageUShort())
        test_support.pickeling(self, gammalib.GFitsTableBitCol())
        test_support.pickeling(self, gammalib.GFitsTableBoolCol())
        test_support.pickeling(self, gammalib.GFitsTableByteCol())
        #test_support.pickeling(self, gammalib.GFitsTableCDoubleCol())
        #test_support.pickeling(self, gammalib.GFitsTableCFloatCol())
        test_support.pickeling(self, gammalib.GFitsTableDoubleCol())
        test_support.pickeling(self, gammalib.GFitsTableFloatCol())
        test_support.pickeling(self, gammalib.GFitsTableLongCol())
        test_support.pickeling(self, gammalib.GFitsTableLongLongCol())
        test_support.pickeling(self, gammalib.GFitsTableShortCol())
        test_support.pickeling(self, gammalib.GFitsTableStringCol())
        test_support.pickeling(self, gammalib.GFitsTableULongCol())
        test_support.pickeling(self, gammalib.GFitsTableUShortCol())

        # Setup for tests
        fits   = gammalib.GFits(os.environ['TEST_DATA']+'/test_cube.fits')
        header = gammalib.GFitsHeader()
        header.append(gammalib.GFitsHeaderCard('key1','"string"','deg','test1'))
        header.append(gammalib.GFitsHeaderCard('key2','1.0','deg','test2'))
        col_bit    = gammalib.GFitsTableBitCol('a',2,3)
        col_bool   = gammalib.GFitsTableBoolCol('a',2,3)
        col_byte   = gammalib.GFitsTableByteCol('a',2,3)
        col_double = gammalib.GFitsTableDoubleCol('a',2,3)
        col_float  = gammalib.GFitsTableFloatCol('a',2,3)
        col_long   = gammalib.GFitsTableLongCol('a',2,3)
        col_llong  = gammalib.GFitsTableLongLongCol('a',2,3)
        col_short  = gammalib.GFitsTableShortCol('a',2,3)
        col_ulong  = gammalib.GFitsTableULongCol('a',2,3)
        col_ushort = gammalib.GFitsTableUShortCol('a',2,3)
        col_string = gammalib.GFitsTableStringCol('a',2,10,3)
        bit        = False
        for row in range(2):
            for col in range(3):
                if bit:
                    bit = False
                else:
                    bit = True
                col_bit[row,col]    = bit
                col_bool[row,col]   = not bit
                col_byte[row,col]   = row*10+col
                col_double[row,col] = row*10+col
                col_float[row,col]  = row*10+col
                col_long[row,col]   = row*10+col
                col_llong[row,col]  = row*10+col
                col_short[row,col]  = row*10+col
                col_ulong[row,col]  = row*10+col
                col_ushort[row,col] = row*10+col
                col_string[row,col] = str(row*10+col)

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GFits(fits))
        test_support.pickeling(self, gammalib.GFitsAsciiTable(3))
        test_support.pickeling(self, gammalib.GFitsBinTable(4))
        test_support.pickeling(self, gammalib.GFitsHeader(header))
        test_support.pickeling(self, gammalib.GFitsHeaderCard('key','value','deg','test'))
        test_support.pickeling(self, gammalib.GFitsImageByte(2,2,[1,2,3,4]))
        test_support.pickeling(self, gammalib.GFitsImageDouble(2,2,[1.0,2.0,3.0,4.0]))
        test_support.pickeling(self, gammalib.GFitsImageFloat(2,2,[1.0,2.0,3.0,4.0]))
        test_support.pickeling(self, gammalib.GFitsImageLong(2,2,[1,2,3,4]))
        test_support.pickeling(self, gammalib.GFitsImageLongLong(2,2,[1,2,3,4]))
        test_support.pickeling(self, gammalib.GFitsImageSByte(2,2,[1,2,3,4]))
        test_support.pickeling(self, gammalib.GFitsImageShort([2,2],[1,2,3,4]))
        test_support.pickeling(self, gammalib.GFitsImageULong(2,2,1,[1,2,3,4]))
        test_support.pickeling(self, gammalib.GFitsImageUShort(3,[1,2,3]))
        test_support.pickeling(self, gammalib.GFitsTableBitCol(col_bit))
        test_support.pickeling(self, gammalib.GFitsTableBoolCol(col_bool))
        test_support.pickeling(self, gammalib.GFitsTableByteCol(col_byte))
        #test_support.pickeling(self, gammalib.GFitsTableCDoubleCol())
        #test_support.pickeling(self, gammalib.GFitsTableCFloatCol())
        test_support.pickeling(self, gammalib.GFitsTableDoubleCol(col_double))
        test_support.pickeling(self, gammalib.GFitsTableFloatCol(col_float))
        test_support.pickeling(self, gammalib.GFitsTableLongCol(col_long))
        test_support.pickeling(self, gammalib.GFitsTableLongLongCol(col_llong))
        test_support.pickeling(self, gammalib.GFitsTableShortCol(col_short))
        test_support.pickeling(self, gammalib.GFitsTableStringCol(col_string))
        test_support.pickeling(self, gammalib.GFitsTableULongCol(col_ulong))
        test_support.pickeling(self, gammalib.GFitsTableUShortCol(col_ushort))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('fits')

        # Append tests
        self.append(self._test_fits, 'Test GFits')
        self.append(self._test_fits_access, 'Test GFits member access')
        self.append(self._test_fits_slicing, 'Test GFits slicing')
        self.append(self._test_header_access, 'Test GFitsHeader member access')
        self.append(self._test_header_slicing, 'Test GFitsHeader slicing')
        self.append(self._test_fits_image, 'Test GFitsImage')
        self.append(self._test_fits_table, 'Test GFitsTable')
        self.append(self._test_fits_table_columns, 'Test GFitsTableCol')
        self.append(self._test_pickeling, 'Test Fits class pickeling')

        # Return
        return
