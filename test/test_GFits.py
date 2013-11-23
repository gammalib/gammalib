# ==========================================================================
# This module performs unit tests for the GammaLib FITS module.
#
# Copyright (C) 2012 Juergen Knoedlseder
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
from gammalib import *
from math import *
import os
import sys


# =================================== #
# Test class for GammaLib FITS module #
# =================================== #
class Test(GPythonTestSuite):
    """
    Test class for GammaLib FITS module.
    """
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        GPythonTestSuite.__init__(self)

        # Return
        return

    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("fits")

        # Append tests
        self.append(self.test_fits, "Test GFits")

        # Return
        return

    def test_fits(self):
        """
        Test FITS interface.
        """
        # Set filenames
        file1 = "test_python_fits_v1.fits"
        file2 = "test_python_fits_v2.fits"

        # Remove test files
        try:
            os.remove(file1)
            os.remove(file2)
        except:
            pass

        # Create FITS file
        fits = GFits(file1, True)
        sys.stdout.write(".")

        # Create images
        nx = 10
        ny = 10
        img1 = GFitsImageByte(nx, ny)
        img2 = GFitsImageDouble(nx, ny)
        img3 = GFitsImageFloat(nx, ny)
        img4 = GFitsImageLong(nx, ny)
        img5 = GFitsImageLongLong(nx, ny)
        img6 = GFitsImageSByte(nx, ny)
        img7 = GFitsImageShort(nx, ny)
        img8 = GFitsImageULong(nx, ny)
        img9 = GFitsImageUShort(nx, ny)
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
        img1.extname("Byte")
        img2.extname("Double")
        img3.extname("Float")
        img4.extname("Long")
        img5.extname("LongLong")
        img6.extname("SByte")
        img7.extname("Short")
        img8.extname("ULong")
        img9.extname("UShort")
        sys.stdout.write(".")

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
        sys.stdout.write(".")

        # Set header keywords
        img_byte = fits.image(0)
        img_byte.card("test", "test-value", "this is for testing")
        img_byte.card("real", 3.1415, "a real value")
        img_byte.card("int", 41, "an integer value")
        sys.stdout.write(".")

        # Create table columns
        nrows = 10
        col1 = GFitsTableBitCol("BIT", nrows)
        col2 = GFitsTableBoolCol("BOOLEAN", nrows)
        col3 = GFitsTableByteCol("BYTE", nrows)
        col4 = GFitsTableDoubleCol("DOUBLE", nrows)
        col5 = GFitsTableFloatCol("FLOAT", nrows)
        col6 = GFitsTableLongCol("LONG", nrows)
        col7 = GFitsTableLongLongCol("LONGLONG", nrows)
        col8 = GFitsTableShortCol("SHORT", nrows)
        col9 = GFitsTableStringCol("STRING", nrows, 20)
        col10 = GFitsTableULongCol("ULONG", nrows)
        col11 = GFitsTableUShortCol("USHORT", nrows)
        for i in range(nrows):
            col1[i] = i % 2
            col2[i] = i % 2
            col3[i] = i
            col4[i] = i * 0.01
            col5[i] = i * 0.01
            col6[i] = i * 100
            col7[i] = i * 10000
            col8[i] = i * 100
            col9[i] = str(i * 100)
            col10[i] = i * 100
            col11[i] = i * 100
        sys.stdout.write(".")

        # Set ASCII table
        tbl_ascii = GFitsAsciiTable()
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
        tbl_ascii.extname("ASCII table")
        fits.append(tbl_ascii)
        sys.stdout.write(".")

        # Set binary table
        tbl_bin = GFitsBinTable()
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
        tbl_bin.extname("Binary table")
        fits.append(tbl_bin)
        sys.stdout.write(".")

        # Save FITS file
        # sys.stdout.write(fits+"\n")
        fits.save()
        sys.stdout.write(".")

        # Close FITS file
        fits.close()
        sys.stdout.write(".")

        # Re-open FITS file
        fits = GFits(file1)
        sys.stdout.write(".")

        # Get double precision image, take square root of pixel and save in
        # another file
        img_double = fits.image("Double")
        for x in range(nx):
            for y in range(ny):
                img_double[x, y] = sqrt(img_double[x, y])
        sys.stdout.write(".")

        # Save into another FITS file
        fits.saveto(file2)
        sys.stdout.write(".")

        # Close FITS file
        fits.close()
        sys.stdout.write(".")

        # Return
        return
