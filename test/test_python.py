#! /usr/bin/env python

from gammalib import *


#================#
# Test FITS file #
#================#
def test_fits():
    """
    Test GammaLib GFits interface.
    """
    # Create FITS file
    fits = GFits()
    fits.open("test_python.fits")

    # Print number of HDUs
    print "Number of HDUs ..: " + str(fits.num_hdus())

    # Allocate HDU "Testing"
    hdu = GFitsHDU()
    hdu.extname("Testing")
    print "HDU name ........: " + hdu.extname()
    print "HDU number ......: " + str(hdu.extno())
    print "HDU type ........: " + str(hdu.exttype())

    # Set header keywords
    hdr = hdu.header()
    hdr.update(GFitsHeaderCard("test", "test-value", "this is for testing"))
    hdr.update(GFitsHeaderCard("real", 3.1415, "a real value"))
    hdr.update(GFitsHeaderCard("int", 41, "an integer value"))
    

    # Append HDU to FITS file
    fits.append_hdu(hdu)

    # Set double image
    img1 = GFitsDblImage(100)
    img2 = GFitsDblImage(10, 20)
    img3 = GFitsDblImage(10, 10, 10)
    img4 = GFitsDblImage(5, 9, 20, 20)
    print img1
    print img2
    print img3
    print img4
    print img1.bitpix()
    print img3.naxis()
    print img2.naxes(0),  img2.naxes(1)
    print img4.num_pixels()
    fits.append_hdu(GFitsHDU(img1))
    fits.append_hdu(GFitsHDU(img2))
    fits.append_hdu(GFitsHDU(img3))
    fits.append_hdu(GFitsHDU(img4))
    
    # Set ASCII table
    tbl1 = GFitsAsciiTable(10)
    col1 = GFitsTableDblCol("DOUBLE", 10)
    col2 = GFitsTableFltCol("FLOAT", 10)
    col3 = GFitsTableLngCol("LONG", 10)
    col4 = GFitsTableShtCol("SHORT", 10)
    col5 = GFitsTableStrCol("STRING", 10, 20)
    for i in range(10):
        col1.set(i, i*0.01)
        col2.set(i, i*0.01)
        col3.set(i, i*100)
        col4.set(i, i*100)
        col5.set(i, str(i*100))
    tbl1.append_column(col1)
    #tbl1.append_column(col2) # FLOAT: datatype conversion overflow (status=412)
    #tbl1.append_column(col3) # LONG: unknown TFORM datatype code (status=262)
    #tbl1.append_column(col4) # SHORT: datatype conversion overflow (status=412)
    tbl1.append_column(col5)
    fits.append_hdu(GFitsHDU(tbl1))

    # Set binary table
    tbl2 = GFitsBinTable(10)
    col1 = GFitsTableDblCol("DOUBLE", 10)
    col2 = GFitsTableFltCol("FLOAT", 10)
    col3 = GFitsTableLngCol("LONG", 10)
    col4 = GFitsTableShtCol("SHORT", 10)
    col5 = GFitsTableStrCol("STRING", 10, 20)
    for i in range(10):
        col1.set(i, i*0.01)
        col2.set(i, i*0.01)
        col3.set(i, i*100)
        col4.set(i, i*100)
        col5.set(i, str(i*100))
    tbl2.append_column(col1)
    tbl2.append_column(col2)
    tbl2.append_column(col3)
    tbl2.append_column(col4)
    tbl2.append_column(col5)
    fits.append_hdu(GFitsHDU(tbl2))

    #print fits
    
    # Save FITS file
    fits.save()

    # Save into another FITS file
    fits.saveto("test_python_2.fits")


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    test_fits()
