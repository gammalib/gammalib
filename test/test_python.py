#! /usr/bin/env python

from gammalib import *
from math import *
import os

#================#
# Test FITS file #
#================#
def test_fits():
    """
    Test GFits interface.
    """
    # Dump testing
    print "Test GFits."

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
    fits = GFits(file1)

    # Create images
    nx   = 10
    ny   = 10
    img1 = GFitsImageByte(nx,ny)
    img2 = GFitsImageDouble(nx,ny)
    img3 = GFitsImageFloat(nx,ny)
    img4 = GFitsImageLong(nx,ny)
    img5 = GFitsImageLongLong(nx,ny)
    img6 = GFitsImageSByte(nx,ny)
    img7 = GFitsImageShort(nx,ny)
    img8 = GFitsImageULong(nx,ny)
    img9 = GFitsImageUShort(nx,ny)
    for x in range(nx):
        for y in range(ny):
            img1[x,y] = x+y*nx
            img2[x,y] = x+y*nx
            img3[x,y] = x+y*nx
            img4[x,y] = x+y*nx
            img5[x,y] = x+y*nx
            img6[x,y] = x+y*nx
            img7[x,y] = x+y*nx
            img8[x,y] = x+y*nx
            img9[x,y] = x+y*nx
    img1.extname("Byte")
    img2.extname("Double")
    img3.extname("Float")
    img4.extname("Long")
    img5.extname("LongLong")
    img6.extname("SByte")
    img7.extname("Short")
    img8.extname("ULong")
    img9.extname("UShort")
    
    # Append images to FITS file
    fits.append(img1)
    fits.append(img2)
    fits.append(img3)
    fits.append(img4)
    fits.append(img5)
    fits.append(img6)
    fits.append(img7)
    fits.append(img8)
    fits.append(img9)

    # Set header keywords
    img_byte = fits.image(0)
    img_byte.card("test", "test-value", "this is for testing")
    img_byte.card("real", 3.1415, "a real value")
    img_byte.card("int", 41, "an integer value")

    # Create table columns
    nrows = 10
    col1  = GFitsTableBitCol("BIT", nrows)
    col2  = GFitsTableBoolCol("BOOLEAN", nrows)
    col3  = GFitsTableByteCol("BYTE", nrows)
    col4  = GFitsTableDoubleCol("DOUBLE", nrows)
    col5  = GFitsTableFloatCol("FLOAT", nrows)
    col6  = GFitsTableLongCol("LONG", nrows)
    col7  = GFitsTableLongLongCol("LONGLONG", nrows)
    col8  = GFitsTableShortCol("SHORT", nrows)
    col9  = GFitsTableStringCol("STRING", nrows, 20)
    col10 = GFitsTableULongCol("ULONG", nrows)
    col11 = GFitsTableUShortCol("USHORT", nrows)
    for i in range(nrows):
        col1[i] = i % 2
        col2[i] = i % 2
        col3[i] = i
        col4[i] = i*0.01
        col5[i] = i*0.01
        col6[i] = i*100
        col7[i] = i*10000
        col8[i] = i*100
        col9[i] = str(i*100)
        col10[i] = i*100
        col11[i] = i*100

    # Set ASCII table
    tbl_ascii = GFitsAsciiTable()
    #tbl_ascii.append_column(col1) # Need to implement ?/!
    #tbl_ascii.append_column(col2) # Need to implement ?/!
    tbl_ascii.append_column(col3)
    tbl_ascii.append_column(col4)
    tbl_ascii.append_column(col5)
    tbl_ascii.append_column(col6)
    tbl_ascii.append_column(col7)
    tbl_ascii.append_column(col8)
    tbl_ascii.append_column(col9)
    tbl_ascii.append_column(col10)
    tbl_ascii.append_column(col11)
    tbl_ascii.extname("ASCII table")
    fits.append(tbl_ascii)

    # Set binary table
    tbl_bin = GFitsBinTable()
    tbl_bin.append_column(col1)
    tbl_bin.append_column(col2)
    tbl_bin.append_column(col3)
    tbl_bin.append_column(col4)
    tbl_bin.append_column(col5)
    tbl_bin.append_column(col6)
    tbl_bin.append_column(col7)
    tbl_bin.append_column(col8)
    tbl_bin.append_column(col9)
    tbl_bin.append_column(col10)
    tbl_bin.append_column(col11)
    tbl_bin.extname("Binary table")
    fits.append(tbl_bin)

    # Save FITS file
    fits.save()

    # Close FITS file
    fits.close()

    # Re-open FITS file
    fits = GFits(file1)

    # Get double precision image, take square root of pixel and save in
    # another file
    img_double = cast_double(fits.image("Double"))
    for x in range(nx):
        for y in range(ny):
            img_double[x,y] = sqrt(img_double[x,y])
    #img_byte = cast_byte(fits.image("Double"))

    # Save into another FITS file
    fits.saveto(file2)

    # Close FITS file
    fits.close()


#==============#
# Test GSkymap #
#==============#
def test_skymap():
    """
    Test GSkymap interface.
    """
    # Dump testing
    print "Test GSkymap."

    # Set filenames
    file1 = "test_python_skymap_hpx_v1.fits"
    file2 = "test_python_skymap_hpx_v2.fits"

    # Remove test files
    try:
        os.remove(file1)
        os.remove(file2)
    except:
        pass
    
    # Create HEALPix skymap
    pixels = GSkymap("HPX", "GAL", 4, "RING", 2)
    for i, pixel in enumerate(pixels):
        pixels[i]   = i+1.0
        pixels[i:1] = i+1.0 + 1000.0
    pixels.save(file1)
    
    # Load HEALPix skymap
    pixels = GSkymap(file1)

    # Dump skymap
    #print pixels

    # Control pix2dir / dir2pix
    for i in range(pixels.npix()):
        dir  = pixels.pix2dir(i)
        ipix = pixels.dir2pix(dir)
        if (i != ipix):
            print "GSkymap Trouble with pixel "+str(i)+" ("+str(ipix)+ \
                  "), RA="+str(dir.ra()*180/pi)+", Dec="+str(dir.dec()*180/pi)

    # Control SkyDir coordinate transformation for all pixels
    for i in range(pixels.npix()):
        dir     = pixels.pix2dir(i)
        dir_new = GSkyDir()
        dir_new.lb(dir.l(),dir.b())
        dra  = abs(dir.ra()  - dir_new.ra())
        if (dra >= 5.0):
            dra -= 2.0*pi
        ddec = abs(dir.dec() - dir_new.dec())
        if (dra > 1.0e-9 or ddec > 1.0e-9):
            print "GSkyDir Trouble with pixel "+str(i)+" ("+str(dra)+ \
                  ","+str(ddec)+")"
    
    # Save HEALPix skymap twice. The second saving should fail.
    try:
        pixels.save(file2, True)
        pixels.save(file2)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("*** TEST ERROR: FITS file overwritten!")
    pixels.save(file2, True)

    # Load again HEALPix skymap
    pixels = GSkymap(file1)

    # Dump skymap
    #print pixels


#=================#
# Test Node array #
#=================#
def test_node_array():
    """
    Test GNodeArray interface.
    """
    # Set-up vector and data array
    vector = GVector(20)
    data   = GVector(20)
    for i in range(20):
        vector[i] = 10.0 + i*5.0
        data[i]   = sin(0.15*(vector[i]-10.0))

    # Set-up node array
    array = GNodeArray()
    array.nodes(vector)

    # Get values
    x_val = []
    y_val = []
    for i in range(100):
        x = i-10
        array.set_value(x)
        inx_left  = array.inx_left()
        inx_right = array.inx_right()
        wgt_left  = array.wgt_left()
        wgt_right = array.wgt_right()
        y         = wgt_left*data[inx_left] + wgt_right*data[inx_right]
        x_val.append(x)
        y_val.append(y)


#===================#
# Test LAT response #
#===================#
def test_lat_response():
    """
    Test GLATResponse interface.
    """
    # Remove test file
    try:
        os.remove("test_rsp.fits")
    except:
        pass

    # Allocate LAT response
    rsp = GLATResponse()

    # Set calibration database
    rsp.set_caldb("irf/lat")

    # Load response
    rsp.load("Pass5_v0", "front")

    # Save response
    rsp.save("test_rsp.fits")


#======================#
# Test LAT observation #
#======================#
def test_lat_observation():
    """
    Test GLATObservation interface.
    """
    # Allocate LAT Observation
    obs = GLATObservation("data/FT1_253582800.fits.gz", "data/FT2_253582800.fits.gz")

    #
    #print obs.ft1()
    #print obs.ft2()


#===================#
# Test optimisation #
#===================#
def test_optimise():
    """
    """
    # Load CTA run
    run_cta = GCTAObservation()
    run_cta.load_binned("../inst/cta/test/data/run_00006028_cntmap.fits.gz")
    run_cta.response("kb_E_50h_v3", "../inst/cta/test/irf")

    # Build observation list
    obs = GObservations()
    obs.append(run_cta)

    # Setup model for optimizing
    dir          = GSkyDir()
    dir.radec_deg(117.0, -33.0)
    point_source = GModelSpatialPtsrc(dir)
    power_law    = GModelSpectralPlaw(1.0e-7, -2.1)
    power_law.par(0).min(1.0e-12)
    crab         = GModel(point_source, power_law)
    crab.name("Crab")
    models       = GModels()
    models.append(crab)
    obs.models(models)

    # Perform LM optimization
    opt = GOptimizerLM()
    opt.max_iter(1000)
    obs.optimize(opt)
    print opt
    print obs


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    # Dump result
    print
    print "****************************"
    print "* Python interface testing *"
    print "****************************"

    # Initialise success counter
    tests   = 0
    success = 0
    
    # Perform tests
    test_fits()
    test_skymap()
#    test_node_array()
#    test_lat_response()
#    test_lat_observation()
#    test_optimise()
