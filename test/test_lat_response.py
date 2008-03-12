#! /usr/bin/env python

from gammalib import *
from pylab import *

#==============#
# Test LAT PSF #
#==============#
def test_psf():
    """
    Test GammaLib LAT PSF interface.
    """
    # Allocate LAT response
    front = GLATResponse()
    back  = GLATResponse()
    
    # Set calibration database
    front.set_caldb("irf/lat")
    back.set_caldb("irf/lat")
    
    # Load response
    front.load("Pass5_v0", "front")
    back.load("Pass5_v0",  "back")
    
    # Save response
    front.save("test_rsp_pass5_front.fits")
    back.save("test_rsp_pass5_back.fits")

    # Build offset array
    nbins = 100000
    bin   = 0.0174532 * 0.001
    delta = GVector(nbins)
    for i in range(nbins):
        delta[i] = i*bin

    # Set cos theta
    ctheta = 0.70

    # Figure 1: Front PSF using vector
    figure(1)
    xlabel('Angular distance (deg)')
    ylabel('Value')

    # Get PSF
    x_int = []
    y_int = []
    print "Front - vector access"
    for k in range(17):
        logE   = 1.2+k*0.3  # [1.2 - 5.7]
        psf    = front.psf(delta, logE, ctheta)
        sum    = psf * delta.sin()
        sum   *= 2.0 * 3.1415926535897931159979635 * bin
        print "E=" + str(int(pow(10.0,logE)))+ " MeV, Sum=" + str(sum)
        x_int.append(pow(10.0,logE))
        y_int.append(sum)
    
        # Build plot
        x_val = []
        y_val = []
        for i in range(nbins):
            x_val.append(delta[i]/0.0174532)
            y_val.append(psf[i])
        loglog(x_val, y_val)
    
    # Plot integral versus energy
    figure(2)
    semilogx(x_int, y_int)
    xlabel('Energy (MeV)')
    ylabel('integral')

    # Figure 3: Back PSF using vector
    figure(3)
    xlabel('Angular distance (deg)')
    ylabel('Value')

    # Get PSF
    x_int = []
    y_int = []
    print "Back - vector access"
    for k in range(17):
        logE   = 1.2+k*0.3  # [1.2 - 5.7]
        psf    = back.psf(delta, logE, ctheta)
        sum    = psf * delta.sin()
        sum   *= 2.0 * 3.1415926535897931159979635 * bin
        print "E=" + str(int(pow(10.0,logE)))+ " MeV, Sum=" + str(sum)
        x_int.append(pow(10.0,logE))
        y_int.append(sum)
    
        # Build plot
        x_val = []
        y_val = []
        for i in range(nbins):
            x_val.append(delta[i]/0.0174532)
            y_val.append(psf[i])
        loglog(x_val, y_val)
    
    # Plot integral versus energy
    figure(4)
    semilogx(x_int, y_int)
    xlabel('Energy (MeV)')
    ylabel('integral')

    # Figure 5: Front PSF using vector
    figure(5)
    xlabel('Angular distance (deg)')
    ylabel('Value')

    # Get PSF
    x_int = []
    y_int = []
    print "Front - direct access"
    for k in range(17):
        logE   = 1.2+k*0.3  # [1.2 - 5.7]
        for i in range(nbins):
            psf[i] = front.psf(delta[i], logE, ctheta)
        sum    = psf * delta.sin()
        sum   *= 2.0 * 3.1415926535897931159979635 * bin
        print "E=" + str(int(pow(10.0,logE)))+ " MeV, Sum=" + str(sum)
        x_int.append(pow(10.0,logE))
        y_int.append(sum)
    
        # Build plot
        x_val = []
        y_val = []
        for i in range(nbins):
            x_val.append(delta[i]/0.0174532)
            y_val.append(psf[i])
        loglog(x_val, y_val)
    
    # Plot integral versus energy
    figure(6)
    semilogx(x_int, y_int)
    xlabel('Energy (MeV)')
    ylabel('integral')


    # Show plots
    show()


#===============#
# Test LAT Aeff #
#===============#
def test_aeff():
    """
    Test GammaLib LAT Aeff interface.
    """
    # Allocate LAT response
    front = GLATResponse()
    back  = GLATResponse()
    
    # Set calibration database
    front.set_caldb("irf/lat")
    back.set_caldb("irf/lat")
    
    # Load response
    front.load("Pass5_v0", "front")
    back.load("Pass5_v0",  "back")
    
    # Save response
    front.save("test_rsp_pass5_front.fits")
    back.save("test_rsp_pass5_back.fits")
    
    # Restrict to 70 deg zenith angle
    front.aeff_ctheta_min(cos(70.0*pi/180))
    back.aeff_ctheta_min(cos(70.0*pi/180))

    # Figure 1: front
    figure(1)
    make_aeff_panel(front)

    # Figure 2: back
    figure(2)
    make_aeff_panel(back)

    # Show figures
    show()

def make_aeff_panel(rsp):
    """
    Make 6 figures.
    """
    # logE=1.25527 (18 MeV)
    subplot(231)
    make_aeff_image(1.25527, rsp)
    title("Front, 18 MeV")

    # logE=1.47712
    subplot(232)
    make_aeff_image(1.47712, rsp)
    title("Front, 30 MeV")

    # logE=2.0
    subplot(233)
    make_aeff_image(2.0, rsp)
    title("Front, 100 MeV")

    # logE=3.0
    subplot(234)
    make_aeff_image(3.0, rsp)
    title("Front, 1 GeV")

    # logE=4.0
    subplot(235)
    make_aeff_image(4.0, rsp)
    title("Front, 10 GeV")

    # logE=5.0
    subplot(236)
    make_aeff_image(5.0, rsp)
    title("Front, 100 GeV")


def make_aeff_image(logE, rsp):
    """
    Make Aeff image.
    """
    # Setup results
    aeff_list = []
    x  = arange(-90.0, 90.0, 1.0)
    y  = arange(-90.0, 90.0, 1.0)
    nx = len(x)
    ny = len(y)
    for xval in x:
        for yval in y:
            angle  = sqrt(xval*xval+yval*yval)
            ctheta = cos(angle*pi/180.0)
            value  = rsp.aeff(logE, ctheta)
            aeff_list.append(value)
    
    # Build image array
    aeff = array(aeff_list)
    aeff.shape = nx, ny
    im = imshow(aeff)
    #colorbar()
    

#=================#
# Test Node array #
#=================#
def test_node_array():
    """
    Test GammaLib GNodeArray interface.
    """
    # Set-up linear vector and data array
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
        if (wgt_left+wgt_right != 1.0):
            print "WARNING (linear): Weights do not sum to 1.0 ("+str(wgt_left+wgt_right)+")"
        x_val.append(x)
        y_val.append(y)
    
    # Plot Figure 1
    figure(1)
    plot(x_val, y_val)
    
    # Set-up non-linear vector and data array
    vector = GVector(20)
    data   = GVector(20)
    for i in range(20):
        vector[i] = 10.0 + i*i*1.0
        data[i]   = sin(0.7*(i))
    
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
        if (wgt_left+wgt_right != 1.0):
            print "WARNING (bisection): Weights do not sum to 1.0 ("+str(wgt_left+wgt_right)+")"
        x_val.append(x)
        y_val.append(y)
    
    # Plot Figure 2
    figure(2)
    plot(x_val, y_val)
    
    # Show figures
    show()
    

#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    #test_psf()
    test_aeff()
    #test_node_array()