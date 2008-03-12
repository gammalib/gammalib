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
    

#=================#
# Test Node array #
#=================#
def test_node_array():
    """
    Test GammaLib GNodeArray interface.
    """
    # Set-up vector and data array
    vector = GVector(20)
    data   = GVector(20)
    for i in range(20):
        vector[i] = 10.0 + i*5.0
        data[i]   = sin(0.15*(vector[i]-10.0))
    #print vector
    
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
    
    plot(x_val, y_val)
    show()
    

#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    test_psf()
    #test_node_array()