# ==========================================================================
# This module performs unit tests for the GammaLib numerics module.
#
# Copyright (C) 2012-2016 Juergen Knoedlseder
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
import gammalib


# ======================================= #
# Test class for GammaLib numerics module #
# ======================================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib numerics module
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

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('numerics')

        # Append tests
        self.append(self.test_fft, 'Test GFft')

        # Return
        return

    # Test Fast Fourier Transform
    def test_fft(self):
        """
        Test function
        """
        # Allocate and set 2-d array
        array = gammalib.GNdarray(10, 10)
        for i in range(3,7):
            array[i,i] = 1.0
        ref = array.sum()
        
        # Allocate and set 2-d kernel. The kernel needs to be normalized
        # to unity and the centre of the kernel needs to be at pixel [0,0],
        # and it needs to be wraped around to negative indices
        kernel = gammalib.GNdarray(10, 10)
        kernel[0,0] = 0.4
        # 0.4
        kernel[0,1] = 0.1
        kernel[1,0] = 0.1
        kernel[0,9] = 0.1
        kernel[9,0] = 0.1
        # 0.2
        kernel[1,1] = 0.05
        kernel[9,1] = 0.05
        kernel[9,9] = 0.05
        kernel[1,9] = 0.05

        # Smooth 2-d array using FFT
        fft_array  = gammalib.GFft(array)
        fft_kernel = gammalib.GFft(kernel)
        fft_smooth = fft_array * fft_kernel

        # Backtransform
        smooth = fft_smooth.backward()

        # Test sum
        sum = smooth.sum()
        self.test_value(sum, ref)

        # Store in sky map
        map = gammalib.GSkyMap('CAR','CEL',0.0,0.0,-1.0,1.0,10,10)
        for iy in range(10):
            for ix in range(10):
                map[ix+iy*10] = smooth[ix,iy]
        map.save('test_fft.fits', True)
        
        # Return
        return
