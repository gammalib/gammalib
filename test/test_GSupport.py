# ==========================================================================
# This module performs unit tests for the GammaLib support module.
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
import math


# =============================== #
# Test class for GammaLib support #
# =============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib support.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set all test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('support')

        # Append tests
        self.append(self.test_node_array, 'Test GNodeArray')
        self.append(self.test_url_file,   'Test GUrlFile')
        self.append(self.test_url_string, 'Test GUrlString')
        self.append(self.test_filename,   'Test GFilename')
        self.append(self.test_csv,        'Test GCsv')

        # Return
        return

    # Test GNodeArray class
    def test_node_array(self):
        """
        Test GNodeArray class.
        """
        # Set-up vector and data array. Test all vector elements.
        vector = gammalib.GVector(20)
        data   = gammalib.GVector(20)
        for i in range(20):
            vector[i] = 10.0 + i * 5.0
            data[i]   = math.sin(0.15 * (vector[i] - 10.0))
            self.test_value(data[i], math.sin(0.15 * i * 5.0))

        # Set-up node array
        array = gammalib.GNodeArray()
        array.nodes(vector)

        # Get values
        x_val = []
        y_val = []
        for i in range(100):
            x = i - 10
            array.set_value(x)
            inx_left  = array.inx_left()
            inx_right = array.inx_right()
            wgt_left  = array.wgt_left()
            wgt_right = array.wgt_right()
            y = wgt_left * data[inx_left] + wgt_right * data[inx_right]
            x_val.append(x)
            y_val.append(y)

        # Return
        return

    # Test GUrlFile class
    def test_url_file(self):
        """
        Test GUrlFile class.
        """
        # Test file writing
        url = gammalib.GUrlFile('test_url.dat', 'w')
        self.test_value(url.write('abcd', 4), 4)
        url.put_char(ord('e'))
        url.close()

        # Test file reading
        buffer = ''
        url    = gammalib.GUrlFile('test_url.dat', 'r')
        buffer = url.read(99)
        self.test_value(buffer, 'abcde')
        url.close()

        # Return
        return
 
    # Test GUrlString class
    def test_url_string(self):
        """
        Test GUrlString class.
        """
        # Test writing
        url = gammalib.GUrlString()
        self.test_value(url.write('abcd', 4), 4)
        url.put_char(ord('e'))

        # Test reading
        url.rewind()
        buffer = ''
        buffer = url.read(99)
        self.test_value(buffer, 'abcde')
        url.close()

        # Return
        return
       
    # Test GFilename class
    def test_filename(self):
        """
        Test GFilename class.
        """
        # Test assigning file names
        filename = gammalib.GFilename('myfile.fits')
        self.test_value(filename.url(), 'myfile.fits')
        #
        filename = gammalib.GFilename('myfile.fits[EVENTS]')
        self.test_value(filename.url(), 'myfile.fits')
        self.test_value(filename.extname(), 'EVENTS')
        #
        filename = gammalib.GFilename('myfile.fits[0]')
        self.test_value(filename.url(), 'myfile.fits')
        self.test_value(filename.extno(), 0)
        #
        filename = gammalib.GFilename('myfile.fits[0,2]')
        self.test_value(filename.url(), 'myfile.fits')
        self.test_value(filename.extno(), 0)
        self.test_value(filename.extver(), 2)
        
        # Return
        return

    # Test GCsv class
    def test_csv(self):
        """
        Test GCsv class.
        """
        # Create CSV file using the append method
        csv = gammalib.GCsv()
        csv.append(['ra','dec','duration'])
        csv.append(['20.0','-70.0','1800.0'])

        # Test created file
        self.test_value(csv.size(), 6)
        self.test_value(csv.ncols(), 3)
        self.test_value(csv.nrows(), 2)
        self.test_value(csv[0,0], 'ra')
        self.test_value(csv[0,1], 'dec')
        self.test_value(csv[0,2], 'duration')
        self.test_value(csv[1,0], '20.0')
        self.test_value(csv[1,1], '-70.0')
        self.test_value(csv[1,2], '1800.0')
        
        # Test saving and loading of file
        csv.save('test_csv_py.dat', ',', True)
        csv.load('test_csv_py.dat', ',')

        # Test loaded file
        self.test_value(csv.size(), 6)
        self.test_value(csv.ncols(), 3)
        self.test_value(csv.nrows(), 2)
        self.test_value(csv[0,0], 'ra')
        self.test_value(csv[0,1], 'dec')
        self.test_value(csv[0,2], 'duration')
        self.test_value(csv[1,0], '20.0')
        self.test_value(csv[1,1], '-70.0')
        self.test_value(csv[1,2], '1800.0')

        # Return
        return
