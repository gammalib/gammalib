# ==========================================================================
# This module performs unit tests for the GammaLib support module.
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
import gammalib
import math
import test_support


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

    # Setup GNodeArray container
    def _setup_nodes(self):
        """
        Setup GNodeArray

        Returns
        -------
        nodes : `~gammalib.GNodeArray`
            Node array
        """
        # Setup time container
        nodes = gammalib.GNodeArray()
        for i in range(10):
            nodes.append(float(i))

        # Return GNodeArray
        return nodes

    # Test GNodeArray class
    def _test_node_array(self):
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

    # Test GNodeArray access
    def _test_node_array_access(self):
        """
        Test GNodeArray access
        """
        # Setup nodes container
        nodes = self._setup_nodes()

        # Loop over all elements using the container iterator and count the
        # number of iterations
        nnodes    = 0
        reference = 0.0
        for node in nodes:
            self.test_value(node, reference)
            nnodes    += 1
            reference += 1.0

        # Check that looping was successful
        self.test_value(nodes.size(), 10, 'Check container size')
        self.test_value(nnodes, 10, 'Check container iterator')

        # Test access from start
        self.test_value(nodes[3], 3.0)

        # Check access from end
        self.test_value(nodes[-2], 8.0)

        # Return
        return

    # Test GNodeArray slicing
    def _test_node_array_slicing(self):
        """
        Test GNodeArray slicing
        """
        # Setup nodes container
        nodes = self._setup_nodes()

        # Test nodes[start:end]
        self.test_value(len(nodes[3:5]), 2)
        self.test_value(nodes[3:5][0], 3.0)
        self.test_value(nodes[3:5][1], 4.0)

        # Test nodes[start:]
        self.test_value(len(nodes[7:]), 3)
        self.test_value(nodes[7:][0], 7.0)
        self.test_value(nodes[7:][1], 8.0)
        self.test_value(nodes[7:][2], 9.0)

        # Test nodes[:end]
        self.test_value(len(nodes[:2]), 2)
        self.test_value(nodes[:2][0], 0.0)
        self.test_value(nodes[:2][1], 1.0)

        # Test nodes[:]
        self.test_value(len(nodes[:]), 10)
        for i in range(10):
            self.test_value(nodes[:][i], float(i))

        # Test nodes[start:end:step]
        self.test_value(len(nodes[3:7:2]), 2)
        self.test_value(nodes[3:7:2][0], 3.0)
        self.test_value(nodes[3:7:2][1], 5.0)

        # Test nodes[start:end:step]
        self.test_value(len(nodes[6:3:-2]), 2)
        self.test_value(nodes[6:3:-2][0], 6.0)
        self.test_value(nodes[6:3:-2][1], 4.0)

        # Test nodes[-start:]
        self.test_value(len(nodes[-2:]), 2)
        self.test_value(nodes[-2:][0], 8.0)
        self.test_value(nodes[-2:][1], 9.0)

        # Test nodes[:-end]
        self.test_value(len(nodes[:-7]), 3)
        self.test_value(nodes[:-7][0], 0.0)
        self.test_value(nodes[:-7][1], 1.0)
        self.test_value(nodes[:-7][2], 2.0)

        # Return
        return

    # Test GUrlFile class
    def _test_url_file(self):
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
    def _test_url_string(self):
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
    def _test_filename(self):
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
    def _test_csv(self):
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

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        #test_support._pickeling(self, gammalib.GBilinear())
        #test_support._pickeling(self, gammalib.GCsv())
        test_support._pickeling(self, gammalib.GFilename())
        #test_support._pickeling(self, gammalib.GNodeArray())
        #test_support._pickeling(self, gammalib.GRan())
        #test_support._pickeling(self, gammalib.GUrl())
        #test_support._pickeling(self, gammalib.GUrlFile())
        #test_support._pickeling(self, gammalib.GUrlString())

        # Perform pickeling tests of filled classes
        #test_support._pickeling(self, gammalib.GBilinear())
        #test_support._pickeling(self, gammalib.GCsv())
        test_support._pickeling(self, gammalib.GFilename('test.fits[URL]'))
        #test_support._pickeling(self, gammalib.GNodeArray())
        #test_support._pickeling(self, gammalib.GRan())
        #test_support._pickeling(self, gammalib.GUrl())
        #test_support._pickeling(self, gammalib.GUrlFile())
        #test_support._pickeling(self, gammalib.GUrlString())

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
        self.append(self._test_node_array, 'Test GNodeArray')
        self.append(self._test_node_array_access, 'Test GNodeArray node access')
        self.append(self._test_node_array_slicing, 'Test GNodeArray slicing')
        self.append(self._test_url_file,   'Test GUrlFile')
        self.append(self._test_url_string, 'Test GUrlString')
        self.append(self._test_filename,   'Test GFilename')
        self.append(self._test_csv,        'Test GCsv')
        self.append(self._test_pickeling, 'Test support class pickeling')

        # Return
        return

