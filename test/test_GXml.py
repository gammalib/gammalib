# ==========================================================================
# This module performs unit tests for the GammaLib xml module.
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
import test_support


# ================================== #
# Test class for GammaLib xml module #
# ================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib xml module
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

    # Setup GXml container
    def _setup_xml(self):
        """
        Setup GXml container

        Returns
        -------
        xml : `~gammalib.GXml`
            GXml container
        """
        # Setup XML container
        xml     = gammalib.GXml()
        element = gammalib.GXmlElement()
        for i in range(10):
            element.name('%s' % i)
            xml.append(element)

        # Return XML container
        return xml

    # Setup GXmlElement container
    def _setup_elements(self):
        """
        Setup GXmlElement container

        Returns
        -------
        elements : `~gammalib.GXmlElement`
            GXml container
        """
        # Setup GXmlElement container
        elements = gammalib.GXmlElement()
        element  = gammalib.GXmlElement()
        for i in range(10):
            element.name('%s' % i)
            elements.append(element)

        # Return GXmlElement container
        return elements

    # Test GXml class access operators
    def _test_xml_access(self):
        """
        Test GXml class parameter access
        """
        # Setup XML container and element
        xml     = self._setup_xml()
        element = gammalib.GXmlElement()

        # Perform GXml access tests
        test_support._container_access_index(self, xml)

        # Check parameter setting by index from start
        element.name('98')
        xml[3] = element
        self.test_value(xml[3].name(), '98')

        # Check parameter setting by index from end
        element.name('99')
        xml[-2] = element
        self.test_value(xml[-2].name(), '99')

        # Return
        return

    # Test GXml class slicing
    def _test_xml_slicing(self):
        """
        Test GXml class slicing
        """
        # Setup XML container
        xml = self._setup_xml()

        # Perform slicing tests
        test_support._container_slicing(self, xml)

        # Return
        return

    # Test GXmlElement class access operators
    def _test_elements_access(self):
        """
        Test GXmlElement class parameter access
        """
        # Setup GXmlElement container and element
        elements = self._setup_elements()
        element  = gammalib.GXmlElement()

        # Perform GXmlElement access tests
        test_support._container_access_index(self, elements)

        # Check parameter setting by index from start
        element.name('98')
        elements[3] = element
        self.test_value(elements[3].name(), '98')

        # Check parameter setting by index from end
        element.name('99')
        elements[-2] = element
        self.test_value(elements[-2].name(), '99')

        # Return
        return

    # Test GXmlElement class slicing
    def _test_elements_slicing(self):
        """
        Test GXmlElement class slicing
        """
        # Setup GXmlElement container
        elements = self._setup_elements()

        # Perform slicing tests
        test_support._container_slicing(self, elements)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support._pickeling(self, gammalib.GXml())
        test_support._pickeling(self, gammalib.GXmlAttribute())
        test_support._pickeling(self, gammalib.GXmlComment())
        test_support._pickeling(self, gammalib.GXmlDocument())
        test_support._pickeling(self, gammalib.GXmlElement())
        test_support._pickeling(self, gammalib.GXmlPI())
        test_support._pickeling(self, gammalib.GXmlText())

        # Perform pickeling tests of filled classes
        test_support._pickeling(self, gammalib.GXml('<?xml version="1.0" standalone="no"?><a>Text</a>'))
        test_support._pickeling(self, gammalib.GXmlAttribute('name','value'))
        test_support._pickeling(self, gammalib.GXmlComment('comment'))
        test_support._pickeling(self, gammalib.GXmlDocument('test.xml','a','b','c'))
        test_support._pickeling(self, gammalib.GXmlElement('a file="x"'))
        test_support._pickeling(self, gammalib.GXmlPI('<?PI?>'))
        test_support._pickeling(self, gammalib.GXmlText('Text'))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('xml')

        # Append tests
        self.append(self._test_xml_access, 'Test GXml parameter access')
        self.append(self._test_xml_slicing, 'Test GXml slicing')
        self.append(self._test_elements_access, 'Test GXmlElement parameter access')
        self.append(self._test_elements_slicing, 'Test GXmlElement slicing')
        self.append(self._test_pickeling, 'Test Xml class pickeling')

        # Return
        return
