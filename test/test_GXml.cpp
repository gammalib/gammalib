/***************************************************************************
 *                    test_GXml.cpp - Test xml module                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2022 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file test_GXml.cpp
 * @brief Implementation of unit tests for XML module
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <cstdlib>     // getenv
#include "test_GXml.hpp"
#include "GTools.hpp"

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGXml::set(void)
{
    // Test name
    name("GXml");

    // Set datadir
    std::string datadir = gammalib::getenv("TEST_DATA");

    // Set XML filename
    m_xml_file = datadir + "/test.xml";

    // Append tests
    append(static_cast<pfunction>(&TestGXml::test_GXml_attributes),
           "Test XML attributes");
    append(static_cast<pfunction>(&TestGXml::test_GXml_elements),
           "Test XML elements");
    append(static_cast<pfunction>(&TestGXml::test_GXml_construct),
           "Test XML constructors");
    append(static_cast<pfunction>(&TestGXml::test_GXml_load),
           "Test XML load");
    append(static_cast<pfunction>(&TestGXml::test_GXml_access),
           "Test XML access");

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGXml* TestGXml::clone(void) const
{
    // Clone test suite
    return new TestGXml(*this);
}


/***********************************************************************//**
 * @brief Test XML attributes
 **************************************************************************/
void TestGXml::test_GXml_attributes(void)
{
    // Test attribute name-value constructors
    GXmlAttribute attr1("test", "1.0");
    test_assert(attr1.value() == "1.0", "Test if value()= 1.0",
                "Unexpected attribute "+attr1.value());
    GXmlAttribute attr2("test", "\"1.0\"");
    test_assert(attr2.value() == "1.0", "Test if value()= 1.0",
                "Unexpected attribute "+attr2.value());
    GXmlAttribute attr3("test", "'1.0'");
    test_assert(attr3.value() == "1.0", "Test if value()= 1.0",
                "Unexpected attribute "+attr3.value());
    GXmlAttribute attr4("test", "''1.0'");
    test_assert(attr4.value() == "'1.0", "Test if value()= '1.0",
                "Unexpected attribute "+attr4.value());
    GXmlAttribute attr5("test", "'1.0");
    test_assert(attr5.value() == "'1.0", "Test if value()= '1.0",
                "Unexpected attribute "+attr5.value());
    GXmlAttribute attr6("test", "1.0'");
    test_assert(attr6.value() == "1.0'", "Test if value()= 1.0'",
                "Unexpected attribute "+attr6.value());
    GXmlAttribute attr7("test", "\"1.0");
    test_assert(attr7.value() == "\"1.0", "Test if value()= \"1.0",
                "Unexpected attribute "+attr7.value());
    GXmlAttribute attr8("test", "1.0\"");
    test_assert(attr8.value() == "1.0\"", "Test if value()= 1.0\"",
                "Unexpected attribute "+attr8.value());
    GXmlAttribute attr9("test", "\"\"1.0\"");
    test_assert(attr9.value() == "\"1.0", "Test if value()= \"1.0",
                "Unexpected attribute "+attr9.value());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test XML elements
 **************************************************************************/
void TestGXml::test_GXml_elements(void)
{
    // Test void constructor
    GXmlElement element1;
    test_value(element1.name(), "", "Check empty name of void constructor");

    // Test integer value constructor
    GXmlElement element2("integer", 41);
    test_value(element2.name(), "integer", "Check name of integer constructor");
    test_value(element2.value(), "41", "Check string value of integer constructor");
    test_value(element2.integer(), 41, "Check integer value of integer constructor");

    // Test floating point value constructor
    GXmlElement element3("double", 41.123);
    test_value(element3.name(), "double", "Check name of double constructor");
    test_value(element3.value(), "41.123", "Check string value of double constructor");
    test_value(element3.real(), 41.123, "Check double value of double constructor");

    // Test string value constructor
    GXmlElement element4("string", "41");
    test_value(element4.name(), "string", "Check name of string constructor");
    test_value(element4.value(), "41", "Check string value of string constructor");
    test_value(element4.real(), 41.0, "Check double value of string constructor");
    test_value(element4.integer(), 41, "Check integer value of string constructor");

    // Perform element attributes
    GXmlElement element;
    element.attribute("test", "1.0");
    test_assert(element.attribute("test") == "1.0",
                "Test if element.attribute(\"test\")= 1.0",
                "Unexpected attribute "+element.attribute("test"));

    element.attribute("test", "2.0");
    test_assert(element.attribute("test") == "2.0",
                "Test if element.attribute(\"test\")= 2.0",
                "Unexpected attribute "+element.attribute("test"));

    element.attribute("test2", "1.0");
    test_assert(element.attribute("test2") == "1.0",
                "Test if element.attribute(\"test2\")= 1.0",
                "Unexpected attribute "+element.attribute("test2"));

    test_assert(element.attribute("test")  == "2.0" &&
                element.attribute("test2") == "1.0",
                "Test if  element.attribute(\"test\")= 2.0 and if element.attribute(\"test2\")= 1.0",
                "Unexpected attributes "+element.attribute("test")+" "+element.attribute("test2"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test XML constructors
 **************************************************************************/
void TestGXml::test_GXml_construct(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GXml xml;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test load constructor
    test_try("Test load constructor");
    try {
        GXml xml(m_xml_file);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test copy constructor
    test_try("Test copy constructor");
    try {
        GXml xml1(m_xml_file);
        GXml xml2 = xml1;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML file creation
    test_try("Test XML file creation");
    try {
        GXml xml;
        xml.append(GXmlComment("This is a comment."));
        xml.append(GXmlElement("source_library"));
        GXmlElement* lib = xml.element("source_library", 0);
        lib->append(GXmlElement("source name=\"LMC\" type=\"DiffuseSource\""));
        GXmlNode* src = lib->element("source", 0);
        src->append(GXmlElement("spectrum type=\"PLSuperExpCutoff\""));
        GXmlNode* spec = src->element("spectrum", 0);
        spec->append(GXmlElement("parameter free=\"1\" max=\"1000\" min=\"1e-07\""
                " name=\"Prefactor\" scale=\"1e-07\""
                " value=\"0.02754520844\""));
        spec->append(GXmlElement("parameter free=\"1\" max=\"5\" min=\"-5\""
                " name=\"Index1\" scale=\"1\" value=\"-2.0458781\""));
        GXmlElement* par = spec->element("parameter", 0);
        par->attribute("value", "1.01");
        par->attribute("error", "3.145");
        par = spec->element("parameter", 1);
        par->attribute("value", "-2.100");
        par->attribute("error", "9.876");
        src->append(GXmlElement("spatialModel file=\"LMC.fits\" type=\"SpatialMap\""));
        GXmlNode* spat = src->element("spatialModel", 0);
        spat->append(GXmlElement("parameter free=\"0\" max=\"1000\" min=\"0.001\""
                " name=\"Prefactor\" scale=\"1\" value=\"1\""));
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test XML loading/saving
 **************************************************************************/
void TestGXml::test_GXml_load(void) 
{
    // Test loading
    test_try("Test loading");
    try {
        GXml xml;
        xml.load(m_xml_file);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test saving
    test_try("Test saving");
    try {
        GXml xml;
        xml.load(m_xml_file);
        xml.save("test.xml");
        xml.load("test.xml");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test loading of saved XML document
    test_try("Test loading of saved XML document");
    try {
        GXml xml;
        xml.load("test.xml");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test XML element access
 **************************************************************************/
void TestGXml::test_GXml_access(void)
{
    // Test root document access
    GXml xml;
    xml.load(m_xml_file);
    test_value(xml.size(), 3, "Check number of child elements");

    // Test node access
    for (int i = 0; i < xml.size(); ++i) {
        GXmlNode* ptr = xml[i];
        test_assert(ptr != 0, "Check XML node access");
    }
    test_value(xml.elements(), 1, "Check number of child elements");

    // Test node access
    for (int i = 0; i < xml.elements(); ++i) {
        GXmlNode* ptr = xml.element(i);
        test_assert(ptr != 0, "Check XML element access");
    }
    test_value(xml.elements("source_library"), 1, "Check number of child elements");

    // Test element access
    for (int i = 0; i < xml.elements("source_library"); ++i) {
        GXmlElement* ptr = xml.element("source_library", i);
        test_value(ptr->name(), "source_library", "Check element name");
    }

    // Test hierarchy access
    GXmlElement* ptr = NULL;
    ptr = xml.element("source_library");
    test_value(ptr->name(), "source_library", "Check hierarchy level 1");
    ptr = xml.element("source_library > source");
    test_value(ptr->name(),  "source", "Check hierarchy level 2");
    ptr = xml.element("source_library > source > spectrum");
    test_value(ptr->name(), "spectrum", "Check hierarchy level 3");
    ptr = xml.element("source_library > source > spectrum > parameter[2]");
    test_value(ptr->name(), "parameter", "Check hierarchy level 4");
    test_value(ptr->attribute("name"), "PivotEnergy", "Check for attribute");

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("XML module");

    // Create a test suite
    TestGXml test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    bool success = testsuites.run();

    // Save test report
    testsuites.save("reports/GXml.xml");

    // Return success status
    return (success ? 0 : 1);
}
