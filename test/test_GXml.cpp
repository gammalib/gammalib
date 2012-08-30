/***************************************************************************
 *                  test_GXml.cpp  -  test GXml classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
 * @brief Testing of XML module.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GammaLib.hpp"
#include "test_GXml.hpp"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

void TestGXml::set(void)
{
    // Test name
    name("GXml");

    // xml file
    m_xml_file="data/test.xml";

    //Add Tests
    add_test(static_cast<pfunction>(&TestGXml::test_GXml_attributes),"xml attributes");
    add_test(static_cast<pfunction>(&TestGXml::test_GXml_elements),"xml elements");
    add_test(static_cast<pfunction>(&TestGXml::test_GXml_construct),"xml constructors");
    add_test(static_cast<pfunction>(&TestGXml::test_GXml_load),"xml load");
    add_test(static_cast<pfunction>(&TestGXml::test_GXml_access),"xml access");

    return; 
}

/***************************************************************************
 * Test: GXml element attributes                                           *
 ***************************************************************************/
void TestGXml::test_GXml_attributes(void)
{
    // Test valid attributes
    test_try("Test valid attributes");
    try {
        GXmlAttribute attr("test", "1.0");
        test_assert(attr.value() == "1.0","Test if value()= 1.0","Unexpected attribute "+attr.value());

        attr.value("\"1.0\"");
        test_assert(attr.value() == "1.0","Test if value()= 1.0","Unexpected attribute "+attr.value());

        attr.value("'1.0'");
        test_assert(attr.value() == "1.0","Test if value()= 1.0","Unexpected attribute "+attr.value());

        attr.value("''1.0'");
        test_assert(attr.value() == "''1.0'",
                              "Test if value()= ''1.0'",
                              "Unexpected attribute "+attr.value());

        attr.value("'1.0");
        test_assert(attr.value() == "'1.0",
                              "Test if value()= '1.0",
                              "Unexpected attribute "+attr.value());

        attr.value("1.0'");
        test_assert(attr.value() == "1.0'",
                              "Test if value()= 1.0'",
                              "Unexpected attribute "+attr.value());

        attr.value("\"1.0");
        test_assert(attr.value() == "\"1.0",
                              "Test if value()= \"1.0",
                              "Unexpected attribute "+attr.value());


        attr.value("1.0\"");
        test_assert(attr.value() == "1.0\"",
                              "Test if value()= 1.0\"",
                              "Unexpected attribute "+attr.value());

        attr.value("\"\"1.0\"");
        test_assert(attr.value() == "\"\"1.0\"",
                              "Test if value()= \"\"1.0\"",
                              "Unexpected attribute "+attr.value());
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test invalid attributes
    test_try("Test invalid attributes 1");
    try {
        GXmlAttribute attr("test", "\"\"1.0'\"");
        test_try_failure();
    }
    catch (GException::xml_attribute_value &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    test_try("Test invalid attributes 2");
    try {
        GXmlAttribute attr("test", "''1.0\"");
        test_try_failure();
    }
    catch (GException::xml_attribute_value &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    test_try("Test invalid attributes 3");
    try {
        GXmlAttribute attr("test", "\"1.0'");
        test_try_failure();
    }
    catch (GException::xml_attribute_value &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    return;
}


/***************************************************************************
 * Test: GXml elements                                                     *
 ***************************************************************************/
void TestGXml::test_GXml_elements(void)
{
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

    test_assert(element.attribute("test")  == "2.0"&& element.attribute("test2") == "1.0",
                            "Test if  element.attribute(\"test\")= 2.0 and if element.attribute(\"test2\")= 1.0",
                            "Unexpected attributes "+element.attribute("test")+" "+element.attribute("test2"));


    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml constructors                                                 *
 ***************************************************************************/
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
        xml.append(new GXmlComment("This is a comment."));
        xml.append(new GXmlElement("source_library"));
        GXmlElement* lib = xml.element("source_library", 0);
        lib->append(new GXmlElement("source name=\"LMC\" type=\"DiffuseSource\""));
        GXmlNode* src = lib->element("source", 0);
        src->append(new GXmlElement("spectrum type=\"PLSuperExpCutoff\""));
        GXmlNode* spec = src->element("spectrum", 0);
        spec->append(new GXmlElement("parameter free=\"1\" max=\"1000\" min=\"1e-07\""
                " name=\"Prefactor\" scale=\"1e-07\""
                " value=\"0.02754520844\""));
        spec->append(new GXmlElement("parameter free=\"1\" max=\"5\" min=\"-5\""
                " name=\"Index1\" scale=\"1\" value=\"-2.0458781\""));
        GXmlElement* par = (GXmlElement*)spec->element("parameter", 0);
        par->attribute("value", "1.01");
        par->attribute("error", "3.145");
        par = (GXmlElement*)spec->element("parameter", 1);
        par->attribute("value", "-2.100");
        par->attribute("error", "9.876");
        src->append(new GXmlElement("spatialModel file=\"LMC.fits\" type=\"SpatialMap\""));
        GXmlNode* spat = src->element("spatialModel", 0);
        spat->append(new GXmlElement("parameter free=\"0\" max=\"1000\" min=\"0.001\""
                " name=\"Prefactor\" scale=\"1\" value=\"1\""));
        //xml.save("test2.xml");
        //std::cout << xml << std::endl;

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml loading/saving                                               *
 ***************************************************************************/
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

    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml element access                                               *
 ***************************************************************************/
void TestGXml::test_GXml_access(void)
{
    // Test root document access
    GXml xml;
    xml.load(m_xml_file);

    test_assert(xml.children() == 3,
                            "Test if xml.children()==3",
                            "Unexpected number of children in document "+str(xml.children()));

    for (int i = 0; i < xml.children(); ++i)
        GXmlNode* ptr = xml.child(i);

    test_assert(xml.elements() == 1,
                            "Test if xml.elements()==1",
                            "Unexpected number of child elements in document "+str(xml.elements()));

    for (int i = 0; i < xml.elements(); ++i)
        GXmlNode* ptr = xml.element(i);

    test_assert( xml.elements("source_library") == 1,
                            "Test if the source_library = 1",
                            "Unexpected number of child elements in document "+str(xml.elements("source_library")));

    for (int i = 0; i < xml.elements("source_library"); ++i) {
        GXmlElement* ptr = xml.element("source_library", i);
        test_assert(ptr->name() == "source_library","Test name","Unexpected element name "+ptr->name());
    }
    //std::cout << xml.elements("source_library") << std::endl;


    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GXml");

    bool was_successful=true;

    //Create a test suite
    TestGXml test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GXml.xml");

    // Return
    return was_successful ? 0:1;
}
