/***************************************************************************
 *                  test_GXml.cpp  -  test GXml classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

/* __ Globals ____________________________________________________________ */
const std::string xml_file = "data/test.xml";


/***************************************************************************
 * Test: GXml element attributes                                           *
 ***************************************************************************/
void test_GXml_attributes(void)
{
    // Dump header
    std::cout << "Test XML attributes: ";

    // Test valid attributes
    try {
        GXmlAttribute attr("test", "1.0");
        if (attr.value() != "\"1.0\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("\"1.0\"");
        if (attr.value() != "\"1.0\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("'1.0'");
        if (attr.value() != "'1.0'") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("''1.0'");
        if (attr.value() != "\"''1.0'\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("'1.0");
        if (attr.value() != "\"'1.0\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("1.0'");
        if (attr.value() != "\"1.0'\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("\"1.0");
        if (attr.value() != "'\"1.0'") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("1.0\"");
        if (attr.value() != "'1.0\"'") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
        attr.value("\"\"1.0\"");
        if (attr.value() != "'\"\"1.0\"'") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << attr.value() << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to assign valid XML attributes."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test invalid attributes
    try {
        GXmlAttribute attr("test", "\"\"1.0'\"");
    }
    catch (GException::xml_attribute_value &e) {
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to catch invalid XML attribute."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    try {
        GXmlAttribute attr("test", "''1.0\"");
    }
    catch (GException::xml_attribute_value &e) {
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to catch invalid XML attribute."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    try {
        GXmlAttribute attr("test", "\"1.0'");
    }
    catch (GException::xml_attribute_value &e) {
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to catch invalid XML attribute."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml elements                                                     *
 ***************************************************************************/
void test_GXml_elements(void)
{
    // Dump header
    std::cout << "Test XML elements: ";

    // Test attribute setting
    try {
        GXmlElement   element;
        element.attribute("test", "1.0");
        if (element.attribute("test") != "\"1.0\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << element.attribute("test") << std::endl;
            throw;
        }
        std::cout << ".";
        element.attribute("test", "2.0");
        if (element.attribute("test") != "\"2.0\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attribute "
                      << element.attribute("test") << std::endl;
            throw;
        }
        std::cout << ".";
        element.attribute("test2", "1.0");
        if (element.attribute("test")  != "\"2.0\"" ||
            element.attribute("test2") != "\"1.0\"") {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected attributes "
                      << element.attribute("test")
                      << " " << element.attribute("test2") << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to assign XML attributes."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml constructors                                                 *
 ***************************************************************************/
void test_GXml_construct(void)
{
    // Dump header
    std::cout << "Test XML constructors: ";

    // Test void constructor
    try {
        GXml xml;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to construct empty XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test load constructor
    try {
        GXml xml(xml_file);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to construct empty XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test copy constructor
    try {
        GXml xml1(xml_file);
        GXml xml2 = xml1;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to copy construct XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml loading/saving                                               *
 ***************************************************************************/
void test_GXml_load(void)
{
    // Dump header
    std::cout << "Test XML loading and saving: ";

    // Test loading
    try {
        GXml xml;
        xml.load(xml_file);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to load XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test saving
    try {
        GXml xml;
        xml.load(xml_file);
        xml.save("test.xml");
        xml.load("test.xml");
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to save XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test loading of saved XML document
    try {
        GXml xml;
        xml.load("test.xml");
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to load previously saved XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 * Test: GXml element access                                               *
 ***************************************************************************/
void test_GXml_access(void)
{
    // Dump header
    std::cout << "Test XML element access: ";

    // Test root document access
    try {
        GXml xml;
        xml.load(xml_file);
        std::cout << ".";
        if (xml.children() != 3) {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected number of children in document "
                      << xml.children() << std::endl;
            throw;
        }
        std::cout << ".";
        for (int i = 0; i < xml.children(); ++i)
            GXmlNode* ptr = xml.child(i);
        std::cout << ".";
        if (xml.elements() != 1) {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected number of child elements in document "
                      << xml.elements() << std::endl;
            throw;
        }
        std::cout << ".";
        for (int i = 0; i < xml.elements(); ++i)
            GXmlNode* ptr = xml.element(i);
        std::cout << ".";
        if (xml.elements("source_library") != 1) {
            std::cout << std::endl
                      << "TEST ERROR: Unexpected number of child elements in document "
                      << xml.elements("source_library") << std::endl;
            throw;
        }
        std::cout << ".";
        for (int i = 0; i < xml.elements("source_library"); ++i) {
            GXmlElement* ptr = xml.element("source_library", i);
            if (ptr->name() != "source_library") {
                std::cout << std::endl
                        << "TEST ERROR: Unexpected element name "
                        << ptr->name() << std::endl;
                throw;
            }
        }
        std::cout << ".";
        //std::cout << xml.elements("source_library") << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to access XML elements."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "************************" << std::endl;
    std::cout << "* GXml classes testing *" << std::endl;
    std::cout << "************************" << std::endl;

    // Execute XML tests
    test_GXml_attributes();
    test_GXml_elements();
    test_GXml_construct();
    test_GXml_load();
    test_GXml_access();

    // Return
    return 0;
}
