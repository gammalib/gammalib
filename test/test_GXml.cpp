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
    test_GXml_construct();
    test_GXml_load();

    // Return
    return 0;
}
