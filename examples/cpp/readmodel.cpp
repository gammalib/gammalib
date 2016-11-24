/***************************************************************************
 *             readmodel.cpp - Illustrates how to read a model             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Juergen Knoedlseder                         *
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
 * @file readmodel.cpp
 * @brief Illustrates how to read a model
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <fstream>
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Create XML file
 *
 * Creates XML file that should be read back to avoid shipping of a test
 * XML file.
 ***************************************************************************/
void create_xml_file(void)
{
    // Open ASCII file
    std::ofstream file("example_source.xml");

    // Write file if file is open
    if (file.is_open()) {

        // Write XML file
        file << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
        file << "<source_library title=\"source library\">\n";
        file << "  <source name=\"1FGL J0005.7+3815\" type=\"PointSource\">\n";
        file << "    <spectrum type=\"PowerLaw\">\n";
        file << "      <parameter scale=\"1e-07\" name=\"Prefactor\"";
        file << " min=\"1e-07\" max=\"1000.0\" value=\"1.73\" free=\"1\"/>\n";
        file << "      <parameter scale=\"1.0\" name=\"Index\" min=\"-5.0\"";
        file << " max=\"+5.0\" value=\"-2.1\" free=\"1\"/>\n";
        file << "    <parameter scale=\"1.0\" name=\"PivotEnergy\" min=\"10.0\"";
        file << " max=\"1000000.0\" value=\"100.0\" free=\"0\"/>\n";
        file << "    </spectrum>\n";
        file << "    <spatialModel type=\"PointSource\">\n";
        file << "      <parameter free=\"0\" max=\"360\" min=\"-360\"";
        file << " name=\"RA\" scale=\"1\" value=\"83.6331\"/>\n";
        file << "      <parameter free=\"0\" max=\"90\" min=\"-90\"";
        file << " name=\"DEC\" scale=\"1\" value=\"22.0145\" />\n";
        file << "    </spatialModel>\n";
        file << "  </source>\n";
        file << "</source_library>\n";

        // Close file
        file.close();
        
    } // endif: file was open

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read XML model
 ***************************************************************************/
int main(void)
{
    // Create XML file
    create_xml_file();

    // Read XML model
    GModels models("example_source.xml");
    std::cout << models << std::endl;

    // Return
    return 0;
}

