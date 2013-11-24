/***************************************************************************
 *                    createxml.cpp - Create XML file                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file createxml.cpp
 * @brief Create XML file
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Create XML file
 *
 * This code illustrates the creation of a XML file.
 ***************************************************************************/
int main(void) {

    // Allocate XML object
    GXml xml;

    // Create base elements
    GXmlComment comment("Now 2 nodes with spatial and spectral info");
    GXmlElement spatial("spatial type=\"Position\"");
    GXmlElement spectral("spectrum type=\"PowerLaw\"");
    GXmlElement text("text");
    GXmlPI      pi("<?process now?>");

    // Append spatial parameters
    spatial.append(GXmlElement("parameter ra=\"83.0\""));
    spatial.append(GXmlElement("parameter dec=\"22.0\""));

    // Append spectral parameters
    spectral.append(GXmlElement("parameter prefactor=\"1e-7\""));
    spectral.append(GXmlElement("parameter index=\"-2.1\""));

    // Append text
    text.append(GXmlText("Finish with text"));

    // Append elements
    xml.append(comment);
    xml.append(spatial);
    xml.append(spectral);
    xml.append(text);
    xml.append(pi);

    // Manipulate source element
    xml.element("spatial", 0)->append(GXmlElement("parameter ra=\"83.0\""));
    xml.element("spatial", 0)->append(GXmlElement("parameter dec=\"22.0\""));

    // Show alternative access
    GXmlNode* node = xml.append(spatial);
    node->append(GXmlElement("parameter ra=\"83.0\""));
    node->append(GXmlElement("parameter dec=\"22.0\""));

    // Save XML document to file
    xml.save("my_xml_file.xml");

    // Exit
    return 0;
}

