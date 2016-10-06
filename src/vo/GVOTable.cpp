/***************************************************************************
 *                      GVOTable.cpp - VO table class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Thierry Louge                               *
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
 * @file GVOTable.cpp
 * @brief VO table class implementation
 * @author Thierry Louge
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GVOTable.hpp"
#include "GXmlElement.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GVOTable::GVOTable(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS table constructor
 *
 * @param[in] table FITS table
 ***************************************************************************/
GVOTable::GVOTable(const GFitsTable& table)
{
    // Initialise members
    init_members();
    
    // Read VO table from FITS table
    read(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] votable VO table
 ***************************************************************************/
GVOTable::GVOTable(const GVOTable& votable)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(votable);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GVOTable::~GVOTable(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] VO table.
 * @return VO table.
 ***************************************************************************/
GVOTable& GVOTable::operator=(const GVOTable& votable)
{
    // Execute only if object is not identical
    if (this != &votable) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(votable);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object.
 *
 * Reset object to a clean initial state.
 ***************************************************************************/
void GVOTable::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GVOTable* GVOTable::clone(void) const
{
    // Clone client
    return new GVOTable(*this);
}


/***********************************************************************//**
 * @brief Read VO table from FITS table
 ***************************************************************************/
void GVOTable::read(const GFitsTable& table)
{
    // Clear VO table
    m_xml.clear();

    // Set name and resource string
    m_name     = table.extname();
    m_resource = table.classname();

    // Set VOTABLE element string
    std::string s1("VOTABLE version=\"1.3\" xmlns:xsi=\"http://www.w3.org/"
                   "2001/XMLSchema-instance\" xmlns=\"http://www.ivoa.net/"
                   "xml/VOTable/v1.3\" xmlns:stc=\"http://www.ivoa.net/xml/"
                   "STC/v1.30\"");
    std::string s2("RESOURCE name=\""+m_resource+"\"");
    std::string s3("TABLE name=\""+m_name+"\"");

    // Create table element
    GXmlElement* votable = m_xml.append(s1)->append(s2)->append(s3);

    // Append description
    votable->append("DESCRIPTION")->append(GXmlText(m_description));

    // Append FIELD elements by extracting the information from the FITS
    // columns
    for (int i = 0; i < table.ncols(); ++i) {
        votable->append(field_from_fits_column(*table[i]));
    }

    // Append data
    votable->append(data_from_fits_table(table));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print VO Table information
 *
 * @param[in] chatter Chattiness.
 * @return String containing VO table information.
 ***************************************************************************/
std::string GVOTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append VO table
        result.append("=== GVOTable ===");
        result.append(m_xml.print(chatter, 0));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GVOTable::init_members(void)
{
    // Initialise members
    m_xml.clear();
    m_name.clear();
    m_resource.clear();
    m_description.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] table VO table.
 ***************************************************************************/
void GVOTable::copy_members(const GVOTable& table)
{
    // Copy members
    m_xml         = table.m_xml;
    m_name        = table.m_name;
    m_resource    = table.m_resource;
    m_description = table.m_description;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVOTable::free_members(void)
{ 
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return FIELD element with column description
 *
 * @param[in] column FITS column
 * @return XML FIELD element with column description
 ***************************************************************************/
GXmlElement GVOTable::field_from_fits_column(const GFitsTableCol& column) const
{
    // Create FIELD element
    GXmlElement field("FIELD");

    // Set datatype
    std::string datatype;
    switch (column.type()) {
        case __TLOGICAL:
            datatype = "boolean";
            break;
        case __TBIT:
            datatype = "bit";
            break;
        case __TBYTE:
        case __TSBYTE:
            datatype = "unsignedByte";
            break;
        case __TSHORT:
        case __TUSHORT:
            datatype = "short";
            break;
        case __TINT:
        case __TUINT:
            datatype = "int";
            break;
        case __TLONG:
        case __TULONG:
        case __TLONGLONG:
            datatype = "long";
            break;
        case __TSTRING:
            datatype = "char";
            break;
        case __TFLOAT:
            datatype = "float";
            break;
        case __TDOUBLE:
            datatype = "double";
            break;
        case __TCOMPLEX:
            datatype = "floatComplex";
            break;
        case __TDBLCOMPLEX:
            datatype = "doubleComplex";
            break;
        default:
            datatype = "unknown";
            break;
    }

    // Set FIELD attributes
    field.attribute("name", column.name());
    field.attribute("ID", "col"+gammalib::str(column.colnum()+1));
    field.attribute("datatype", datatype);
    field.attribute("unit", column.unit());

    // Return XML FIELD element
    return field;
}


/***********************************************************************//**
 * @brief Return DATA element with FITS table data
 *
 * @param[in] table FITS table
 * @return XML DATA element with FITS data
 ***************************************************************************/
GXmlElement GVOTable::data_from_fits_table(const GFitsTable& table) const
{
    // Create DATA element and append TABLEDATA element
    GXmlElement  data("DATA");
    GXmlElement* tabledata = data.append("TABLEDATA");

    // Loop over all rows
    for (int row = 0; row < table.nrows(); ++row) {

        // Create TR element
        GXmlElement tr("TR");

        // Append all column data
        for (int col = 0; col < table.ncols(); ++col) {

            // Create TD element
            GXmlElement td("TD");

            // Set column data as text of TD element
            td.append(GXmlText(table[col]->string(row)));

            // Append TD element to TR element
            tr.append(td);

        } // endfor: looped over columns

        // Append TR element to TABLEDATA
        tabledata->append(tr);

    } // endfor: looped over rows

    // Return XML DATA element
    return data;
}
