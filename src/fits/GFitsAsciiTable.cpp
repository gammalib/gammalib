/***************************************************************************
 *              GFitsAsciiTable.cpp - FITS ASCII table class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GFitsAsciiTable.cpp
 * @brief FITS ASCII table class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GFitsAsciiTable.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable(void) : GFitsTable()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Table constructor
 *
 * @param[in] nrows Number of rows in table.
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable(const int& nrows) : GFitsTable(nrows)
{
    // Initialise class members for clean destruction
    init_members();

    // Initialise header
    init_table_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] table ASCII table.
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable(const GFitsAsciiTable& table) : GFitsTable(table)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsAsciiTable::~GFitsAsciiTable()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] table ASCII table.
 * @return ASCII table.
 ***************************************************************************/
GFitsAsciiTable& GFitsAsciiTable::operator=(const GFitsAsciiTable& table)
{
    // Execute only if object is not identical
    if (this != &table) {
  
        // Copy base class members
        this->GFitsTable::operator=(table);

        // Free members
        free_members();
  
        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(table);
	
    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear ASCII table
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GFitsAsciiTable::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GFitsTable::free_members();
    this->GFitsHDU::free_members();

    // Initialise members
    this->GFitsHDU::init_members();
    this->GFitsTable::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone ASCII table
 *
 * @return Pointer to deep copy of ASCII table.
 ***************************************************************************/
GFitsAsciiTable* GFitsAsciiTable::clone(void) const 
{
    return new GFitsAsciiTable(*this);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsAsciiTable::init_members(void)
{
    // Initialise members
    m_type = 1;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param table Table.
 ***************************************************************************/
void GFitsAsciiTable::copy_members(const GFitsAsciiTable& table)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsAsciiTable::free_members(void)
{
    // Return
    return;
}




/***********************************************************************//**
 * @brief Initialise ASCII table header
 *
 * Initialises the table header by setting the default header cards.
 ***************************************************************************/
void GFitsAsciiTable::init_table_header(void)
{
    // Compute total width in Bytes
    int width = 0;
    for (int i = 0; i < ncols(); ++i) {
        width += m_columns[i]->width();
    }

    // Set image header keywords
    m_header.append(GFitsHeaderCard("XTENSION", "TABLE",
                                    "ASCII table extension"));
    m_header.append(GFitsHeaderCard("BITPIX", 8,
                                    "8-bit ASCII characters"));
    m_header.append(GFitsHeaderCard("NAXIS", 2,
                                    "2-dimensional ASCII table"));
    m_header.append(GFitsHeaderCard("NAXIS1", width,
                                    "width of table in characters"));
    m_header.append(GFitsHeaderCard("NAXIS2", nrows(),
                                    "number of rows in table"));
    m_header.append(GFitsHeaderCard("PCOUNT", 0,
                                    "no group parameters (required keyword)"));
    m_header.append(GFitsHeaderCard("GCOUNT", 1,
                                    "one data group (required keyword)"));

    // Return
    return;
}
