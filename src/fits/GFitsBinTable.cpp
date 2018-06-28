/***************************************************************************
 *               GFitsBinTable.cpp - FITS binary table class               *
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
 * @file GFitsBinTable.cpp
 * @brief FITS binary table class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GFitsBinTable.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes of local functions ______________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFitsBinTable::GFitsBinTable(void) : GFitsTable()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Table constructor
 *
 * @param[in] nrows Number of rows in table
 ***************************************************************************/
GFitsBinTable::GFitsBinTable(const int& nrows) : GFitsTable(nrows)
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
 * @param[in] table Binary table.
 ***************************************************************************/
GFitsBinTable::GFitsBinTable(const GFitsBinTable& table) : GFitsTable(table)
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
GFitsBinTable::~GFitsBinTable(void)
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
 * @param[in] table Binary table.
 * @return Binary table.
 ***************************************************************************/
GFitsBinTable& GFitsBinTable::operator=(const GFitsBinTable& table)
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
 * @brief Clear binary table
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GFitsBinTable::clear(void)
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
 * @brief Clone binary table
 *
 * @return Pointer to deep copy of binary table.
 ***************************************************************************/
GFitsBinTable* GFitsBinTable::clone(void) const
{
    return new GFitsBinTable(*this);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsBinTable::init_members(void)
{
    // Initialise members
    m_type = 2;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] table Table to copy
 ***************************************************************************/
void GFitsBinTable::copy_members(const GFitsBinTable& table)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsBinTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise binary table header
 *
 * Initialises the table header by setting the default header cards.
 ***************************************************************************/
void GFitsBinTable::init_table_header(void)
{
    // Compute total width in Bytes
    int width = 0;
    for (int i = 0; i < ncols(); ++i) {
        width += m_columns[i]->width();
    }

    // Set image header keywords
    m_header.append(GFitsHeaderCard("XTENSION", "BINTABLE",
                                    "binary table extension"));
    m_header.append(GFitsHeaderCard("BITPIX", 8,
                                    "8-bit bytes"));
    m_header.append(GFitsHeaderCard("NAXIS", 2,
                                    "2-dimensional binary table"));
    m_header.append(GFitsHeaderCard("NAXIS1", width,
                                    "width of table in bytes"));
    m_header.append(GFitsHeaderCard("NAXIS2", nrows(),
                                    "number of rows in table"));
    m_header.append(GFitsHeaderCard("PCOUNT", 0,
                                    "size of special data area"));
    m_header.append(GFitsHeaderCard("GCOUNT", 1,
                                    "one data group (required keyword)"));

    // Return
    return;
}
