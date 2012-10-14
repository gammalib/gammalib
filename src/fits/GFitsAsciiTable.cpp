/***************************************************************************
 *             GFitsAsciiTable.cpp  - FITS ASCII table class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @brief Constructor
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable(void) : GFitsTable()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] nrows Number of rows in table.
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable(int nrows) : GFitsTable(nrows)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] table Table.
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
 * @param[in] table Table.
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
 * @brief Clear instance
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
