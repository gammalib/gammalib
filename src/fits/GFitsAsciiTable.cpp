/***************************************************************************
 *             GFitsAsciiTable.cpp  - FITS ASCII table class               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GFitsAsciiTable.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                  GFitsAsciiTable constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable() : GFitsTable()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] nrows Number of rows in table
 * @param[in] ncols Number of columns in table
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
 * @param table Table to be used to initialise instance
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
 =                        GFitsAsciiTable operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param table Table to be assigned
 ***************************************************************************/
GFitsAsciiTable& GFitsAsciiTable::operator= (const GFitsAsciiTable& table)
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
 =                      GFitsAsciiTable public methods                     =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                      GFitsAsciiTable private methods                    =
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
 * @param table Table to be copied
 ***************************************************************************/
void GFitsAsciiTable::copy_members(const GFitsAsciiTable& table)
{
    // Copy members
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsAsciiTable::free_members(void)
{
    // Free memory
    
    // Mark memory as freed
    
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
 =                          GFitsAsciiTable friends                        =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsAsciiTable                =
 =                                                                         =
 ==========================================================================*/
