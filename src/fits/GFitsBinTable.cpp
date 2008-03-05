/***************************************************************************
 *              GFitsBinTable.cpp  - FITS binary table class               *
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
#include "GFitsBinTable.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes of local functions ______________________________________ */


/*==========================================================================
 =                                                                         =
 =                   GFitsBinTable constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsBinTable::GFitsBinTable() : GFitsTable()
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
 ***************************************************************************/
GFitsBinTable::GFitsBinTable(int nrows) : GFitsTable(nrows)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] table Table which will be used to construct GFitsBinTable
 *                  instance
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
GFitsBinTable::~GFitsBinTable()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsBinTable operators                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] table Table which will be assigned
 ***************************************************************************/
GFitsBinTable& GFitsBinTable::operator= (const GFitsBinTable& table)
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
 =                       GFitsBinTable public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clone binary table
 ***************************************************************************/
GFitsBinTable* GFitsBinTable::clone(void) const
{
    return new GFitsBinTable(*this);
}


/*==========================================================================
 =                                                                         =
 =                       GFitsBinTable private methods                     =
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
    // Copy members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsBinTable::free_members(void)
{
    // Free memory

    // Mark memory as freed

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GFitsBinTable friends                         =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GFitsBinTable                 =
 =                                                                         =
 ==========================================================================*/
