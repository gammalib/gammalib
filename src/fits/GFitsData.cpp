/***************************************************************************
 *         GFitsData.cpp  - FITS data handling abstract base class         *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 * @brief Void constructor
 ***************************************************************************/
GFitsData::GFitsData(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***************************************************************************
 * @brief Copy constructor
 *
 * @param data Data to use for construction.
 ***************************************************************************/
GFitsData::GFitsData(const GFitsData& data)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(data);

    // Return
    return;
}


/***************************************************************************
 * @brief Destructor
 ***************************************************************************/
GFitsData::~GFitsData(void)
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
 * @param data Data to be assigned.
 ***************************************************************************/
GFitsData& GFitsData::operator= (const GFitsData& data)
{
    // Execute only if object is not identical
    if (this != &data) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(data);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                             Protected methods                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Connect data to FITS file
 *
 * @param[in] vptr Date file pointer.
 *
 * The data is connected to a FITS file by copying the FITS file pointer that
 * points towards the FITS data structure into the m_fitsfile member.
 ***************************************************************************/
void GFitsData::connect(void* vptr)
{
    // Connect data by copying the data file pointer
    FPTR_COPY(m_fitsfile, vptr);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsData::init_members(void)
{
    // Initialise members
    m_fitsfile = new __fitsfile;
    FPTR(m_fitsfile)->HDUposition = 0;
    FPTR(m_fitsfile)->Fptr        = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] data Data to be copied.
 ***************************************************************************/
void GFitsData::copy_members(const GFitsData& data)
{
    // Copy members
    FPTR_COPY(m_fitsfile, data.m_fitsfile);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsData::free_members(void)
{
    // Free memory
    if (m_fitsfile != NULL) delete FPTR(m_fitsfile);

    // Properly mark as free
    m_fitsfile = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
