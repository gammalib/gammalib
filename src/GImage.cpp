/***************************************************************************
 *                 GImage.cpp  -  Image abstract base class                *
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
/**
 * @file GImage.cpp
 * @brief GImage abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GImage.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                      GImage constructors/destructors                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GImage::GImage()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param image GImage instance which should be used for construction
 ***************************************************************************/
GImage::GImage(const GImage& image)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GImage::~GImage()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GImage operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] image GImage instance to be assigned
 ***************************************************************************/
GImage& GImage::operator= (const GImage& image)
{
    // Execute only if object is not identical
    if (this != &image) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(image);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GImage public methods                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                          GImage private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GImage::init_members(void)
{
    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] image GImage instance from which members should be copied
 ***************************************************************************/
void GImage::copy_members(const GImage& image)
{
    // Copy attributes
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GImage::free_members(void)
{
    // Free memory

    // Mark memory as free

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GImage friends                             =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GImage                     =
 =                                                                         =
 ==========================================================================*/
