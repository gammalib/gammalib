/***************************************************************************
 *       GSkyPixel.cpp  -  Class that implements a 2D sky pixel index      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
#include "GSkyPixel.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototype __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                     GSkyPixel constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GSkyPixel::GSkyPixel()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Index constructor
 *
 * @param[in] x X index.
 * @param[in] y Y index.
 ***************************************************************************/
GSkyPixel::GSkyPixel(const double& x, const double& y)
{
    // Set members
    m_x = x;
    m_y = y;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pixel Sky pixel from which class should be instantiated.
 ***************************************************************************/
GSkyPixel::GSkyPixel(const GSkyPixel& pixel)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(pixel);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyPixel::~GSkyPixel()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GSkyPixel operators                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] pixel Sky pixel to be assigned.
 ***************************************************************************/
GSkyPixel& GSkyPixel::operator= (const GSkyPixel& pixel)
{
    // Execute only if object is not identical
    if (this != &pixel) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(pixel);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GSkyPixel public methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set x value of sky pixel
 *
 * @param[in] x X value to set.
 ***************************************************************************/
inline void GSkyPixel::x(const double& x)
{
    // Set x value
    m_x = x;
}


/***********************************************************************//**
 * @brief Set y value of sky pixel
 *
 * @param[in] y Y value to set.
 ***************************************************************************/
inline void GSkyPixel::y(const double& y)
{
    // Set y value
    m_y = y;
}


/***********************************************************************//**
 * @brief Return x value of sky pixel
 ***************************************************************************/
inline double GSkyPixel::x(void) const
{
    // Return x value
    return m_x;
}


/***********************************************************************//**
 * @brief Return x value of sky pixel
 ***************************************************************************/
inline double GSkyPixel::y(void) const
{
    // Return y value
    return m_y;
}


/*==========================================================================
 =                                                                         =
 =                         GSkyPixel private methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyPixel::init_members(void)
{
    // Initialise members
    m_x = 0.0;
    m_y = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pixel Sky pixel from which members should be copied.
 ***************************************************************************/
void GSkyPixel::copy_members(const GSkyPixel& pixel)
{
    // Copy attributes
    m_x = pixel.m_x;
    m_y = pixel.m_y;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyPixel::free_members(void)
{
    // Free memory

    // Signal free pointers

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GSkyPixel friends                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] pixel Sky pixel to put in output stream.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSkyPixel& pixel)
{
    // Put pixel in output stream
    os << "(" << pixel.x() << "," << pixel.y() << ")";
    
    // Return output stream
    return os;
}
