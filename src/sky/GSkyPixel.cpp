/***************************************************************************
 *       GSkyPixel.cpp  -  Class that implements a 2D sky pixel index      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GSkyPixel.hpp
 * @brief Sky pixel class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSkyPixel.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototype __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GSkyPixel::GSkyPixel(void)
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
GSkyPixel::GSkyPixel(double x, double y)
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
 * @param[in] pixel Sky pixel.
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
GSkyPixel::~GSkyPixel(void)
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
 * @param[in] pixel Sky pixel.
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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GSkyPixel::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GSkyPixel* GSkyPixel::clone(void) const
{
    // Clone this image
    return new GSkyPixel(*this);
}


/***********************************************************************//**
 * @brief Set x value of sky pixel
 *
 * @param[in] x X value.
 ***************************************************************************/
void GSkyPixel::x(const double& x)
{
    // Set x value
    m_x = x;
}


/***********************************************************************//**
 * @brief Set y value of sky pixel
 *
 * @param[in] y Y value.
 ***************************************************************************/
void GSkyPixel::y(const double& y)
{
    // Set y value
    m_y = y;
}


/***********************************************************************//**
 * @brief Return x value of sky pixel
 ***************************************************************************/
double GSkyPixel::x(void) const
{
    // Return x value
    return m_x;
}


/***********************************************************************//**
 * @brief Return x value of sky pixel
 ***************************************************************************/
double GSkyPixel::y(void) const
{
    // Return y value
    return m_y;
}


/***********************************************************************//**
 * @brief Print pixel
 ***************************************************************************/
std::string GSkyPixel::print(void) const
{
    // Initialise result string
    std::string result = "(";

    // Append pixel
    result.append(str(x()));
    result.append(",");
    result.append(str(y()));
    result.append(")");
    
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
 * @param[in] pixel Sky pixel.
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
    // Return
    return;
}
