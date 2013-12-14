/***************************************************************************
 *                     GSkyPixel.cpp - Sky map pixel class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @brief Sky map pixel class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSkyPixel.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INT                                     "GSkyPixel::operator int()"
#define G_DOUBLE                               "GSkyPixel::operator double()"

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
 * @brief 1D pixel constructor (integer version)
 *
 * @param[in] index Pixel index.
 ***************************************************************************/
GSkyPixel::GSkyPixel(const int& index)
{
    // Set members
    m_size = 1;
    m_x    = double(index);
    m_y    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief 1D pixel constructor (double precision version)
 *
 * @param[in] index Pixel index.
 ***************************************************************************/
GSkyPixel::GSkyPixel(const double& index)
{
    // Set members
    m_size = 1;
    m_x    = index;
    m_y    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief 2D pixel constructor (integer version)
 *
 * @param[in] x Pixel X index.
 * @param[in] y Pixel Y index.
 ***************************************************************************/
GSkyPixel::GSkyPixel(const int& x, const int& y)
{
    // Set members
    m_size = 2;
    m_x    = double(x);
    m_y    = double(y);

    // Return
    return;
}


/***********************************************************************//**
 * @brief 2D pixel constructor (double precision version)
 *
 * @param[in] x Pixel X index.
 * @param[in] y Pixel Y index.
 ***************************************************************************/
GSkyPixel::GSkyPixel(const double& x, const double& y)
{
    // Set members
    m_size = 2;
    m_x    = x;
    m_y    = y;

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
 * @param[in] pixel Sky map pixel.
 * @return Sky map pixel.
 ***************************************************************************/
GSkyPixel& GSkyPixel::operator=(const GSkyPixel& pixel)
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


/***********************************************************************//**
 * @brief To integer type conversion
 *
 * @return Pixel index.
 *
 * Converts the sky map pixel into an integer value.
 ***************************************************************************/
GSkyPixel::operator int() const
{
    // Throw an exception if pixel is not 1D
    if (!is1D()) {
        std::string msg = "Sky map pixel is not 1-dimensional.\n"
                          "Conversion from GSkyPixel to int is only allowed"
                          " for 1-dimensional sky map pixels.";
        throw GException::invalid_value(G_INT, msg);
    }

    // Round pixel to integer value
    int value = int(m_x + 0.5);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief To double type conversion
 *
 * @return Pixel index.
 *
 * Converts the sky map pixel into a double precision value.
 ***************************************************************************/
GSkyPixel::operator double() const
{
    // Throw an exception if pixel is not 1D
    if (!is1D()) {
        std::string msg = "Sky map pixel is not 1-dimensional.\n"
                          "Conversion from GSkyPixel to double is only allowed"
                          " for 1-dimensional sky map pixels.";
        throw GException::invalid_value(G_DOUBLE, msg);
    }

    // Return value
    return m_x;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * Set sky map pixel to a clean initial state.
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
 * @brief Clone sky map pixel
 *
 * @return Pointer to deep copy of sky map pixel.
 *
 * Returns a pointer to a deep copy of a sky map pixel.
 ***************************************************************************/
GSkyPixel* GSkyPixel::clone(void) const
{
    // Clone pixel
    return new GSkyPixel(*this);
}


/***********************************************************************//**
 * @brief Print pixel
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing pixel information.
 ***************************************************************************/
std::string GSkyPixel::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append pixel
        if (is1D()) {
            result.append(gammalib::str(x()));
        }
        else if (is2D()) {
            result.append("(");
            result.append(gammalib::str(x()));
            result.append(",");
            result.append(gammalib::str(y()));
            result.append(")");
        }

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
void GSkyPixel::init_members(void)
{
    // Initialise members
    m_size = 0;
    m_x    = 0.0;
    m_y    = 0.0;

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
    m_size = pixel.m_size;
    m_x    = pixel.m_x;
    m_y    = pixel.m_y;

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
