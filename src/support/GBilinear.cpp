/***************************************************************************
 *                GBilinear.cpp - Bilinear interpolator class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Juergen Knoedlseder                              *
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
 * @file GBilinear.cpp
 * @brief Bilinear interpolator class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GBilinear.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GBilinear::GBilinear(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] interpolator Bilinear interpolator.
 ***************************************************************************/
GBilinear::GBilinear(const GBilinear& interpolator)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(interpolator);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GBilinear::~GBilinear(void)
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
 * @param[in] interpolator Bilinear interpolator.
 * @return Bilinear interpolator.
 ***************************************************************************/
GBilinear& GBilinear::operator=(const GBilinear& interpolator)
{
    // Execute only if object is not identical
    if (this != &interpolator) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(interpolator);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Interpolator
 *
 * @param[in] array Bilinear interpolator.
 * @return Interpolated value.
 ***************************************************************************/
double GBilinear::operator()(const double* array)
{
    // Perform interpolation
    double value = m_wgt1 * array[m_inx1] +
                   m_wgt2 * array[m_inx2] +
                   m_wgt3 * array[m_inx3] +
                   m_wgt4 * array[m_inx4];

    // Return interpolated value
    return value;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear node array
 ***************************************************************************/
void GBilinear::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone bilinear interolator
 *
 * @return Pointer to deep copy of bilinear interolator.
 ***************************************************************************/
GBilinear* GBilinear::clone(void) const
{
    return new GBilinear(*this);
}


/***********************************************************************//**
 * @brief Print nodes
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing nodes information.
 ***************************************************************************/
std::string GBilinear::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GBilinear ===");

        // Append indices
        result.append("\n"+gammalib::parformat("Indices"));
        result.append(gammalib::str(m_inx1)+" ");
        result.append(gammalib::str(m_inx2)+" ");
        result.append(gammalib::str(m_inx3)+" ");
        result.append(gammalib::str(m_inx4));

        // Append weights
        result.append("\n"+gammalib::parformat("Weights"));
        result.append(gammalib::str(m_wgt1)+" ");
        result.append(gammalib::str(m_wgt2)+" ");
        result.append(gammalib::str(m_wgt3)+" ");
        result.append(gammalib::str(m_wgt4));

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
void GBilinear::init_members(void)
{
    // Initialise members
    m_inx1 = 0;
    m_inx2 = 0;
    m_inx3 = 0;
    m_inx4 = 0;
    m_wgt1 = 0.0;
    m_wgt2 = 0.0;
    m_wgt3 = 0.0;
    m_wgt4 = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] interpolator Bilinear interpolator.
 ***************************************************************************/
void GBilinear::copy_members(const GBilinear& interpolator)
{
    // Copy members
    m_inx1 = interpolator.m_inx1;
    m_inx2 = interpolator.m_inx2;
    m_inx3 = interpolator.m_inx3;
    m_inx4 = interpolator.m_inx4;
    m_wgt1 = interpolator.m_wgt1;
    m_wgt2 = interpolator.m_wgt2;
    m_wgt3 = interpolator.m_wgt3;
    m_wgt4 = interpolator.m_wgt4;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GBilinear::free_members(void)
{
    // Return
    return;
}
