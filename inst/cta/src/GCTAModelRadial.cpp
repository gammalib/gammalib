/***************************************************************************
 *         GCTAModelRadial.cpp  -  Abstract radial model base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GCTAModelRadial.cpp
 * @brief Abstract radial acceptance model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GCTAModelRadial.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                           "GModelSpatial::operator[](int&)"
#define G_ACCESS2                   "GModelSpatial::operator[](std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelRadial::GCTAModelRadial(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
GCTAModelRadial::GCTAModelRadial(const GCTAModelRadial& model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelRadial::~GCTAModelRadial(void)
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
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
GCTAModelRadial& GCTAModelRadial::operator= (const GCTAModelRadial& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar& GCTAModelRadial::operator[](const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS1, index, 0, size()-1);
    #endif

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
const GModelPar& GCTAModelRadial::operator[](const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS1, index, 0, size()-1);
    #endif

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
 ***************************************************************************/
GModelPar& GCTAModelRadial::operator[](const std::string& name)
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name)
            break;
    }

    // Throw exception if parameter name was not found
    if (index >= size())
        throw GException::par_not_found(G_ACCESS2, name);

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter (const version)
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
 ***************************************************************************/
const GModelPar& GCTAModelRadial::operator[](const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name)
            break;
    }

    // Throw exception if parameter name was not found
    if (index >= size())
        throw GException::par_not_found(G_ACCESS2, name);

    // Return reference
    return *(m_pars[index]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelRadial::init_members(void)
{
    // Initialise members
    m_pars.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
void GCTAModelRadial::copy_members(const GCTAModelRadial& model)
{
    // Copy members
    m_pars = model.m_pars;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRadial::free_members(void)
{
    // Return
    return;
}
