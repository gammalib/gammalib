/***************************************************************************
 *          GModelSpatial.cpp - Abstract spatial model base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GModelSpatial.cpp
 * @brief Abstract spatial model base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelSpatial.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                    "GModelSpatial::operator[](std::string&)"
#define G_AT                             "GModelPar& GModelSpatial::at(int&)"

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
GModelSpatial::GModelSpatial(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spatial model.
 ***************************************************************************/
GModelSpatial::GModelSpatial(const GModelSpatial& model)
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
GModelSpatial::~GModelSpatial(void)
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
 * @param[in] model Spatial model.
 * @return Spatial model.
 ***************************************************************************/
GModelSpatial& GModelSpatial::operator= (const GModelSpatial& model)
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
 * @param[in] name Parameter name.
 * @return Model parameter reference.
 *
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
 *
 * Returns reference to the model parameter of specified @p name.
 ***************************************************************************/
GModelPar& GModelSpatial::operator[](const std::string& name)
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        throw GException::par_not_found(G_ACCESS, name);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] name Parameter name.
 * @return Model parameter reference.
 *
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
 *
 * Returns reference to the model parameter of specified @p name.
 ***************************************************************************/
const GModelPar& GModelSpatial::operator[](const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        throw GException::par_not_found(G_ACCESS, name);
    }

    // Return reference
    return *(m_pars[index]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Returns model parameter with @p index range checking.
 ***************************************************************************/
GModelPar& GModelSpatial::at(const int& index)
{
    // Compile option: raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Returns model parameter with @p index range checking.
 ***************************************************************************/
const GModelPar& GModelSpatial::at(const int& index) const
{
    // Compile option: raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Autoscale parameters
 *
 * Sets the scale factors for all parameters so that the values are unity.
 ***************************************************************************/
void GModelSpatial::autoscale(void)
{
    // Loop over all parameters
    for (int i = 0; i < m_pars.size(); ++i) {
        if (m_pars[i] != NULL) {
            m_pars[i]->autoscale();
        }
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatial::init_members(void)
{
    // Initialise members
    m_pars.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spatial model.
 ***************************************************************************/
void GModelSpatial::copy_members(const GModelSpatial& model)
{
    // Copy members
    m_pars = model.m_pars;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatial::free_members(void)
{
    // Return
    return;
}
