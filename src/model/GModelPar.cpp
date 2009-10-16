/***************************************************************************
 *                   GModelPar.cpp  -  Model parameter class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
 * @file GModelPar.cpp
 * @brief GModelPar class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GModelPar.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_VALUE                             "GModelPar::value(const double&)"
#define G_MIN                                 "GModelPar::min(const double&)"
#define G_MAX                                 "GModelPar::max(const double&)"
#define G_RANGE               "GModelPar::range(const double&,const double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                     GModelPar constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GModelPar::GModelPar(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Parameter from which the instance should be built.
 ***************************************************************************/
GModelPar::GModelPar(const GModelPar& par)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelPar::~GModelPar()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GModelPar operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] par Parameter which should be assigned.
 ***************************************************************************/
GModelPar& GModelPar::operator= (const GModelPar& par)
{ 
    // Execute only if object is not identical
    if (this != &par) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(par);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GModelPar public methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 ***************************************************************************/
void GModelPar::value(const double& value)
{
    // If there is a minimum boundary and if value is below this boundary
    // then throw an error
    if (m_hasmin && value < m_min)
        throw GException::out_of_range(G_VALUE, value, m_min, m_max);

    // If there is a maximum boundary and if value is above this boundary
    // then throw an error
    if (m_hasmax && value > m_max)
        throw GException::out_of_range(G_VALUE, value, m_min, m_max);

    // Assign value
    m_value = value;
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter minimum
 *
 * @param[in] min Parameter minimum.
 ***************************************************************************/
void GModelPar::min(const double& min)
{
    // If minimum is above actual value then throw error
    if (m_value < min)
        throw GException::out_of_range(G_MIN, m_value, min, m_max);

    // Assign minimum
    m_min = min;
    
    // Flag that minimum was set
    m_hasmin = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter maximum
 *
 * @param[in] max Parameter maximum.
 ***************************************************************************/
void GModelPar::max(const double& max)
{
    // If maximum is below value then throw error
    if (m_value > max)
        throw GException::out_of_range(G_MAX, m_value, m_min, max);

    // Assign maximum
    m_max = max;
    
    // Flag that maximum was set
    m_hasmax = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter minimum and maximum
 *
 * @param[in] min Parameter minimum.
 * @param[in] max Parameter maximum.
 ***************************************************************************/
void GModelPar::range(const double& min, const double& max)
{
    // If maximum is below value then throw error
    if (m_value < min || m_value > max || min > max)
        throw GException::out_of_range(G_RANGE, m_value, min, max);

    // Assign range
    m_min = min;
    m_max = max;
    
    // Flag that range was set
    m_hasmin = true;
    m_hasmax = true;
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        GModelPar private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelPar::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_unit.clear();
    m_value    = 0.0;
    m_error    = 0.0;
    m_gradient = 0.0;
    m_min      = 0.0;
    m_max      = 0.0;
    m_scale    = 1.0;
    m_hasmin   = false;
    m_hasmax   = false;
    m_free     = true;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par GModelPar members which should be copied.
 ***************************************************************************/
void GModelPar::copy_members(const GModelPar& par)
{
    // Copy members
    m_name    = par.m_name;
    m_unit    = par.m_unit;
    m_value   = par.m_value;
    m_error   = par.m_error;
    m_gradient = par.m_gradient;
    m_min     = par.m_min;
    m_max     = par.m_max;
    m_scale   = par.m_scale;
    m_hasmin  = par.m_hasmin;
    m_hasmax  = par.m_hasmax;
    m_free    = par.m_free;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelPar::free_members(void)
{
  
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GModelPar friends                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put parameter in output stream
 *
 * @param[in] os Output stream into which the parameter will be dumped
 * @param[in] par Parameter to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelPar& par)
{
    // Put parameter in stream
    os << std::scientific;
    os << par.m_name << ": ";
    os << par.real_value();
    os << " +/- ";
    os << par.real_error();
    if (par.m_hasmin && par.m_hasmax)
        os << " [" << par.real_min() << "," << par.real_max() << "]";
    else if (par.m_hasmin)
        os << " [" << par.real_min() << ",infty[";
    else if (par.m_hasmax)
        os << " ]-infty," << par.real_max() << "]";
    os << " " << par.m_unit;
    if (par.m_free)
        os << " (free,";
    else
        os << " (fixed,";
    os << "scale=" << par.m_scale << ")";

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GModelPar                   =
 =                                                                         =
 ==========================================================================*/
