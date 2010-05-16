/***************************************************************************
 *                    GPar.cpp - Application parameter                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GPar.cpp
 * @brief Application parameter class implementation
 * @author Jurgen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GPar.hpp"
//#include "GTools.hpp"
//#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GPar::GPar(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 ***************************************************************************/
GPar::GPar(const std::string& name, const std::string& type,
           const std::string& mode, const std::string& value,
           const std::string& min, const std::string& max, 
           const std::string& desc)
{
    // Initialise private members for clean destruction
    init_members();
    
    // Set parameter attributes
    m_name  = name;
    m_type  = type;
    m_mode  = mode;
    m_value = value;
    m_min   = min;
    m_max   = max;
    m_desc  = desc;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Object from which the instance should be built.
 ***************************************************************************/
GPar::GPar(const GPar& par)
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
GPar::~GPar(void)
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
 * @param[in] par Object which should be assigned.
 ***************************************************************************/
GPar& GPar::operator= (const GPar& par)
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
void GPar::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_type.clear();
    m_mode.clear();
    m_value.clear();
    m_min.clear();
    m_max.clear();
    m_desc.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par Object from which members which should be copied.
 ***************************************************************************/
void GPar::copy_members(const GPar& par)
{
    // Copy attributes
    m_name  = par.m_name;
    m_type  = par.m_type;
    m_mode  = par.m_mode;
    m_value = par.m_value;
    m_min   = par.m_min;
    m_max   = par.m_max;
    m_desc  = par.m_desc;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPar::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
***************************************************************************/
GPar* GPar::clone(void) const
{
    return new GPar(*this);
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] GPar Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GPar& par)
{
    // Put object in stream
    os << par.m_name << "=" << par.m_value;
    if (par.m_min.length() > 0 && par.m_max.length() > 0)
        os << " (" << par.m_min << "-" << par.m_max << ")";
    else if (par.m_min.length() > 0)
        os << " (>" << par.m_min << ")";
    else if (par.m_max.length() > 0)
        os << " (<" << par.m_max << ")";
    os << " t=" << par.m_type << " m=" << par.m_mode;
    os << " " << par.m_desc;

    // Return output stream
    return os;
}


