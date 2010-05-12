/***************************************************************************
 *             GApplication.cpp - GammaLib application base class          *
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
 * @file GApplication.cpp
 * @brief GammaLib application base class
 * @author Jurgen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GApplication.hpp"

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
GApplication::GApplication(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Save the execution start time
    time(&m_tstart);
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Application constructor
 ***************************************************************************/
GApplication::GApplication(std::string name)
{
    // Initialise private members for clean destruction
    init_members();
    
    // Set application name
    m_name = name;

    // Save the execution start time
    time(&m_tstart);
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Object from which the instance should be built.
 ***************************************************************************/
GApplication::GApplication(const GApplication& app)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(app);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GApplication::~GApplication(void)
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
 * @param[in] app Object which should be assigned.
 ***************************************************************************/
GApplication& GApplication::operator= (const GApplication& app)
{ 
    // Execute only if object is not identical
    if (this != &app) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(app);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return application name
 ***************************************************************************/
std::string GApplication::name(void) const
{
    // Return name
    return m_name;
}


/***********************************************************************//**
 * @brief Return application elapsed time in CPU seconds
 ***************************************************************************/
double GApplication::telapse(void) const
{
    // Get actual time
    time_t acttime;
    time(&acttime);
    
    // Compute elapsed time
    double telapse = difftime(acttime, m_tstart);

    // Return elapsed time
    return telapse;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GApplication::init_members(void)
{
    // Initialise members
    m_name.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] c Object from which members which should be copied.
 ***************************************************************************/
void GApplication::copy_members(const GApplication& app)
{
    // Copy attributes
    m_name   = app.m_name;
    m_tstart = app.m_tstart;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplication::free_members(void)
{
    // Return
    return;
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
 * @param[in] app Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GApplication& app)
{
    // Put object in stream
    os << "=== GApplication ===" << std::endl;

    // Return output stream
    return os;
}


