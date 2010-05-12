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
    m_name  = app.m_name;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplication::free_members(void)
{
    // Free memory
    //if (m_par      != NULL) delete [] m_par;

    // Signal free pointers
    //m_par      = NULL;
  
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


