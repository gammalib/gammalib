/***************************************************************************
 *                       GLog.cpp - Information logger                     *
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
 * @file GLog.hpp
 * @brief Information logger class implementation
 * @author Jurgen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GLog.hpp"

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
GLog::GLog(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] log Object from which the instance should be built.
 ***************************************************************************/
GLog::GLog(const GLog& log)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(log);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLog::~GLog(void)
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
 * @param[in] log Object which should be assigned.
 ***************************************************************************/
GLog& GLog::operator= (const GLog& log)
{ 
    // Execute only if object is not identical
    if (this != &log) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(log);

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
void GLog::init_members(void)
{
    // Initialise members
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] log Object from which members which should be copied.
 ***************************************************************************/
void GLog::copy_members(const GLog& log)
{
    // Copy attributes
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLog::free_members(void)
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
 * @param[in] log Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLog& log)
{
    // Put object in stream
    os << "=== GLog ===" << std::endl;

    // Return output stream
    return os;
}


