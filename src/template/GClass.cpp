/***************************************************************************
 *                       GClass.hpp - <brief descriptor>                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-20xx by <author>                                    *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GClass.cpp
 * @brief <brief descriptor>
 * @author <author>
 */

/* __ Includes ___________________________________________________________ */
#include "GClass.hpp"

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
GClass::GClass(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] c Object from which the instance should be built.
 ***************************************************************************/
GClass::GClass(const GClass& c)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(c);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GClass::~GClass(void)
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
 * @param[in] c Object which should be assigned.
 ***************************************************************************/
GClass& GClass::operator= (const GClass& c)
{ 
    // Execute only if object is not identical
    if (this != &c) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(c);

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
void GClass::init_members(void)
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
void GClass::copy_members(const GClass& c)
{
    // Copy attributes
    m_name  = c.m_name;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GClass::free_members(void)
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
 * @param[in] c Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GClass& c)
{
    // Put object in stream
    os << "=== GClass ===" << std::endl;

    // Return output stream
    return os;
}


