/***************************************************************************
 *                GEvents.cpp  -  Events container class                   *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
 * @file GEvents.cpp
 * @brief GEvents container class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GEvents.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                      GEvents constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEvents::GEvents()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] events Events from which the instance should be built.
 ***************************************************************************/
GEvents::GEvents(const GEvents& events)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(events);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEvents::~GEvents()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GEvents operators                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] events Events which should be assigned.
 ***************************************************************************/
GEvents& GEvents::operator= (const GEvents& events)
{
    // Execute only if object is not identical
    if (this != &events) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(events);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GEvents event iterator                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Iterator Constructor
 ***************************************************************************/
GEvents::iterator::iterator()
{
    // Initialise iterator
    m_num   = 0;
    m_index = 0;
    m_base  = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Iterator Constructor
 ***************************************************************************/
GEvents::iterator::iterator(GEvents *events)
{
    // Initialise iterator
    m_num   = (events != NULL) ? events->elements() : 0;
    m_index = 0;
    m_base  = events;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Iterator Destructor
 ***************************************************************************/
GEvents::iterator::~iterator()
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Iterator prefix operator
 ***************************************************************************/
GEvents::iterator& GEvents::iterator::operator++(void)
{
    // Get next event
    m_index++;
    
    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Iterator postfix operator
 ***************************************************************************/
GEvents::iterator GEvents::iterator::operator++(int junk)
{
    // Save actual iterator
    GEvents::iterator actual = *this;
    
    // Increment using prefix iterator
    ++(*this);
    
    // Return actual iterator
    return actual;
}


/***********************************************************************//**
 * @brief Iterator == operator
 ***************************************************************************/
bool GEvents::iterator::operator==(const iterator& it) const
{
    // Return result
    return (m_index == it.m_index);
}


/***********************************************************************//**
 * @brief Iterator != operator
 ***************************************************************************/
bool GEvents::iterator::operator!=(const iterator& it) const
{
    // Return result
    return (m_index != it.m_index);
}


/***********************************************************************//**
 * @brief Iterator dereference operator
 ***************************************************************************/
GEvent& GEvents::iterator::operator*(void)
{
    // Return event
    return *(m_base->pointer(m_index));
}


/***********************************************************************//**
 * @brief Iterator dereference operator
 ***************************************************************************/
GEvent* GEvents::iterator::operator->(void)
{
    // Return event
    return m_base->pointer(m_index);
}


/***********************************************************************//**
 * @brief Get iterator for first event
 ***************************************************************************/
GEvents::iterator GEvents::begin(void)
{
    // Allocate iterator object
    GEvents::iterator iter(this);

    // Set iterator for first event
    iter.m_index = 0;
	
    // Return
    return iter;
}


/***********************************************************************//**
 * @brief Get iterator after last event
 ***************************************************************************/
GEvents::iterator GEvents::end(void)
{
    // Allocate iterator object
    GEvents::iterator iter(this);
    
    // Set iterator beyond last event
    iter.m_index = iter.m_num;
	
    // Return
    return iter;
}


/*==========================================================================
 =                                                                         =
 =                         GEvents public methods                          =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                         GEvents private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEvents::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] events GEvents members which should be copied.
 ***************************************************************************/
void GEvents::copy_members(const GEvents& events)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEvents::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GEvents friends                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream into which the events will be dumped
 * @param[in] events Events to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEvents& events)
{
    // Put header in stream
    os << "=== GEvents ===" << std::endl;
    os << " Number of elements ........: " << events.elements();
        
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GEvents                     =
 =                                                                         =
 ==========================================================================*/
