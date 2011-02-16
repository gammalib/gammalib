/***************************************************************************
 *              GEvents.cpp  -  Abstract Event container class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEvents.cpp
 * @brief Abstract event container class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GEvents.hpp"

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
GEvents::GEvents(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] events Event container.
 ***************************************************************************/
GEvents::GEvents(const GEvents& events)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(events);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEvents::~GEvents(void)
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
 * @param[in] events Event container.
 ***************************************************************************/
GEvents& GEvents::operator=(const GEvents& events)
{
    // Execute only if object is not identical
    if (this != &events) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(events);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

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


/***********************************************************************//**
 * @brief Set energy boundaries
 *
 * @param[in] ebounds Energy boundaries.
 ***************************************************************************/
void GEvents::ebounds(const GEbounds& ebounds)
{
    // Store energy boundaries
    m_ebounds = ebounds;

    // Call (optional) energy boundary update method
    set_energies();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set Good Time Intervals
 *
 * @param[in] gti Good Time Intervals.
 ***************************************************************************/
void GEvents::gti(const GGti& gti)
{
    // Store Good Time Intervals
    m_gti = gti;

    // Call (optional) good time interval update method
    set_times();

    // Return
    return;
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
    m_num   = (events != NULL) ? events->size() : 0;
    m_index = 0;
    m_base  = events;

    // Return
    return;
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


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEvents::init_members(void)
{
    // Initialise members
    m_ebounds.clear();
    m_gti.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] events Event container.
 ***************************************************************************/
void GEvents::copy_members(const GEvents& events)
{
    // Copy members
    m_ebounds = events.m_ebounds;
    m_gti     = events.m_gti;

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
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] events Event container.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEvents& events)
{
     // Write events in output stream
    os << events.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] events Event container.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GEvents& events)
{
    // Write events into logger
    log << events.print();

    // Return logger
    return log;
}
