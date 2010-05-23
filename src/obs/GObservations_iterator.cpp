/***************************************************************************
 *   GObservations_iterator.cpp  -  Iterator class of observations class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservations_iterator.cpp
 * @brief GObservations::iterator class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservations.hpp"

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
 * @brief Void iterator constructor
 ***************************************************************************/
GObservations::iterator::iterator() 
{
    // Initialise iterator
    m_data  = NULL;
    m_index = 0;
    m_obs   = NULL;
    m_event = GEvents::iterator();
    m_end   = GEvents::iterator();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Iterator constructor based on specific object
 *
 * @param[in] obs Observations container to use for iterator.
 ***************************************************************************/
GObservations::iterator::iterator(GObservations *obs)
{
    // Initialise iterator
    m_data  = obs;
    m_index = 0;
    m_obs   = NULL;
    m_event = GEvents::iterator();
    m_end   = GEvents::iterator();
    
    // Return
    return;
}



/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Iterator prefix operator
 *
 * The iterator skips any empty observations (i.e. NULL observation pointer)
 * and all observations with 0 events.
 ***************************************************************************/
GObservations::iterator& GObservations::iterator::operator++(void)
{
    // Only iterate if we have a valid data pointer
    if (m_data != NULL) {
        
        // Increment event iterator
        m_event++;

        // If end of observation is reached then set iterator to first event of
        // next valid observation
        if (m_event == m_end) {

            // Find next valid observation
            while (m_index < m_data->m_num) {
        
                // Go to next observation
                m_index++;
            
                // If we still have an observation then set now the iterator
                // to the first event from this observation ...
                if (m_index < m_data->m_num) {
            
                    // Get next observation. Skip if empty
                    m_obs = m_data->m_obs[m_index];
                    if (m_obs == NULL)
                        continue;
                
                    // Skip observation if it contains no events 
                    if (m_obs->events() == NULL)
                        continue;
                    
                    // Set iterator to first event of current observation. Exit
                    // only if we have events in that observation
                    m_event = m_obs->events()->begin();
                    m_end   = m_obs->events()->end();
                    if (m_event != m_end)
                        break;
                
                } // endif: we still had an observation
            
            } // endwhile: searched for next valid observation
        
            // If observations are exhausted then set iterator to signal the end
            if (m_index >= m_data->m_num) {
                m_index = m_data->m_num;
                m_obs   = NULL;
                m_event = GEvents::iterator();
                m_end   = GEvents::iterator();
            }
        
        } // endif: end of observation has been reached
    } // endif: we had a valid data pointer
    
    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Iterator postfix operator
 ***************************************************************************/
GObservations::iterator GObservations::iterator::operator++(int junk)
{
    // Save actual iterator
    GObservations::iterator actual = *this;
    
    // Increment using prefix iterator
    ++(*this);
    
    // Return actual iterator
    return actual;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
