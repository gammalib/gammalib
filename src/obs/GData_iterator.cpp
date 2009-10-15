/***************************************************************************
 *             GData_iterator.cpp  -  Iterator class of data class         *
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
 * @file GData_iterator.cpp
 * @brief GData::iterator class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GData.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GData::iterator constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void iterator constructor
 ***************************************************************************/
GData::iterator::iterator() 
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
 * @brief Iterator constructor based on specific GData object
 ***************************************************************************/
GData::iterator::iterator(GData *data)
{
    // Initialise iterator
    m_data  = data;
    m_index = 0;
    m_obs   = NULL;
    m_event = GEvents::iterator();
    m_end   = GEvents::iterator();
    
    // Return
    return;
}



/*==========================================================================
 =                                                                         =
 =                        GData::iterator operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Iterator prefix operator
 *
 * The iterator skips any empty observations (i.e. NULL observation pointer)
 * and all observations with 0 events.
 ***************************************************************************/
GData::iterator& GData::iterator::operator++(void)
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
GData::iterator GData::iterator::operator++(int junk)
{
    // Save actual iterator
    GData::iterator actual = *this;
    
    // Increment using prefix iterator
    ++(*this);
    
    // Return actual iterator
    return actual;
}


/*==========================================================================
 =                                                                         =
 =                      GData::iterator public methods                     =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                      GData::iterator private methods                    =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                          GData::iterator friends                        =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                  Other functions used by GData::iterator                =
 =                                                                         =
 ==========================================================================*/
