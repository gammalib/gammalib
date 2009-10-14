/***************************************************************************
 *                        GData.cpp  -  Data class                         *
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
 * @file GData.cpp
 * @brief GData class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GData.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ADD                                     "GData::add(GObservation&)"
#define G_COPY_MEMBERS                    "GData::copy_members(const GData&)"
#define G_OBSERVATION                         "GData::observation(int) const"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       GData constructors/destructors                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GData::GData()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param data GData instance which should be used for construction
 ***************************************************************************/
GData::GData(const GData& data)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(data);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GData::~GData()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GData operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] data GData instance to be assigned
 ***************************************************************************/
GData& GData::operator= (const GData& data)
{
    // Execute only if object is not identical
    if (this != &data) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(data);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                           GData event iterator                          =
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


/***********************************************************************//**
 * @brief Iterator Destructor
 ***************************************************************************/
GData::iterator::~iterator()
{
    // Return
    return;
}


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


/***********************************************************************//**
 * @brief Iterator == operator
 ***************************************************************************/
bool GData::iterator::operator==(const iterator& it) const
{
    // Return result
    return ((m_index == it.m_index) && (m_event == it.m_event));
}


/***********************************************************************//**
 * @brief Iterator != operator
 ***************************************************************************/
bool GData::iterator::operator!=(const iterator& it) const
{
    // Return result
    return ((m_index != it.m_index) || (m_event != it.m_event));
}


/***********************************************************************//**
 * @brief Iterator dereference operator
 ***************************************************************************/
GEvent& GData::iterator::operator*(void)
{
    // Return event
    return *m_event;
}


/***********************************************************************//**
 * @brief Iterator dereference operator
 ***************************************************************************/
GEvent* GData::iterator::operator->(void)
{
    // Return event
    return &(*m_event);
}


/***********************************************************************//**
 * @brief Get iterator for first observation and event
 *
 * Any empty observations (NULL pointer) or observations with 0 events are
 * skipped.
 ***************************************************************************/
GData::iterator GData::begin(void)
{
    // Allocate iterator object
    GData::iterator iter(this);
    
    // Get first valid observation
    if (iter.m_data != NULL) {
        while (iter.m_index < iter.m_data->m_num) {
            if (iter.m_data->m_obs[iter.m_index] != NULL) {
                iter.m_obs         = iter.m_data->m_obs[iter.m_index];
                break;
            }
        }
    }
    
    // Initialise event iterator
    if (iter.m_obs != NULL) {
        iter.m_event = iter.m_obs->events()->begin();
        iter.m_end   = iter.m_obs->events()->end();
    }
	
    // Return
    return iter;
}


/***********************************************************************//**
 * @brief Get iterator after last observation and event
 ***************************************************************************/
GData::iterator GData::end(void)
{
    // Allocate iterator object
    GData::iterator iter(this);
    
    // Set obeservation number beyond last observation
    iter.m_index = iter.m_data->m_num;
    iter.m_obs   = NULL;
	
    // Return
    return iter;
}


/*==========================================================================
 =                                                                         =
 =                           GData public methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Add observation to data
 ***************************************************************************/
void GData::add(GObservation& obs)
{
	// Allocate new observation pointers
    GObservation** new_obs = new GObservation*[m_num+1];
	if (new_obs == NULL)
		throw GException::mem_alloc(G_ADD, m_num+1);

	// If we have already observation pointers then copy them over to the
    // new pointer array
    if (m_num > 0) {
        for (int i = 0; i < m_num; ++i)
			new_obs[i] = m_obs[i];
	}
	
	// Create a copy of the observation that should be added and store the
    // pointer to this copy as last element of the pointer array
	new_obs[m_num] = obs.clone();
	
	// Release old pointer array
	if (m_obs != NULL) delete [] m_obs;
	
	// Attach new pointer array to this object
	m_obs = new_obs;
	
	// Increment number of observations
	m_num++;
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer on observation
 *
 * @param[in] index Index of observation (0,1,2,...)  
 ***************************************************************************/
GObservation* GData::observation(int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num)
        throw GException::out_of_range(G_OBSERVATION, index, 0, m_num-1);

    // Return observation pointer
    return m_obs[index];
}

/*==========================================================================
 =                                                                         =
 =                           GData private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GData::init_members(void)
{
    // Initialise members
    m_num    = 0;
	m_obs    = NULL;
    m_models = GModels();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] data GData instance from which members should be copied
 *
 * Copy observations from a GData object into the present object by invoking
 * the observation clone method of each observation.
 ***************************************************************************/
void GData::copy_members(const GData& data)
{
    // Copy attributes
    m_num    = data.m_num;
    m_models = data.m_models;
	
	// Copy observations
    if (m_num > 0) {
        m_obs = new GObservation*[m_num];
        if (m_obs == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_num);
        for (int i = 0; i < m_num; ++i) {
			if (data.m_obs[i] != NULL)
				m_obs[i] = (data.m_obs[i])->clone();
			else
				m_obs[i] = NULL;
		}
	}
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GData::free_members(void)
{
    // Free memory
    if (m_obs != NULL) {
        for (int i = 0; i < m_num; ++i) {
			if (m_obs[i] != NULL) delete m_obs[i];
        }
		delete [] m_obs;
	}
	
    // Mark memory as free
	m_obs = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               GData friends                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream into which the data will be dumped
 * @param[in] data Data to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GData& data)
{
    // Put header in stream
    os << "=== GData ===" << std::endl;
    os << " Number of observations ....: " << data.m_num << std::endl;

	// Put observations in stream
	for (int i = 0; i < data.m_num; ++i)
		os << *(data.m_obs[i]);

    // Add models to stream
    os << data.m_models;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GData                     =
 =                                                                         =
 ==========================================================================*/
