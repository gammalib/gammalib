/***************************************************************************
 *               GObservations.cpp  -  Observation container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
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
 * @file GObservations.cpp
 * @brief GObservations class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GObservations.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS1                        "GObservations::operator() (int)"
#define G_OP_ACCESS2                  "GObservations::operator() (int) const"
#define G_OBSERVATION                 "GObservations::observation(int) const"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GObservations::GObservations(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param obs Instance which should be used for construction
 ***************************************************************************/
GObservations::GObservations(const GObservations& obs)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GObservations::~GObservations(void)
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
 * @param[in] obs Instance to be assigned
 ***************************************************************************/
GObservations& GObservations::operator= (const GObservations& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Observation access operator
 *
 * @param[in] index Index of observation (0,1,2,...)
 *
 * @exception GException::out_of_range
 *            Operation index is out of range.
 ***************************************************************************/
GObservation& GObservations::operator() (int index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num)
        throw GException::out_of_range(G_OP_ACCESS1, index, 0, m_num-1);

    // Return observation pointer
    return *(m_obs[index]);
}


/***********************************************************************//**
 * @brief Observation access operator
 *
 * @param[in] index Index of observation (0,1,2,...)
 *
 * @exception GException::out_of_range
 *            Operation index is out of range.
 ***************************************************************************/
const GObservation& GObservations::operator() (int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num)
        throw GException::out_of_range(G_OP_ACCESS2, index, 0, m_num-1);

    // Return observation pointer
    return *(m_obs[index]);
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Append observation to container
 *
 * @param[in] obs Observation to be appended 
 ***************************************************************************/
void GObservations::append(GObservation& obs)
{
    // Allocate new observation pointers
    GObservation** new_obs = new GObservation*[m_num+1];

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
 * @brief Optimize model parameters using optimizer
 *
 * @param[in] opt Optimiser to be used for 
 ***************************************************************************/
void GObservations::optimize(GOptimizer& opt)
{
    // Set optimizer function
    GObservations::optimizer fct(this);

    // Optimise model parameters
    m_models = opt(fct, m_models);

    // Store total number of predicted events
    m_npred = fct.npred();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get iterator for first observation and event
 *
 * Any empty observations (NULL pointer) or observations with 0 events are
 * skipped.
 ***************************************************************************/
GObservations::iterator GObservations::begin(void)
{
    // Allocate iterator object
    GObservations::iterator iter(this);

    // Get first valid observation
    if (iter.m_this != NULL) {
        while (iter.m_index < iter.m_this->m_num) {
            if (iter.m_this->m_obs[iter.m_index] != NULL) {
                iter.m_obs = iter.m_this->m_obs[iter.m_index];
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
GObservations::iterator GObservations::end(void)
{
    // Allocate iterator object
    GObservations::iterator iter(this);

    // Set obeservation number beyond last observation
    iter.m_index = iter.m_this->m_num;
    iter.m_obs   = NULL;

    // Return
    return iter;
}


/*==========================================================================
 =                                                                         =
 =                           GData private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * @todo Implement GModels::clear() method.
 ***************************************************************************/
void GObservations::init_members(void)
{
    // Initialise members
    m_num    = 0;
    m_obs    = NULL;
    m_models = GModels();
    m_npred  = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Instance from which members should be copied
 *
 * Copy observations from a GData object into the present object by invoking
 * the observation clone method of each observation.
 ***************************************************************************/
void GObservations::copy_members(const GObservations& obs)
{
    // Copy attributes
    m_num    = obs.m_num;
    m_models = obs.m_models;
    m_npred  = obs.m_npred;

    // Copy observations
    if (m_num > 0) {
        m_obs = new GObservation*[m_num];
        for (int i = 0; i < m_num; ++i) {
            if (obs.m_obs[i] != NULL)
                m_obs[i] = (obs.m_obs[i])->clone();
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
void GObservations::free_members(void)
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
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream into which the data will be dumped
 * @param[in] data Data to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GObservations& obs)
{
    // Put header in stream
    os << "=== GObservations ===" << std::endl;
    os << " Number of observations ....: " << obs.m_num << std::endl;

    // Put observations in stream
    for (int i = 0; i < obs.m_num; ++i)
        os << *(obs.m_obs[i]);

    // Add models to stream
    os << " Number of predicted events : " << obs.m_npred << std::endl;
    os << obs.m_models;

    // Return output stream
    return os;
}
