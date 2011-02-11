/***************************************************************************
 *               GObservations.cpp  -  Observation container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservations.cpp
 * @brief Observations container class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GObservations.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS                         "GObservations::operator[](int&)"
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
 * @brief Void constructor
 ***************************************************************************/
GObservations::GObservations(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param obs Observation container.
 ***************************************************************************/
GObservations::GObservations(const GObservations& obs)
{
    // Initialise members
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
 * @param[in] obs Observation container.
 ***************************************************************************/
GObservations& GObservations::operator=(const GObservations& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to observation
 *
 * @param[in] index Index of observation [0,...,size()-1]
 *
 * @exception GException::out_of_range
 *            Operation index is out of range.
 ***************************************************************************/
GObservation& GObservations::operator[](const int& index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_OP_ACCESS, index, 0, size()-1);

    // Return reference
    return *(m_obs[index]);
}


/***********************************************************************//**
 * @brief Return reference to observation (const version)
 *
 * @param[in] index Index of observation [0,...,size()-1]
 *
 * @exception GException::out_of_range
 *            Operation index is out of range.
 ***************************************************************************/
const GObservation& GObservations::operator[](const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_OP_ACCESS, index, 0, size()-1);

    // Return reference
    return *(m_obs[index]);
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear container
 ***************************************************************************/
void GObservations::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append observation to container
 *
 * @param[in] obs Observation.
 *
 * This method appends an observation to the container by cloning it.
 ***************************************************************************/
void GObservations::append(GObservation& obs)
{
    // Clone observation and append to list
    m_obs.push_back(obs.clone());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load models from XML file
 *
 * @param[in] filename XML filename. 
 ***************************************************************************/
void GObservations::models(const std::string& filename)
{
    // Load models
    m_models.load(filename);

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
 * @brief Print observations information
 ***************************************************************************/
std::string GObservations::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GObservations ===\n");
    result.append(parformat("Number of observations")+str(size())+"\n");
    result.append(parformat("Number of predicted events")+str(npred()));

    // Append observations
    for (int i = 0; i < size(); ++i) {
        result.append("\n");
        result.append((*this)[i].print());
    }

    // Append models
    result.append("\n"+m_models.print());

    // Return result
    return result;
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
        while (iter.m_index < iter.m_this->size()) {
            if (iter.m_this->m_obs[iter.m_index] != NULL) {
                iter.m_obs = iter.m_this->m_obs[iter.m_index];
                break;
            }
        }
    }

    // Initialise event iterator (circumvent const correctness)
    if (iter.m_obs != NULL) {
        GEvents* events = (GEvents*)(iter.m_obs->events());
        iter.m_event    = events->begin();
        iter.m_end      = events->end();
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
    iter.m_index = iter.m_this->size();
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
 ***************************************************************************/
void GObservations::init_members(void)
{
    // Initialise members
    m_obs.clear();
    m_models.clear();
    m_npred = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation container.
 ***************************************************************************/
void GObservations::copy_members(const GObservations& obs)
{
    // Copy attributes
    m_obs    = obs.m_obs;
    m_models = obs.m_models;
    m_npred  = obs.m_npred;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservations::free_members(void)
{
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
 * @param[in] os Output stream.
 * @param[in] obs Observation container.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GObservations& obs)
{
     // Write observations in output stream
    os << obs.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] obs Observation container.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GObservations& obs)
{
    // Write observations into logger
    log << obs.print();

    // Return logger
    return log;
}
