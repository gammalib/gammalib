/***************************************************************************
 *        GModelTemporalConst.cpp  -  Temporal constant model class        *
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
 * @file GModelTemporalConst.cpp
 * @brief GModelTemporalConst class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelTemporalRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelTemporalConst    g_temporal_const_seed;
const GModelTemporalRegistry g_temporal_const_registry(&g_temporal_const_seed);

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelTemporalConst::GModelTemporalConst(void) : GModelTemporal()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model from which the instance should be built.
 ***************************************************************************/
GModelTemporalConst::GModelTemporalConst(const GModelTemporalConst& model) :
    GModelTemporal(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelTemporalConst::~GModelTemporalConst(void)
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
 * @param[in] model Model which should be assigned.
 ***************************************************************************/
GModelTemporalConst& GModelTemporalConst::operator= (const GModelTemporalConst& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelTemporal::operator=(model);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GModelTemporalConst::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelTemporal::free_members();

    // Initialise members
    this->GModelTemporal::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelTemporalConst* GModelTemporalConst::clone(void) const
{
    return new GModelTemporalConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcTime True photon arrival time (not used).
 *
 * This method implements the temporal component of a constant model. It
 * simply returns the normalization constant (which typically is set to 1).
 ***************************************************************************/
double GModelTemporalConst::eval(const GTime& srcTime)
{
    // Return
    return (norm());
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcTime True photon arrival time (not used).
 *
 * This method implements the temporal component of a constant model. It
 * simply returns the normalization constant (which typically is set to 1)
 * and a gradient of 0 (as it is not expected to fit this parameter).
 ***************************************************************************/
double GModelTemporalConst::eval_gradients(const GTime& srcTime)
{
    // Set gradient to 0
    m_norm.gradient(0.0);

    // Return
    return (norm());
}


/***********************************************************************//**
 * @brief Returns vector of random event times
 *
 * @param[in] rate Mean event rate (events per second).
 * @param[in] tmin Minimum event time.
 * @param[in] tmax Maximum event time.
 * @param[in] ran Random number generator.
 *
 * This method returns a vector of random event times assuming a constant
 * event rate that is specified by the rate parameter.
 ***************************************************************************/
GTimes GModelTemporalConst::mc(const double& rate, const GTime&  tmin,
                               const GTime&  tmax, GRan& ran)
{
    // Allocates empty vector of times
    GTimes times;

    // Compute event rate (in events per seconds)
    double lambda = rate * norm();

    // Initialise start and stop times (using MET)
    double time  = tmin.met();
    double tstop = tmax.met();

    // Generate events until maximum event time is exceeded
    while (time <= tstop) {

        // Simulate next event time
        time += ran.exp(lambda);

        // Add time if it is not beyod the stop time
        if (time <= tstop) {
            GTime event;
            event.met(time);
            times.push_back(event);
        }

    } // endwhile: loop until stop time is reached

    // Return vector of times
    return times;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
void GModelTemporalConst::read(const GXmlElement& xml)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
void GModelTemporalConst::write(GXmlElement& xml) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print constant information
 ***************************************************************************/
std::string GModelTemporalConst::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelTemporalConst ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_par[i]->print());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelTemporalConst::init_members(void)
{
    // Initialise parameters
    m_npars  = 1;
    m_par[0] = &m_norm;

    // Initialise normalisation parameter
    m_norm = GModelPar();
    m_norm.name("Constant");
    m_norm.unit("(relative value)");
    m_norm.value(1.0);
    m_norm.fix();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelTemporalConst members which should be copied.
 ***************************************************************************/
void GModelTemporalConst::copy_members(const GModelTemporalConst& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_norm  = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelTemporalConst::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
