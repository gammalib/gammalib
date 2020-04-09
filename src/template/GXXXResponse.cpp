/***************************************************************************
 *                GXXXResponse.cpp - [INSTRUMENT] response class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GXXXResponse.hpp
 * @brief [INSTRUMENT] instrument response function class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <typeinfo>
#include "GException.hpp"
#include "GSource.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GXXXResponse.hpp"
#include "GXXXObservation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF           "GXXXResponse::irf(GEvent&, GSource&, GObservation&)"
#define G_NROI            "GXXXResponse::nroi(GModelSky&, GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_EBOUNDS                           "GXXXResponse::ebounds(GEnergy&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty [INSTRUMENT] response.
 ***************************************************************************/
GXXXResponse::GXXXResponse(void) : GResponse()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp [INSTRUMENT] response.
 **************************************************************************/
GXXXResponse::GXXXResponse(const GXXXResponse& rsp) : GResponse(rsp)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXXXResponse::~GXXXResponse(void)
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
 * @param[in] rsp [INSTRUMENT] response.
 * @return [INSTRUMENT] response.
 *
 * Assign [INSTRUMENT] response to this object. The assignment performs
 * a deep copy of all information, hence the original object from which the
 * assignment was made can be destroyed after this operation without any loss
 * of information.
 ***************************************************************************/
GXXXResponse& GXXXResponse::operator=(const GXXXResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(rsp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * Clears [INSTRUMENT] response by resetting all class members to an initial
 * state. Any information that was present before will be lost.
 ***************************************************************************/
void GXXXResponse::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GResponse::free_members();

    // Initialise members
    this->GResponse::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of [INSTRUMENT] response.
 ***************************************************************************/
GXXXResponse* GXXXResponse::clone(void) const
{
    return new GXXXResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of [INSTRUMENT] instrument response for a photon
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 * @return Instrument response $\f(cm^2 sr^{-1})$\f
 *
 * @exception GException::invalid_argument
 *            Observation is not a [INSTRUMENT] observation.
 *
 * Returns the instrument response function for a given observed photon
 * direction as function of the assumed true photon direction. The result
 * is given by
 *
 * \f[
 *    R(p'|p) = 
 * \f]
 *
 * @todo Write down formula
 * @todo Describe in detail how the response is computed.
 * @todo Implement method.
 ***************************************************************************/
double GXXXResponse::irf(const GEvent&       event,
                         const GPhoton&      photon,
                         const GObservation& obs) const
{
    // Extract [INSTRUMENT] observation
    const GXXXObservation* observation = dynamic_cast<const GXXXObservation*>(&obs);
    if (observation == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Observation of type \""+cls+"\" is not a [INSTRUMENT] "
                          "observation. Please specify a [INSTRUMENT] observation "
                          "as argument.";
        throw GException::invalid_argument(G_IRF, msg);
    }

    // TODO: Now comes the real magic, this is your job !!!

    // Initialise IRF
    double irf = 0.0;

    // Compile option: Check for NaN/Inf
    // TODO: Add any relevant information that may help you to debug the
    // response code. You want to make sure that any NaN is signalled.
    // This can be switched off upon configuration.
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GXXXResponse::irf:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ")";
        std::cout << std::endl;
    }
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return integral of event probability for a given sky model over ROI
 *
 * @param[in] model Sky model.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] obs Observation.
 * @return 0.0
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GXXXResponse::nroi(const GModelSky&    model,
                          const GEnergy&      obsEng,
                          const GTime&        obsTime,
                          const GObservation& obs) const
{
    // Method is not implemented
    std::string msg = "Spatial integration of sky model over the data space "
                      "is not implemented.";
    throw GException::feature_not_implemented(G_NROI, msg);

    // Return Npred
    return (0.0);
}


/***********************************************************************//**
 * @brief Return true energy boundaries for a specific observed energy
 *
 * @param[in] obsEnergy Observed Energy.
 * @return True energy boundaries for given observed energy.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 *
 * @todo Implement this method if you need energy dispersion.
 ***************************************************************************/
GEbounds GXXXResponse::ebounds(const GEnergy& obsEnergy) const
{
    // Initialise an empty boundary object
    GEbounds ebounds;

    // Throw an exception
    std::string msg = "Energy dispersion not implemented.";
    throw GException::feature_not_implemented(G_EBOUNDS, msg);

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Print [INSTRUMENT] response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing [INSTRUMENT] response information.
 ***************************************************************************/
std::string GXXXResponse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GXXXResponse ===");

        // Append information
        // TODO: Add any relevant information

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GXXXResponse::init_members(void)
{
    // Initialise members
    // TODO: Initialise all data members
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp [INSTRUMENT] response function.
 ***************************************************************************/
void GXXXResponse::copy_members(const GXXXResponse& rsp)
{
    // Copy members
    // TODO: Copy all data members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXResponse::free_members(void)
{
    // Return
    return;
}
