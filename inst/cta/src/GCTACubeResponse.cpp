/***************************************************************************
 *      GCTACubeResponse.cpp - CTA cube-style response function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTACubeResponse.cpp
 * @brief CTA cube-style response function class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <string>
#include "GTools.hpp"
#include "GCTACubeResponse.hpp"
#include "GPhoton.hpp"
#include "GEvent.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GObservation.hpp"
#include "GCTAInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF        "GCTACubeResponse::irf(GEvent&, GPhoton& GObservation&)"
#define G_NPRED_DIFFUSE    "GCTACubeResponse::npred(GPhoton&, GObservation&)"

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
 * Constructs void CTA response.
 ***************************************************************************/
GCTACubeResponse::GCTACubeResponse(void) : GResponse()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp CTA response.
 *
 * Constructs CTA cube-style response by making a deep copy of an existing
 * object.
 **************************************************************************/
GCTACubeResponse::GCTACubeResponse(const GCTACubeResponse& rsp) : GResponse(rsp)
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
 *
 * Destroys instance of CTA response object.
 ***************************************************************************/
GCTACubeResponse::~GCTACubeResponse(void)
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
 * @param[in] rsp CTA response.
 * @return CTA response.
 *
 * Assigns CTA response object to another CTA response object. The assignment
 * performs a deep copy of all information, hence the original object from
 * which the assignment has been performed can be destroyed after this
 * operation without any loss of information.
 ***************************************************************************/
GCTACubeResponse& GCTACubeResponse::operator=(const GCTACubeResponse& rsp)
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
 * Clears CTA response object by resetting all members to an initial state.
 * Any information that was present in the object before will be lost.
 ***************************************************************************/
void GCTACubeResponse::clear(void)
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
 * @return Pointer to deep copy of CTA response.
 *
 * Creates a clone (deep copy) of a CTA response object.
 ***************************************************************************/
GCTACubeResponse* GCTACubeResponse::clone(void) const
{
    return new GCTACubeResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 ***************************************************************************/
double GCTACubeResponse::irf(const GEvent&       event,
                             const GPhoton&      photon,
                             const GObservation& obs) const
{
    // Retrieve event instrument direction
    const GCTAInstDir& dir = retrieve_dir(G_IRF, event);

    // Get event attributes
    const GSkyDir& obsDir = dir.dir();
    const GEnergy& obsEng = event.energy();

    // Get photon attributes
    const GSkyDir& srcDir  = photon.dir();
    const GEnergy& srcEng  = photon.energy();
    const GTime&   srcTime = photon.time();

    // Determine angular separation between true and measured photon
    // direction in radians
    double delta = obsDir.dist(srcDir);

    // Get maximum angular separation for PSF (in radians) and add 10%
    // of margin
    double delta_max = 1.1 * psf().delta_max();

    // Initialise IRF value
    double irf = 0.0;

    // Compute only if we're sufficiently close to PSF
    if (delta <= delta_max) {

        // Get exposure
        irf = exposure()(srcDir, srcEng);

        // Multiply-in PSF
        if (irf > 0.0) {

            // Get PSF component
            irf *= psf()(srcDir, delta, srcEng);

            // Divide by ontime as the binned likelihood function is
            // later multiplying by ontime
            irf /= obs.ontime();

            // Apply deadtime correction
            irf *= obs.deadc(srcTime);

        } // endif: Aeff was non-zero

    } // endif: we were sufficiently close to PSF

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCTACubeResponse::irf:";
        std::cout << " NaN/Inf encountered";
        std::cout << " irf=" << irf;
        std::cout << std::endl;
    }
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return spatial integral of point spread function
 *
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method not implemented.
 ***************************************************************************/
double GCTACubeResponse::npred(const GPhoton&      photon,
                               const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_NPRED_DIFFUSE,
          "Npred computation not implemented for cube-style analysis.");

    // Return Npred
    return 0.0;
}


/***********************************************************************//**
 * @brief Load CTA response
 *
 * @param[in] rspname CTA response name.
 *
 * Loads the CTA response with specified name @p rspname. The method first
 * searchs for an appropriate response in the calibration database. If no
 * appropriate response is found, the method takes the database root path
 * and response name to build the full path to the response file, and tries
 * to load the response from these paths.
 ***************************************************************************/
/*
void GCTACubeResponse::load(const std::string& rspname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // First attempt reading the response using the GCaldb interface
    std::string expr      = "NAME("+rspname+")";
    std::string aeffname  = m_caldb.filename("","","EFF_AREA","","",expr);
    std::string psfname   = m_caldb.filename("","","RPSF","","",expr);
    std::string edispname = m_caldb.filename("","","EDISP","","",expr);
    std::string bgdname   = m_caldb.filename("","","BGD","","",expr);

    // If filenames are empty then build filenames from CALDB root path and
    // response name
    if (aeffname.length() < 1) {
        aeffname = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
    }
    if (psfname.length() < 1) {
        psfname = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
    }
    if (edispname.length() < 1) {
        edispname = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
    }
    if (bgdname.length() < 1) {
        bgdname = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
    }

    // Load effective area
    load_aeff(aeffname);

    // Load point spread function
    load_psf(psfname);

    // Load energy dispersion
    load_edisp(edispname);

    // Load background
    load_background(bgdname);

    // Remove theta cut
    GCTAAeffArf* arf = const_cast<GCTAAeffArf*>(dynamic_cast<const GCTAAeffArf*>(m_aeff));
    if (arf != NULL) {
        arf->remove_thetacut(*this);
    }

    // Store response name
    m_rspname = rspname;

    // Return
    return;
}
*/


/***********************************************************************//**
 * @brief Print CTA response information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing CTA response information.
 ***************************************************************************/
std::string GCTACubeResponse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTACubeResponse ===");

        // Append exposure cube information
        result.append("\n"+m_exposure.print(chatter));

        // Append point spread function information
        result.append("\n"+m_psf.print(chatter));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =              Model type dependent CTA response methods                  =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                    Low-level CTA response methods                       =
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
void GCTACubeResponse::init_members(void)
{
    // Initialise members
    m_exposure.clear();
    m_psf.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GCTACubeResponse::copy_members(const GCTACubeResponse& rsp)
{
    // Copy members
    m_exposure = rsp.m_exposure;
    m_psf      = rsp.m_psf;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTACubeResponse::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Retrieve CTA observation from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs is not a CTA observations.
 *
 * Dynamically casts generic observation into a CTA observation. If the
 * generic observation is not a CTA observation, an exception is thrown.
 ***************************************************************************/
const GCTAObservation& GCTACubeResponse::retrieve_obs(const std::string& origin,
                                                      const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);

    // If pointer is not valid then throw an exception
    if (cta == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n"
                          "Please specify a CTA observation when calling"
                          " this method.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return reference
    return *cta;
}


/***********************************************************************//**
 * @brief Retrieve CTA instrument direction from generic event
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] event Generic event.
 *
 * @exception GException::invalid_argument
 *            @p event does not contain a CTA instrument direction.
 *
 * Extract CTA Instrument Direction from an event.
 ***************************************************************************/
const GCTAInstDir& GCTACubeResponse::retrieve_dir(const std::string& origin,
                                                  const GEvent&      event) const
{
    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&(event.dir()));

    // If pointer is not valid then throw an exception
    if (dir == NULL) {
        std::string msg = "Specified event is not a CTA event.\n"
                          "Please specify a CTA event when calling this method.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return reference
    return *dir;
}
