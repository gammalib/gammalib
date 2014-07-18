/***************************************************************************
 *      GCTAResponseCube.cpp - CTA cube-style response function class      *
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
 * @file GCTAResponseCube.cpp
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
#include "GCTAResponseCube.hpp"
#include "GPhoton.hpp"
#include "GEvent.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GObservation.hpp"
#include "GCTAInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF        "GCTAResponseCube::irf(GEvent&, GPhoton& GObservation&)"
#define G_NPRED_DIFFUSE    "GCTAResponseCube::npred(GPhoton&, GObservation&)"
#define G_READ                         "GCTAResponseCube::read(GXmlElement&)"

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
GCTAResponseCube::GCTAResponseCube(void) : GCTAResponse()
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
GCTAResponseCube::GCTAResponseCube(const GCTAResponseCube& rsp) :
                  GCTAResponse(rsp)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * @param[in] exposure CTA exposure.
 * @param[in] psf CTA mean point spread function.
 *
 * Constructs CTA cube-style response from a CTA exposure and a mean point
 * spread function.
 **************************************************************************/
GCTAResponseCube::GCTAResponseCube(const GCTAExposure& exposure,
                                   const GCTAMeanPsf&  psf) :
                  GCTAResponse()
{
    // Initialise members
    init_members();

    // Set members
    m_exposure = exposure;
    m_psf      = psf;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of CTA response object.
 ***************************************************************************/
GCTAResponseCube::~GCTAResponseCube(void)
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
GCTAResponseCube& GCTAResponseCube::operator=(const GCTAResponseCube& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GCTAResponse::operator=(rsp);

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
void GCTAResponseCube::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAResponse::free_members();
    this->GResponse::free_members();

    // Initialise members
    this->GResponse::init_members();
    this->GCTAResponse::init_members();
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
GCTAResponseCube* GCTAResponseCube::clone(void) const
{
    return new GCTAResponseCube(*this);
}


/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 ***************************************************************************/
double GCTAResponseCube::irf(const GEvent&       event,
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
        std::cout << "*** ERROR: GCTAResponseCube::irf:";
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
double GCTAResponseCube::npred(const GPhoton&      photon,
                               const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_NPRED_DIFFUSE,
          "Npred computation not implemented for cube-style analysis.");

    // Return Npred
    return 0.0;
}


/***********************************************************************//**
 * @brief Read response
 *
 * @param[in] xml XML element.
 *
 * @exception GException::xml_invalid_parnames
 *            Invalid parameter names found in XML element.
 *
 * Reads information for a CTA observation from an XML element. The exposure
 * and PSF cubes are specified using
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="ExposureCube" file="..."/>
 *       <parameter name="PsfCube"      file="..."/>
 *     </observation>
 *
 ***************************************************************************/
void GCTAResponseCube::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Extract parameters
    int npar[] = {0,0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle ExposureCube
        if (par->attribute("name") == "ExposureCube") {

            // Get filename
            std::string filename = gammalib::strip_whitespace(par->attribute("file"));

            // Load exposure cube
            m_exposure.load(filename);

            // Increment parameter counter
            npar[0]++;
        }

        // Handle PsfCube
        else if (par->attribute("name") == "PsfCube") {

            // Get filename
            std::string filename = gammalib::strip_whitespace(par->attribute("file"));

            // Load PSF cube
            m_psf.load(filename);

            // Increment parameter counter
            npar[1]++;
        }

    } // endfor: looped over observation parameters

    // Verify that all required parameters were found
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::xml_invalid_parnames(G_READ, xml,
              "Require \"ExposureCube\" and \"PsfCube\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response
 *
 * @param[in] xml XML element.
 *
 * Writes information for a CTA observation into an XML element. The exposure
 * and PSF cubes are specified using
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="ExposureCube" file="..."/>
 *       <parameter name="PsfCube"      file="..."/>
 *     </observation>
 *
 ***************************************************************************/
void GCTAResponseCube::write(GXmlElement& xml) const
{
    // Add exposure cube filename
    std::string filename = gammalib::strip_whitespace(m_exposure.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::parameter(xml, "ExposureCube");
        par->attribute("file", filename);
    }

    // Add PSF cube filename
    filename = gammalib::strip_whitespace(m_psf.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::parameter(xml, "PsfCube");
        par->attribute("file", filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA response information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing CTA response information.
 ***************************************************************************/
std::string GCTAResponseCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAResponseCube ===");

        // Append response information
        result.append("\n"+gammalib::parformat("Energy dispersion"));
        if (use_edisp()) {
            result.append("Used");
        }
        else {
            if (apply_edisp()) {
                result.append("Not available");
            }
            else {
                result.append("Not used");
            }
        }

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
void GCTAResponseCube::init_members(void)
{
    // Initialise members
    m_exposure.clear();
    m_psf.clear();
    m_apply_edisp = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GCTAResponseCube::copy_members(const GCTAResponseCube& rsp)
{
    // Copy members
    m_exposure    = rsp.m_exposure;
    m_psf         = rsp.m_psf;
    m_apply_edisp = rsp.m_apply_edisp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponseCube::free_members(void)
{
    // Return
    return;
}
