/***************************************************************************
 *     GCTAResponseCube.hpp - CTA cube analysis response function class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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
#include "GCTAResponse_helpers.hpp"
#include "GCTACubeSourceDiffuse.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAEventBin.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GPhoton.hpp"
#include "GEvent.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GIntegral.hpp"
#include "GObservation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF        "GCTAResponseCube::irf(GEvent&, GPhoton& GObservation&)"
#define G_IRF_PTSRC          "GCTAResponseCube::irf_ptsrc(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_RADIAL        "GCTAResponseCube::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE      "GCTAResponseCube::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_NPRED            "GCTAResponseCube::npred(GPhoton&, GObservation&)"
#define G_READ                         "GCTAResponseCube::read(GXmlElement&)"
#define G_WRITE                       "GCTAResponseCube::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_RADIAL_PSF_BASED  //!< Use Psf-based integration for radial model

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
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Construct CTA response from XML element.
 ***************************************************************************/
GCTAResponseCube::GCTAResponseCube(const GXmlElement& xml) : GCTAResponse()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * @param[in] exposure CTA cube analysis exposure.
 * @param[in] psf CTA cube analysis point spread function.
 *
 * Constructs CTA cube analysis response from a cube analysis exposure and
 * point spread function.
 **************************************************************************/
GCTAResponseCube::GCTAResponseCube(const GCTACubeExposure& exposure,
                                   const GCTACubePsf&      psf) :
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
 * @brief Return instrument response
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 * @return Instrument response.
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

    // Get maximum angular separation for PSF (in radians)
    double delta_max = psf().delta_max();

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
 * @brief Return instrument response
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response.
 *
 * Returns the instrument response for a given event, source and observation.
 ***************************************************************************/
double GCTAResponseCube::irf(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Select IRF depending on the spatial model type
    switch (source.model()->code()) {
        case GMODEL_SPATIAL_POINT_SOURCE:
            irf = irf_ptsrc(event, source, obs);
            break;
        case GMODEL_SPATIAL_RADIAL:
            irf = irf_radial(event, source, obs);
            break;
        case GMODEL_SPATIAL_ELLIPTICAL:
            irf = irf_elliptical(event, source, obs);
            break;
        case GMODEL_SPATIAL_DIFFUSE:
            irf = irf_diffuse(event, source, obs);
            break;
        default:
            break;
    }

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to point source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to point source.
 *
 * Returns the instrument response to a specified point source.
 ***************************************************************************/
double GCTAResponseCube::irf_ptsrc(const GEvent&       event,
                                   const GSource&      source,
                                   const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to model source model
    const GModelSpatialPointSource* ptsrc = static_cast<const GModelSpatialPointSource*>(source.model());

    // Get point source direction
    GSkyDir srcDir = ptsrc->dir();

    // Get pointer on CTA event bin
    if (!event.is_bin()) {
        std::string msg = "The current event is not a CTA event bin. "
                          "This method only works on binned CTA data. Please "
                          "make sure that a CTA observation containing binned "
                          "CTA data is provided.";
        throw GException::invalid_value(G_IRF_PTSRC, msg);
    }
    const GCTAEventBin* bin = static_cast<const GCTAEventBin*>(&event);
    
    // Determine angular separation between true and measured photon
    // direction in radians
    double delta = bin->dir().dir().dist(srcDir);

    // Get maximum angular separation for PSF (in radians)
    double delta_max = psf().delta_max();

    // Compute only if we're sufficiently close to PSF
    if (delta <= delta_max) {

        // Get exposure
        irf = exposure()(srcDir, source.energy());

        // Multiply-in PSF
        if (irf > 0.0) {

            // Revover effective area from exposure
            irf /= obs.livetime();

            // Get PSF component
            irf *= psf()(srcDir, delta, source.energy());

            // Perform deadtime correction
            irf *= obs.deadc(source.time());

        } // endif: exposure was non-zero

    } // endif: we were sufficiently close to PSF

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCTAResponseCube::irf_ptsrc:";
        std::cout << " NaN/Inf encountered";
        std::cout << " irf=" << irf;
        std::cout << std::endl;
    }
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to radial source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to radial source.
 *
 * Returns the instrument response to a specified radial source.
 ***************************************************************************/
double GCTAResponseCube::irf_radial(const GEvent&       event,
                                    const GSource&      source,
                                    const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to CTA event bin
    if (!event.is_bin()) {
        std::string msg = "The current event is not a CTA event bin. "
                          "This method only works on binned CTA data. Please "
                          "make sure that a CTA observation containing binned "
                          "CTA data is provided.";
        throw GException::invalid_value(G_IRF_RADIAL, msg);
    }
    const GCTAEventBin* bin = static_cast<const GCTAEventBin*>(&event);

    // Get event attribute references
    const GSkyDir& obsDir  = bin->dir().dir();
    const GEnergy& obsEng  = bin->energy();
    const GTime&   obsTime = bin->time();

    // Get pointer to radial model
    const GModelSpatialRadial* model = static_cast<const GModelSpatialRadial*>(source.model());

    // Compute angle between model centre and measured photon direction
    // (radians)
    double rho_obs = model->dir().dist(obsDir);

    // Continue only if we're sufficiently close to the model centre to get a
    // non-zero response
    if (rho_obs <= model->theta_max()+psf().delta_max()) {

        // Get exposure
        irf = exposure()(obsDir, obsEng);

        // Continue only if exposure is positive
        if (irf > 0.0) {
            irf /= obs.livetime(); //!< Recover effective area from exposure
            irf *= psf_radial(model, rho_obs, obsDir, obsEng, obsTime);
            irf *= obs.deadc(obsTime);
        }
        
    } // endif: we were sufficiently close

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCTAResponseCube::irf_radial:";
        std::cout << " NaN/Inf encountered";
        std::cout << " irf=" << irf;
        std::cout << std::endl;
    }
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to elliptical source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to elliptical source.
 *
 * Returns the instrument response to a specified elliptical source.
 ***************************************************************************/
double GCTAResponseCube::irf_elliptical(const GEvent&       event,
                                        const GSource&      source,
                                        const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to CTA event bin
    if (!event.is_bin()) {
        std::string msg = "The current event is not a CTA event bin. "
                          "This method only works on binned CTA data. Please "
                          "make sure that a CTA observation containing binned "
                          "CTA data is provided.";
        throw GException::invalid_value(G_IRF_RADIAL, msg);
    }
    const GCTAEventBin* bin = static_cast<const GCTAEventBin*>(&event);

    // Get event attribute references
    const GSkyDir& obsDir  = bin->dir().dir();
    const GEnergy& obsEng  = bin->energy();
    const GTime&   obsTime = bin->time();

    // Get pointer to elliptical model
    const GModelSpatialElliptical* model = static_cast<const GModelSpatialElliptical*>(source.model());

    // Compute angle between model centre and measured photon direction and
    // position angle (radians)
    double rho_obs      = model->dir().dist(obsDir);
    double posangle_obs = model->dir().posang(obsDir);

    // Continue only if we're sufficiently close to the model centre to get a
    // non-zero response
    if (rho_obs <= model->theta_max()+psf().delta_max()) {

        // Get exposure
        irf = exposure()(obsDir, obsEng);

        // Continue only if exposure is positive
        if (irf > 0.0) {
            irf /= obs.livetime(); //!< Recover effective area from exposure
            irf *= psf_elliptical(model, rho_obs, posangle_obs, obsDir, obsEng, obsTime);
            irf *= obs.deadc(obsTime);
        }
        
    } // endif: we were sufficiently close

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCTAResponseCube::irf_elliptical:";
        std::cout << " NaN/Inf encountered";
        std::cout << " irf=" << irf;
        std::cout << std::endl;
    }
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to diffuse source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to diffuse source.
 *
 * Returns the instrument response to a specified diffuse source.
 *
 * The method uses a pre-computation cache to store the instrument response
 * for the spatial model component. The pre-computation cache is initialised
 * if no cache has yet been allocated, or if at the beginning of a scan over
 * the events, the model parameters have changed. The beginning of a scan is
 * defined by an event bin index of 0.
 ***************************************************************************/
double GCTAResponseCube::irf_diffuse(const GEvent&       event,
                                     const GSource&      source,
                                     const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to CTA event bin
    if (!event.is_bin()) {
        std::string msg = "The current event is not a CTA event bin. "
                          "This method only works on binned CTA data. Please "
                          "make sure that a CTA observation containing binned "
                          "CTA data is provided.";
        throw GException::invalid_value(G_IRF_DIFFUSE, msg);
    }
    const GCTAEventBin* bin = static_cast<const GCTAEventBin*>(&event);

    // Get pointer to source cache. We first search the cache for a model
    // with the source name. If no model was found we initialise a new
    // cache entry for that model. Otherwise, we simply return the actual
    // cache entry.
    GCTACubeSourceDiffuse* cache(NULL);
    int index = cache_index(source.name());
    if (index == -1) {
    
        // No cache entry was found, thus allocate and initialise a new one
        cache = new GCTACubeSourceDiffuse;
        cache->set(source.name(), *source.model(), obs);
        m_cache.push_back(cache);

    } // endif: no cache entry was found
    else {
    
        // Check that the cache entry is of the expected type
        if (m_cache[index]->code() != GCTA_CUBE_SOURCE_DIFFUSE) {
            std::string msg = "Cached model \""+source.name()+"\" is not "
                              "an extended source model. This method only "
                              "applies to extended source models.";
            throw GException::invalid_value(G_IRF_DIFFUSE, msg);
        }
        cache = static_cast<GCTACubeSourceDiffuse*>(m_cache[index]);

    } // endelse: there was a cache entry for this model

    // Determine IRF value
    irf = cache->irf(bin->ipix(), bin->ieng());

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCTAResponseCube::irf_diffuse:";
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
 * @return Spatial integral of point spread function. 
 *
 * @exception GException::feature_not_implemented
 *            Method not implemented.
 *
 * This method is a dummy method that is required for a class derived from
 * GResponse but that needs not to be implemented for a binned response.
 ***************************************************************************/
double GCTAResponseCube::npred(const GPhoton&      photon,
                               const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_NPRED,
          "Npred computation not implemented for cube analysis.");

    // Return Npred
    return 0.0;
}


/***********************************************************************//**
 * @brief Read response information from Xml element
 *
 * @param[in] xml XML element.
 *
 * Reads response information from an XML element. The Exposure and Psf cubes
 * are specified using
 *
 *      <observation name="..." id="..." instrument="...">
 *        ...
 *        <parameter name="ExposureCube" file="..."/>
 *        <parameter name="PsfCube"      file="..."/>
 *      </observation>
 *
 ***************************************************************************/
void GCTAResponseCube::read(const GXmlElement& xml)
{
    // Clear response
    clear();

    // Get exposure cube information and load cube
    const GXmlElement* exppar  = gammalib::xml_getpar(G_READ, xml, "ExposureCube");
    std::string        expname = gammalib::strip_whitespace(exppar->attribute("file"));

    // Get PSF cube information and load cube
    const GXmlElement* psfpar  = gammalib::xml_getpar(G_READ, xml, "PsfCube");
    std::string        psfname = gammalib::strip_whitespace(psfpar->attribute("file"));

    // Load cubes
    m_exposure.load(expname);
    m_psf.load(psfname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response information into Xml element
 *
 * @param[in] xml XML element.
 *
 * Writes response information into an Xml element. The Exposure and Psf
 * cubes are specified using
 *
 *      <observation name="..." id="..." instrument="...">
 *        ...
 *        <parameter name="ExposureCube" file="..."/>
 *        <parameter name="PsfCube"      file="..."/>
 *      </observation>
 *
 ***************************************************************************/
void GCTAResponseCube::write(GXmlElement& xml) const
{
    // Add exposure cube filename
    std::string filename = gammalib::strip_whitespace(m_exposure.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::xml_needpar(G_WRITE, xml, "ExposureCube");
        par->attribute("file", filename);
    }

    // Add PSF cube filename
    filename = gammalib::strip_whitespace(m_psf.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::xml_needpar(G_WRITE, xml, "PsfCube");
        par->attribute("file", filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing response information.
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

    // Initialise cache
    m_cache.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response.
 ***************************************************************************/
void GCTAResponseCube::copy_members(const GCTAResponseCube& rsp)
{
    // Copy members
    m_exposure    = rsp.m_exposure;
    m_psf         = rsp.m_psf;
    m_apply_edisp = rsp.m_apply_edisp;

    // Copy cache
    for (int i = 0; i < rsp.m_cache.size(); ++i) {
        m_cache.push_back((rsp.m_cache[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponseCube::free_members(void)
{
    // Free cache
    for (int i = 0; i < m_cache.size(); ++i) {
        if (m_cache[i] != NULL) delete m_cache[i];
        m_cache[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Determine cache index for model
 *
 * @param[in] name Model name.
 * @return Cache index (-1 if model has not been found).
 *
 * Determines the cache index for a given model @p name.
 ***************************************************************************/
int GCTAResponseCube::cache_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Continue only if there are models in cache
    if (!m_cache.empty()) {

         // Search for model name
         for (int i = 0; i < m_cache.size(); ++i) {
             if (m_cache[i]->name() == name) {
                 index = i;
                 break;
             }
         }

    } // endif: there were models in cache

    // Return index
    return index;
}


#if defined(G_RADIAL_PSF_BASED)
/***********************************************************************//**
 * @brief Integrate Psf over radial model
 *
 * @param[in] model Radial model.
 * @param[in] delta_mod Angle between model centre and measured photon
 *                      direction (radians).
 * @param[in] obsDir Observed event direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * Integrates the product of the spatial model and the point spread
 * function over the true photon arrival direction using
 * 
 * \f[
 *    \int_{\delta_{\rm min}}^{\delta_{\rm max}}
 *    \sin \delta \times {\rm Psf}(\delta) \times
 *    \int_{\phi_{\rm min}}^{\phi_{\rm max}} 
 *    S_{\rm p}(\delta, \phi | E, t) d\phi d\delta
 * \f]
 *
 * where
 * \f$S_{\rm p}(\delta, \phi | E, t)\f$ is the radial spatial model,
 * \f${\rm Psf}(\delta)\f$ is the point spread function,
 * \f$\delta\f$ is angular distance between the true and the measured
 * photon direction, and
 * \f$\phi\f$ is the position angle around the observed photon direction
 * measured counterclockwise from the connecting line between the model
 * centre and the observed photon arrival direction.
 ***************************************************************************/
double GCTAResponseCube::psf_radial(const GModelSpatialRadial* model,
                                    const double&              delta_mod,
                                    const GSkyDir&             obsDir,
                                    const GEnergy&             srcEng,
                                    const GTime&               srcTime) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1291
    static const int iter_delta = 5;
    static const int iter_phi   = 6;

    // Initialise value
    double value = 0.0;

    // Get maximum Psf radius (radians)
    double psf_max = psf().delta_max();

    // Get maximum model radius (radians)
    double theta_max = model->theta_max();

    // Set offset angle integration range (radians)
    double delta_min = (delta_mod > theta_max) ? delta_mod - theta_max : 0.0;
    double delta_max = delta_mod + theta_max;
    if (delta_max > psf_max) {
        delta_max = psf_max;
    }

    // Setup integration kernel. We take here the observed photon arrival
    // direction as the true photon arrival direction because the PSF does
    // not vary significantly over a small region.
    cta_psf_radial_kern_delta integrand(this,
                                        model,
                                        obsDir,
                                        srcEng,
                                        srcTime,
                                        delta_mod,
                                        theta_max,
                                        iter_phi);

    // Integrate over model's zenith angle
    GIntegral integral(&integrand);
    integral.fixed_iter(iter_delta);

    // Setup integration boundaries
    std::vector<double> bounds;
    bounds.push_back(delta_min);
    bounds.push_back(delta_max);

    // If the integration range includes a transition between full
    // containment of Psf within model and partial containment, then
    // add a boundary at this location
    double transition_point = theta_max - delta_mod;
    if (transition_point > delta_min && transition_point < delta_max) {
        bounds.push_back(transition_point);
    }

    // Integrate kernel
    value = integral.romberg(bounds, iter_delta);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAResponseCube::psf_radial:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", delta_min=" << delta_min;
        std::cout << ", delta_max=" << delta_max << ")";
        std::cout << std::endl;
    }
    #endif

    // Return integral
    return value;
}
#else
/***********************************************************************//**
 * @brief Integrate Psf over radial model
 *
 * @param[in] model Radial model.
 * @param[in] rho_obs Angle between model centre and measured photon direction
 *                    (radians).
 * @param[in] obsDir Observed event direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * Integrates the product of the spatial model and the point spread
 * function over the true photon arrival direction using
 * 
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *    \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *    {\rm Psf}(\rho, \omega) d\omega d\rho
 * \f]
 *
 * where
 * \f$S_{\rm p}(\rho | E, t)\f$ is the radial spatial model,
 * \f${\rm Psf}(\rho, \omega)\f$ is the point spread function,
 * \f$\rho\f$ is the radial distance from the model centre, and
 * \f$\omega\f$ is the position angle around the model centre measured
 * counterclockwise from the connecting line between the model centre and
 * the observed photon arrival direction.
 ***************************************************************************/
double GCTAResponseCube::psf_radial(const GModelSpatialRadial* model,
                                    const double&              rho_obs,
                                    const GSkyDir&             obsDir,
                                    const GEnergy&             srcEng,
                                    const GTime&               srcTime) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1299
    static const int iter_rho = 5;
    static const int iter_phi = 5;

    // Initialise value
    double irf = 0.0;

    // Get maximum PSF radius (radians)
    double delta_max = psf().delta_max();

    // Set zenith angle integration range for radial model (radians)
    double rho_min = (rho_obs > delta_max) ? rho_obs - delta_max : 0.0;
    double rho_max = rho_obs + delta_max;
    double src_max = model->theta_max();
    if (rho_max > src_max) {
        rho_max = src_max;
    }

    // Perform zenith angle integration if interval is valid
    if (rho_max > rho_min) {

        // Setup integration kernel. We take here the observed photon arrival
        // direction as the true photon arrival direction because the PSF does
        // not vary significantly over a small region.
        cta_psf_radial_kern_rho integrand(this,
                                          model,
                                          obsDir,
                                          srcEng,
                                          srcTime,
                                          rho_obs,
                                          delta_max,
                                          iter_phi);

        // Integrate over model's zenith angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter_rho);

        // Setup integration boundaries
        std::vector<double> bounds;
        bounds.push_back(rho_min);
        bounds.push_back(rho_max);

        // If the integration range includes a transition between full
        // containment of model within Psf and partial containment, then
        // add a boundary at this location
        double transition_point = delta_max - rho_obs;
        if (transition_point > rho_min && transition_point < rho_max) {
            bounds.push_back(transition_point);
        }

        // If we have a shell model then add an integration boundary for the
        // shell radius as a function discontinuity will occur at this
        // location
        const GModelSpatialRadialShell* shell = dynamic_cast<const GModelSpatialRadialShell*>(model);
        if (shell != NULL) {
            double shell_radius = shell->radius() * gammalib::deg2rad;
            if (shell_radius > rho_min && shell_radius < rho_max) {
                bounds.push_back(shell_radius);
            }
        }

        // Integrate kernel
        irf = integral.romberg(bounds, iter_rho);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GCTAResponseCube::psf_radial:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", rho_min=" << rho_min;
            std::cout << ", rho_max=" << rho_max << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: integration interval is valid

    // Return PSF
    return irf;
}
#endif


/***********************************************************************//**
 * @brief Integrate Psf over elliptical model
 *
 * @param[in] model Elliptical model.
 * @param[in] rho_obs Angle between model centre and measured photon direction
 *                    (radians).
 * @param[in] posangle_obs Position angle of measured photon direction with
 *                         respect to model centre (radians).
 * @param[in] obsDir Observed event direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * Integrates the product of the spatial model and the point spread
 * function over the true photon arrival direction using
 * 
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times 
 *    \int_{\omega} 
 *    S_{\rm p}(\rho, \omega | E, t) \times
 *    {\rm Psf}(\rho, \omega) d\omega d\rho
 * \f]
 *
 * where
 * \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical spatial model,
 * \f${\rm Psf}(\rho, \omega)\f$ is the point spread function,
 * \f$\rho\f$ is the radial distance from the model centre, and
 * \f$\omega\f$ is the position angle around the model centre measured
 * counterclockwise from the connecting line between the model centre and
 * the observed photon arrival direction.
 ***************************************************************************/
double GCTAResponseCube::psf_elliptical(const GModelSpatialElliptical* model,
                                        const double&                  rho_obs,
                                        const double&                  posangle_obs,
                                        const GSkyDir&                 obsDir,
                                        const GEnergy&                 srcEng,
                                        const GTime&                   srcTime) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1299
    static const int iter_rho = 5;
    static const int iter_phi = 5;

    // Initialise value
    double irf = 0.0;

    // Get maximum PSF radius (radians)
    double delta_max = psf().delta_max();

    // Get the ellipse boundary (radians). Note that these are NOT the
    // parameters of the ellipse but of a boundary ellipse that is used
    // for computing the relevant omega angle intervals for a given angle
    // rho. The boundary ellipse takes care of the possibility that the
    // semiminor axis is larger than the semimajor axis
    double semimajor;    // Will be the larger axis
    double semiminor;    // Will be the smaller axis
    double posangle;     // Will be the corrected position angle
    double aspect_ratio; // Ratio between smaller/larger axis of model
    if (model->semimajor() >= model->semiminor()) {
        aspect_ratio = (model->semimajor() > 0.0) ?
                        model->semiminor() / model->semimajor() : 0.0;
        posangle     = model->posangle() * gammalib::deg2rad;
    }
    else {
        aspect_ratio = (model->semiminor() > 0.0) ?
                        model->semimajor() / model->semiminor() : 0.0;
        posangle     = model->posangle() * gammalib::deg2rad + gammalib::pihalf;
    }
    semimajor = model->theta_max();
    semiminor = semimajor * aspect_ratio;

    // Set zenith angle integration range for elliptical model
    double rho_min = (rho_obs > delta_max) ? rho_obs - delta_max : 0.0;
    double rho_max = rho_obs + delta_max;
    if (rho_max > semimajor) {
        rho_max = semimajor;
    }

    // Perform zenith angle integration if interval is valid
    if (rho_max > rho_min) {

        // Setup integration kernel. We take here the observed photon arrival
        // direction as the true photon arrival direction because the PSF does
        // not vary significantly over a small region.
        cta_psf_elliptical_kern_rho integrand(this,
                                              model,
                                              semimajor,
                                              semiminor,
                                              posangle,
                                              obsDir,
                                              srcEng,
                                              srcTime,
                                              rho_obs,
                                              posangle_obs,
                                              delta_max,
                                              iter_phi);

        // Integrate over model's zenith angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter_rho);

        // Setup integration boundaries
        std::vector<double> bounds;
        bounds.push_back(rho_min);
        bounds.push_back(rho_max);

        // Kluge: add this transition point as this allows to fit the test
        // case without any stalls. Not clear why this is the case, maybe
        // simply because the rho integral gets cut down into one more
        // sub-interval which may increase precision and smoothed the
        // likelihood contour
        double transition_point = delta_max - rho_obs;
        if (transition_point > rho_min && transition_point < rho_max) {
            bounds.push_back(transition_point);
        }

        // If the integration range includes the semiminor boundary, then
        // add an integration boundary at that location
        if (semiminor > rho_min && semiminor < rho_max) {
            bounds.push_back(semiminor);
        }

        // Integrate kernel
        irf = integral.romberg(bounds, iter_rho);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GCTAResponseCube::psf_elliptical:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", rho_min=" << rho_min;
            std::cout << ", rho_max=" << rho_max << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: integration interval is valid

    // Return PSF
    return irf;
}
