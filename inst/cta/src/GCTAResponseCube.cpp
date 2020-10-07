/***************************************************************************
 *     GCTAResponseCube.cpp - CTA cube analysis response function class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2020 by Juergen Knoedlseder                         *
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
#include "GPhoton.hpp"
#include "GSource.hpp"
#include "GEvent.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GEnergies.hpp"
#include "GTime.hpp"
#include "GIntegral.hpp"
#include "GIntegrals.hpp"
#include "GObservation.hpp"
#include "GNdarray.hpp"
#include "GMatrixSparse.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialComposite.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GCTAResponseCube.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAEventBin.hpp"
#include "GCTASupport.hpp"
#include "GCTAEventCube.hpp"               // Kludge
#include "GCTAEventBin.hpp"                // Kludge

/* __ Method name definitions ____________________________________________ */
#define G_IRF        "GCTAResponseCube::irf(GEvent&, GPhoton& GObservation&)"
#define G_IRF_PTSRC          "GCTAResponseCube::irf_ptsrc(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_RADIAL        "GCTAResponseCube::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE      "GCTAResponseCube::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_NROI        "GCTAResponseCube::nroi(GModelSky&, GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_EBOUNDS                       "GCTAResponseCube::ebounds(GEnergy&)"
#define G_READ                         "GCTAResponseCube::read(GXmlElement&)"
#define G_WRITE                       "GCTAResponseCube::write(GXmlElement&)"
#define G_IRF_RADIAL2             "GCTAResponseCube::irf_radial(GModelSky&, "\
                                             "GObservation&, GMatrixSparse*)"

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
 * Constructs CTA cube response by making a deep copy of an existing
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
 * @param[in] background CTA cube background response.
 *
 * Constructs CTA cube analysis response from a cube analysis exposure,
 * a point spread function cube and a background cube.
 **************************************************************************/
GCTAResponseCube::GCTAResponseCube(const GCTACubeExposure&   exposure,
                                   const GCTACubePsf&        psf,
                                   const GCTACubeBackground& background) :
                  GCTAResponse()
{
    // Initialise members
    init_members();

    // Set members
    m_exposure   = exposure;
    m_psf        = psf;
    m_background = background;

    // Signal that no energy dispersion was given
    m_has_edisp = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * @param[in] exposure CTA cube analysis exposure.
 * @param[in] psf CTA cube analysis point spread function.
 * @param[in] edisp CTA cube energy dispersion response.
 * @param[in] background CTA cube background response.
 *
 * Constructs CTA cube analysis response from a cube analysis exposure,
 * a point spread function cube, an energy dispersion cube and a background cube.
 **************************************************************************/
GCTAResponseCube::GCTAResponseCube(const GCTACubeExposure&   exposure,
                                   const GCTACubePsf&        psf,
                                   const GCTACubeEdisp&      edisp,
                                   const GCTACubeBackground& background) :
                  GCTAResponse()
{
    // Initialise members
    init_members();

    // Set members
    m_exposure   = exposure;
    m_psf        = psf;
    m_edisp      = edisp;
    m_background = background;

    // Signal that edisp is available
    m_has_edisp = true;

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
 * @param[in] obs Observation (not used).
 * @return Instrument response.
 ***************************************************************************/
double GCTAResponseCube::irf(const GEvent&       event,
                             const GPhoton&      photon,
                             const GObservation& obs) const
{
    // Retrieve event instrument direction
    const GCTAInstDir& dir = gammalib::cta_dir(G_IRF, event);

    // Get event attributes
    const GSkyDir& obsDir = dir.dir();
    const GEnergy& obsEng = event.energy();

    // Get photon attributes
    const GSkyDir& srcDir  = photon.dir();
    const GEnergy& srcEng  = photon.energy();
    //const GTime&   srcTime = photon.time();

    // Determine angular separation between true and measured photon
    // direction in radians
    double delta = obsDir.dist(srcDir);

    // Get maximum angular separation for PSF (in radians)
    double delta_max = psf().delta_max();

    // Initialise IRF value
    double irf = 0.0;

    // Get livetime (in seconds)
    double livetime = exposure().livetime();

    // Continue only if livetime is >0 and if we're sufficiently close
    // to the PSF
    if ((livetime > 0.0) && (delta <= delta_max)) {

        // Get exposure
        irf = exposure()(srcDir, srcEng);

        // Multiply-in PSF
        if (irf > 0.0) {

            // Get PSF component
            irf *= psf()(srcDir, delta, srcEng);

            // Multiply-in energy dispersion
            if (use_edisp() && irf > 0.0) {
                irf *= edisp()(obsEng, srcEng, srcDir);
            }

            // Divide by livetime
            irf /= livetime;

            // Apply deadtime correction
            irf *= exposure().deadc();

        } // endif: Aeff was non-zero

    } // endif: we were sufficiently close to PSF and livetime was >0

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
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t') = \int_{\rm ROI} P(p',E',t') dp'
 * \f]
 *
 * of the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * for a given sky model \f$S(p,E,t)\f$ and response function
 * \f$R(p',E',t'|p,E,t)\f$ over the Region of Interest (ROI).
 *
 * @todo Implement method (is maybe not really needed)
 ***************************************************************************/
double GCTAResponseCube::nroi(const GModelSky&    model,
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
 * @param[in] obsEng Observed photon energy.
 * @return Boundaries in true energy.
 ***************************************************************************/
GEbounds GCTAResponseCube::ebounds(const GEnergy& obsEng) const
{
    // Get energy boundaries from energy dispersion
    GEbounds ebounds = m_edisp.ebounds(obsEng);

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Return instrument response integrated over the spatial model
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to a spatial model.
 *
 * Returns the instrument response for a given event, source and observation
 * integrated over the spatial model component. The method computes
 *
 * \f[
 *    {\tt irf}(p', E', t') = \int_p M_{\rm S}(p | E, t) \,
 *                                   R(p', E', t' | p, E, t) \, d\,p
 * \f]
 *
 * where
 * * \f$M_{\rm S}(p | E, t)\f$ is the spatial model component,
 * * \f$R(p', E', t' | p, E, t)\f$ is the Instrument Response Function (IRF),
 * * \f$p'\f$ is the measured instrument direction,
 * * \f$E'\f$ is the measured or reconstructed energy,
 * * \f$t'\f$ is the measured arrival time,
 * * \f$p\f$ is the true photon arrival direction,
 * * \f$E\f$ is the true photon energy, and
 * * \f$t\f$ is the true trigger time.
 *
 * The integration is done over all relevant true sky directions \f$p\f$.
 *
 * Depending on the type of the source model the method branches to the
 * following methods to perform the actual computations
 *
 *      irf_ptsrc() - for the handling of a point source
 *      irf_radial() - for radial models
 *      irf_elliptical() - for elliptical models
 *      irf_diffuse() - for diffuse models
 *      irf_composite() - for composite models
 *
 * The method implements a caching mechanism for spatial models that have all
 * parameters fixed. For those models the instrument response for a given
 * event and observation is only computed once and then stored in an internal
 * cache from which it is fetched back in case that the method is called
 * again for the same event and observation.
 ***************************************************************************/
double GCTAResponseCube::irf_spatial(const GEvent&       event,
                                     const GSource&      source,
                                     const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Set IRF value attributes
    std::string     name  = obs.id() + "::" + source.name();
    const GInstDir& dir   = event.dir();
    const GEnergy&  ereco = event.energy();
    const GEnergy&  etrue = source.energy();

    // Signal if spatial model has free parameters
    bool has_free_pars = source.model()->has_free_pars();

    // Kludge: the IRF response cache should be used, the model has no
    // free parameters and there is no energy dispersion
    if (m_use_irf_cache && !has_free_pars && !use_edisp()) {

        // Build unique cache name
        std::string name  = obs.id() + "::" + source.name();

        // Get index in cache, returns -1 if name is not found in cache
        int index = -1;
        for (int i = 0; i < m_cache_names.size(); ++i) {
            if (m_cache_names[i] == name) {
                index = i;
                break;
            }
        }

        // If index was not found then allocate a new cache map
        if (index == -1) {

            // Get pointer to event cube
            const GCTAEventCube* cube =
                  static_cast<const GCTAEventCube*>(obs.events());

            // Allocate cache
            GNdarray cache(cube->nx()*cube->ny(), cube->ebins());

            // Initialise all cache values with -1 (not set)
            for (int i = 0; i < cache.size(); ++i) {
                cache(i) = -1.0;
            }

            // Insert cache
            m_cache_names.push_back(name);
            m_cache_values.push_back(cache);

            // Set index
            index = m_cache_names.size()-1;

        } // endif: allocated new cache

        // Get reference to CTA event bin
        const GCTAEventBin& bin = static_cast<const GCTAEventBin&>(event);

        // Get cache value
        irf = m_cache_values[index](bin.ipix(), bin.ieng());

        // If cache value is not valid then copute IRF
        if (irf < 0.0) {

            // Compute IRF for spatial model
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
                case GMODEL_SPATIAL_COMPOSITE:
                    irf = irf_composite(event, source, obs);
                    break;
                default:
                    break;
            }

            // Set cache value
            m_cache_values[index](bin.ipix(), bin.ieng()) = irf;

        } // endif: computed IRF

    } // endif: kludge

    // ... otherwise use release 1.7 response cache
    else {

        // If the spatial model component has free parameters, or the response
        // cache should not be used, or the cache does not contain the requested
        // IRF value then compute the IRF value for the spatial model.
        if (has_free_pars    ||
            !m_use_irf_cache ||
            !m_irf_cache.contains(name, dir, ereco, etrue, &irf)) {

            // Compute IRF for spatial model
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
                case GMODEL_SPATIAL_COMPOSITE:
                    irf = irf_composite(event, source, obs);
                    break;
                default:
                    break;
            }

        } // endif: computed spatial model

        // If the spatial model has no free parameters and the response cache
        // should be used then put the IRF value in the response cache.
        if (!has_free_pars && m_use_irf_cache) {
            m_irf_cache.set(name, dir, ereco, etrue, irf);
        }

    } // endelse: used release 1.7 response cache

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Read response information from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads response information from an XML element. The Exposure, Psf,
 * background cubes, and optionally the energy dispersion cube, are specified
 * using
 *
 *      <observation name="..." id="..." instrument="...">
 *        ...
 *        <parameter name="ExposureCube" file="..."/>
 *        <parameter name="PsfCube"      file="..."/>
 *        <parameter name="EdispCube"    file="..."/>
 *        <parameter name="BkgCube"      file="..."/>
 *      </observation>
 *
 ***************************************************************************/
void GCTAResponseCube::read(const GXmlElement& xml)
{
    // Clear response
    clear();

    // Get exposure cube information and load cube
    const GXmlElement* exppar  = gammalib::xml_get_par(G_READ, xml, "ExposureCube");
    std::string        expname = gammalib::xml_file_expand(xml,
                                           exppar->attribute("file"));

    // Get PSF cube information and load cube
    const GXmlElement* psfpar  = gammalib::xml_get_par(G_READ, xml, "PsfCube");
    std::string        psfname = gammalib::xml_file_expand(xml,
                                           psfpar->attribute("file"));

    // Get background cube information and load cube
    const GXmlElement* bkgpar  = gammalib::xml_get_par(G_READ, xml, "BkgCube");
    std::string        bkgname = gammalib::xml_file_expand(xml,
                                           bkgpar->attribute("file"));

    // Load cubes
    m_exposure.load(expname);
    m_psf.load(psfname);
    m_background.load(bkgname);

    // Optionally load energy dispersion cube
    if (gammalib::xml_has_par(xml, "EdispCube")) {
        const GXmlElement* edisppar  = gammalib::xml_get_par(G_READ, xml, "EdispCube");
        std::string        edispname = gammalib::xml_file_expand(xml,
                                                 edisppar->attribute("file"));
        m_edisp.load(edispname);
        m_has_edisp = true;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response information into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes response information into an XML element. The Exposure, Psf
 * and background cubes, and optionally the energy dispersion cube, are
 * specified using
 *
 *      <observation name="..." id="..." instrument="...">
 *        ...
 *        <parameter name="ExposureCube" file="..."/>
 *        <parameter name="PsfCube"      file="..."/>
 *        <parameter name="EdispCube"    file="..."/>
 *        <parameter name="BkgCube"      file="..."/>
 *      </observation>
 *
 ***************************************************************************/
void GCTAResponseCube::write(GXmlElement& xml) const
{
    // Add exposure cube filename
    std::string filename = gammalib::xml_file_reduce(xml, m_exposure.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "ExposureCube");
        par->attribute("file", filename);
    }

    // Add PSF cube filename
    filename = gammalib::xml_file_reduce(xml, m_psf.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "PsfCube");
        par->attribute("file", filename);
    }

    // Optionally add energy dispersions cube filename
    if (m_has_edisp) {
        filename = gammalib::xml_file_reduce(xml, m_edisp.filename());
        if (!(filename.empty())) {
            GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "EdispCube");
            par->attribute("file", filename);
        }
    }

    // Add background cube filename
    filename = gammalib::xml_file_reduce(xml, m_background.filename());
    if (!(filename.empty())) {
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "BkgCube");
        par->attribute("file", filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response information
 *
 * @param[in] chatter Chattiness.
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

        // Get reduced chatter level
        GChatter reduced_chatter = gammalib::reduce(chatter);

        // Append detailed information
        if (chatter >= NORMAL) {

            // Append exposure cube information
            result.append("\n"+m_exposure.print(reduced_chatter));

            // Append point spread function information
            result.append("\n"+m_psf.print(reduced_chatter));

            // Optionally append energy dispersion information
            if (m_has_edisp) {
                result.append("\n"+m_edisp.print(reduced_chatter));
            }

            // Append background information
            result.append("\n"+m_background.print(reduced_chatter));

            // Append cache information
            result.append("\n"+m_irf_cache.print(reduced_chatter));

        } // endif: appended detailed information

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
    m_edisp.clear();
    m_background.clear();
    m_apply_edisp    = false;
    m_has_edisp      = false;

    // Kludge: Initialise cube response cache
    m_cache_names.clear();
    m_cache_values.clear();

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
    m_edisp       = rsp.m_edisp;
    m_background  = rsp.m_background;
    m_apply_edisp = rsp.m_apply_edisp;
    m_has_edisp   = rsp.m_has_edisp;

    // Kludge: Copy cube response cache
    m_cache_names  = rsp.m_cache_names;
    m_cache_values = rsp.m_cache_values;

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


/***********************************************************************//**
 * @brief Integrate PSF over diffuse model
 *
 * @param[in] model Diffuse spatial model.
 * @param[in] obsDir Observed event direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * Computes the integral
 * 
 * \f[
 *    \int_0^{\delta_{\rm max}}
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm Map}(\delta, \phi) \sin \delta
 *    {\rm d}\phi {\rm d}\delta
 * \f]
 *
 * where \f${\rm Map}(\delta, \phi)\f$ is the diffuse map in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 * \f${\rm PSF}(\delta)\f$ is the azimuthally symmetric point spread
 * function.
 ***************************************************************************/
double GCTAResponseCube::psf_diffuse(const GModelSpatial* model,
                                     const GSkyDir&       obsDir,
                                     const GEnergy&       srcEng,
                                     const GTime&         srcTime) const
{
    // Set minimum and maximum number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/2458
    static const int min_iter_delta = 5;
    static const int min_iter_phi   = 5;
    static const int max_iter_delta = 8;
    static const int max_iter_phi   = 8;

    // Compute rotation matrix to convert from PSF centred coordinate system
    // spanned by delta and phi into the reference frame of the observed
    // arrival direction given in Right Ascension and Declination.
    GMatrix ry;
    GMatrix rz;
    ry.eulery(obsDir.dec_deg() - 90.0);
    rz.eulerz(-obsDir.ra_deg());
    GMatrix rot = (ry * rz).transpose();

    // Get offset angle integration interval in radians
    double delta_min = 0.0;
    double delta_max = 1.1 * this->psf().delta_max();

    // Get resolution of spatial model
    double resolution = gammalib::resolution(model);


    // Setup integration kernel. We take here the observed photon arrival
    // direction as the true photon arrival direction because the PSF does
    // not vary significantly over a small region.
    cta_psf_diffuse_kern_delta integrand(this, model, &rot,
                                         obsDir, srcEng, srcTime,
                                         min_iter_phi, max_iter_phi,
                                         resolution);

    // Set number of radial integration iterations
    int iter  = gammalib::iter_rho(delta_max, resolution,
                                   min_iter_delta, max_iter_delta);

    // Setup integration
    GIntegral integral(&integrand);

    // Set fixed number of iterations
    integral.fixed_iter(iter);

    // Integrate over PSF delta angle
    double psf = integral.romberg(delta_min, delta_max);

    // Return PSF
    return psf;
}


/***********************************************************************//**
 * @brief Return instrument response to point source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation (not used).
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
    const GModelSpatialPointSource* ptsrc =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Get point source direction
    GSkyDir srcDir = ptsrc->dir();

    // Get energy of source model
    GEnergy srcEng = source.energy();

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

    // Get livetime (in seconds)
    double livetime = exposure().livetime();

    // Continue only if livetime is >0 and if we're sufficiently close
    // to the PSF
    if ((livetime > 0.0) && (delta <= delta_max)) {

        // Get exposure
        irf = exposure()(srcDir, srcEng);

        // Multiply-in PSF
        if (irf > 0.0) {

            // Recover effective area from exposure
            irf /= livetime;

            // Get PSF component
            irf *= psf()(srcDir, delta, srcEng);

            // Multiply-in energy dispersion
            if (use_edisp() && irf > 0.0) {
                irf *= edisp()(bin->energy(), srcEng, srcDir);
            }

            // Apply deadtime correction
            irf *= exposure().deadc();

        } // endif: exposure was non-zero

    } // endif: we were sufficiently close to PSF and livetime >0

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
 * @param[in] obs Observation (not used).
 * @return Instrument response to radial source.
 *
 * Returns the instrument response to a specified radial source.
 *
 * @todo Correct assumptions in exposure computation, PSF computation and
 *       energy dispersion computation.
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

    // Get energy of source model
    GEnergy srcEng = source.energy();

    // Get pointer to radial model
    const GModelSpatialRadial* model =
          static_cast<const GModelSpatialRadial*>(source.model());

    // Compute angle between model centre and measured photon direction
    // (radians)
    double rho_obs = model->dir().dist(obsDir);

    // Get livetime (in seconds)
    double livetime = exposure().livetime();

    // Continue only if livetime is >0 and if we're sufficiently close to
    // the model centre to get a non-zero response
    if ((livetime > 0.0) && (rho_obs <= model->theta_max()+psf().delta_max())) {

        // Get exposure at the observed event direction.
        //
        // The current code assumes that the exposure at the observed and
        // true event direction does not vary significantly. In other words,
        // the code assumes that the exposure is constant over the size of
        // the PSF.
        irf = exposure()(obsDir, srcEng);

        // Continue only if exposure is positive
        if (irf > 0.0) {

            // Recover effective area from exposure
            irf /= livetime;

            // Get PSF component
            irf *= psf_radial(model, rho_obs, obsDir, srcEng, obsTime);

            // Multiply-in energy dispersion
            //
            // The current code assumes that the energy dispersion at the
            // observed and true event direction does not vary
            // significantly. In other words, the code assumes that the
            // energy dispersion is constant over the size of the PSF.
            if (use_edisp() && irf > 0.0) {
                irf *= edisp()(bin->energy(), srcEng, obsDir);
            }

            // Apply deadtime correction
            irf *= exposure().deadc();

        } // endif: exposure was positive

    } // endif: we were sufficiently close and livetime >0

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
 * @param[in] obs Observation (not used).
 * @return Instrument response to elliptical source.
 *
 * Returns the instrument response to a specified elliptical source.
 *
 * @todo Correct assumptions in exposure computation, PSF computation and
 *       energy dispersion computation.
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

    // Get energy of source model
    GEnergy srcEng = source.energy();

    // Get pointer to elliptical model
    const GModelSpatialElliptical* model =
          static_cast<const GModelSpatialElliptical*>(source.model());

    // Compute angle between model centre and measured photon direction and
    // position angle (radians)
    double rho_obs      = model->dir().dist(obsDir);
    double posangle_obs = model->dir().posang(obsDir); // Celestial

    // Get livetime (in seconds)
    double livetime = exposure().livetime();

    // Continue only if livetime is >0 and if we're sufficiently close to
    // the model centre to get a non-zero response
    if ((livetime > 0.0) && (rho_obs <= model->theta_max()+psf().delta_max())) {

        // Get exposure
        //
        // The current code assumes that the exposure at the observed and
        // true event direction does not vary significantly. In other words,
        // the code assumes that the exposure is constant over the size of
        // the PSF.
        irf = exposure()(obsDir, srcEng);

        // Continue only if exposure is positive
        if (irf > 0.0) {

            // Recover effective area from exposure
            irf /= livetime;

            // Get PSF component
            irf *= psf_elliptical(model, rho_obs, posangle_obs, obsDir, srcEng, obsTime);

            // Multiply-in energy dispersion
            //
            // The current code assumes that the energy dispersion at the
            // observed and true event direction does not vary
            // significantly. In other words, the code assumes that the
            // energy dispersion is constant over the size of the PSF.
            if (use_edisp() && irf > 0.0) {
                irf *= edisp()(bin->energy(), srcEng, obsDir);
            }

            // Apply deadtime correction
            irf *= exposure().deadc();

        } // endif: exposure was positive

    } // endif: we were sufficiently close and livetime >0

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
        throw GException::invalid_argument(G_IRF_DIFFUSE, msg);
    }
    const GCTAEventBin* bin = static_cast<const GCTAEventBin*>(&event);

    // Get event attribute references
    const GSkyDir& obsDir  = bin->dir().dir();
    const GEnergy& obsEng  = bin->energy();
    const GTime&   obsTime = bin->time();

    // Get energy of source model
    GEnergy srcEng = source.energy();

    // Get pointer to elliptical model
    const GModelSpatialDiffuse* model =
          static_cast<const GModelSpatialDiffuse*>(source.model());

    // Get pointer on CTA response cube
    //const GCTAResponseCube* rsp = gammalib::cta_rsp_cube(G_IRF_DIFFUSE, obs);

    // Get livetime (in seconds) and deadtime correction factor
    double livetime = exposure().livetime();
    double deadc    = exposure().deadc();

    // Get Psf radius (in degrees)
    double delta_max = psf().delta_max() * gammalib::rad2deg;

    // Continue only if livetime is >0 and model contains reconstructed
    // sky direction
    if (livetime > 0.0 && model->contains(obsDir, delta_max))  {

        // Get exposure
        //
        // The current code assumes that the exposure at the observed and
        // true event direction does not vary significantly. In other words,
        // the code assumes that the exposure is constant over the size of
        // the PSF.
        irf = exposure()(obsDir, srcEng);

        // Continue only if exposure is positive
        if (irf > 0.0) {

            // Recover effective area from exposure
            irf /= livetime;

            // Compute product of PSF and diffuse map, integrated over the
            // relevant PSF area. We assume no energy dispersion and thus
            // compute the product using the observed energy.
            irf *= psf_diffuse(model, obsDir, srcEng, obsTime);

            // Multiply-in energy dispersion
            //
            // The current code assumes that the energy dispersion at the
            // observed and true event direction does not vary
            // significantly. In other words, the code assumes that the
            // energy dispersion is constant over the size of the PSF.
            if (use_edisp() && irf > 0.0) {
                irf *= edisp()(bin->energy(), srcEng, obsDir);
            }

            // Apply deadtime correction
            irf *= deadc;

        } // endif: exposure was positive

    } // endif: livetime was positive

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


/*==========================================================================
 =                                                                         =
 =                           New private methods                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return instrument response to radial source
 *
 * @param[in] model Sky model.
 * @param[in] obs Observation.
 * @param[out] gradients Optional spatial model gradients for all events.
 * @return Instrument response to radial source for all events in observation.
 *
 * Returns the instrument response to a specified radial source.
 *
 * @p gradients is an optional sparse matrix where the number of rows
 * corresponds to the number of events in the observation and the number
 * of columns corresponds to the number of spatial model parameters.
 ***************************************************************************/
GVector GCTAResponseCube::irf_radial(const GModelSky&    model,
                                     const GObservation& obs,
                                     GMatrixSparse*      gradients) const
{
    // Get number of events
    int nevents = obs.events()->size();

    // Initialise result
    GVector irfs(nevents);

    // Get livetime (in seconds)
    double livetime = exposure().livetime();

    // Continue only if livetime is positive
    if (livetime > 0.0) {

        // Set gradients flag
        bool grad = (gradients != NULL);

        // Get pointer to radial model
        const GModelSpatialRadial* radial =
              static_cast<const GModelSpatialRadial*>(model.spatial());

        // Get CTA event cube
        const GCTAEventCube& cube = gammalib::cta_event_cube(G_IRF_RADIAL2, obs);

        // Get number of radial model parameters and CTA event cube dimensions
        int npars = radial->size();
        int ndirs = cube.npix();
        int nengs = cube.ebins();

        // Check matrix consistency
        if (grad) {
            if (gradients->rows() != nevents) {
                std::string msg = "Number of "+gammalib::str(gradients->rows())+
                                  " rows in gradient matrix differs from number "
                                  "of "+gammalib::str(nevents)+" events in "
                                  "observation. Please provide a compatible "
                                  "gradient matrix.";
                throw GException::invalid_argument(G_IRF_RADIAL2, msg);
            }
            if (gradients->columns() != npars) {
                std::string msg = "Number of "+gammalib::str(gradients->columns())+
                                  " columns in gradient matrix differs from number "
                                  "of "+gammalib::str(npars)+" parameters in "
                                  "model. Please provide a compatible "
                                  "gradient matrix.";
                throw GException::invalid_argument(G_IRF_RADIAL2, msg);
            }
        }

        // Setup energies container
        GEnergies srcEngs;
        for (int ieng = 0; ieng < nengs; ++ieng) {
            srcEngs.append(cube.energy(ieng));
        }

        // If requested, setup vectors of gradients
        GVector* gradient = NULL;
        if (grad) {
            gradient = new GVector[npars];
            for (int ipar = 0; ipar < npars; ++ipar) {
                gradient[ipar] = GVector(nevents);
            }
        }

        // Setup cube time
        GTime srcTime = cube.time();

        // Set IRF normalisation
        double norm = exposure().deadc() / livetime;

        // Loop over event directions
        for (int idir = 0; idir < ndirs; ++idir) {

            // Get event direction
            GSkyDir obsDir = cube.counts().inx2dir(idir);

            // Compute angle between model centre and measured photon direction
            // (radians)
            double zeta = radial->dir().dist(obsDir);

            // Continue only if we're sufficiently close to the model centre to
            // get a non-zero response
            if (zeta <= radial->theta_max()+psf().delta_max()) {

                // Get radially integrated PSF
                GVector psf = psf_radial(radial, zeta, obsDir, srcEngs, srcTime, grad);

                // Loop over true energies
                for (int ieng = 0, index = idir; ieng < nengs; ++ieng, index += ndirs) {

                    // Get effective area at the observed direction. This
                    // assumes that the exposure at the observed and true event
                    // direction does not vary significantly. In other words,
                    // is is assumed that the exposure is constant over the
                    // size of the PSF.
                    double aeff = norm * exposure()(obsDir, srcEngs[ieng]);

                    // Compute IRF value
                    irfs[index] = aeff * psf[ieng];

                    // Optionally compute gradients
                    if (grad) {
                        for (int ipar = 0, ipsf = ieng+nengs; ipar < npars;
                             ++ipar, ipsf += nengs) {
                            gradient[ipar][index] = aeff * psf[ipsf];
                        }
                    }

                } // endfor: looped over energies

            } // endif: direction close enough to model centre

        } // endfor: looped over event directions

        // If gradients were requested then insert vectors into the gradient
        // matrix
        if (grad) {
            for (int ipar = 0; ipar < npars; ++ipar) {
                gradients->add_to_column(ipar, gradient[ipar]);
            }
        }

        // If needed, free vectors of gradients
        if (gradient != NULL) {
            delete [] gradient;
        }

    } // endif: livetime was positive

    // Return IRF value
    return irfs;
}


/***********************************************************************//**
 * @brief Integrate Psf over radial model
 *
 * @param[in] model Radial model.
 * @param[in] zeta Angle between model centre and event direction (radians).
 * @param[in] obsDir Observed event direction.
 * @param[in] srcEngs True photon energies.
 * @param[in] srcTime True photon arrival time.
 * @param[in] grad Compute gradients?
 ***************************************************************************/
GVector GCTAResponseCube::psf_radial(const GModelSpatialRadial* model,
                                     const double&              zeta,
                                     const GSkyDir&             obsDir,
                                     const GEnergies            srcEngs,
                                     const GTime&               srcTime,
                                     const bool&                grad) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1291
    static const int iter_delta = 5;
    static const int iter_phi   = 6;

    // Get maximum Psf radius (radians)
    double psf_max = psf().delta_max();

    // Get maximum model radius (radians)
    double theta_max = model->theta_max();

    // Set offset angle integration range (radians)
    double delta_min = (zeta > theta_max) ? zeta - theta_max : 0.0;
    double delta_max = zeta + theta_max;
    if (delta_max > psf_max) {
        delta_max = psf_max;
    }

    // Setup integration kernel. We take here the observed photon arrival
    // direction as the true photon arrival direction because the PSF does
    // not vary significantly over a small region.
    cta_psf_radial_kerns_delta integrand(this,
                                         model,
                                         obsDir,
                                         srcEngs,
                                         zeta,
                                         theta_max,
                                         iter_phi,
                                         grad);

    // Integrate over model's zenith angle
    GIntegrals integral(&integrand);
    integral.fixed_iter(iter_delta);

    // Setup integration boundaries
    std::vector<double> bounds;
    bounds.push_back(delta_min);
    bounds.push_back(delta_max);

    // If the integration range includes a transition between full
    // containment of Psf within model and partial containment, then
    // add a boundary at this location
    /*
    double transition_point = theta_max - zeta;
    if (transition_point > delta_min && transition_point < delta_max) {
        bounds.push_back(transition_point);
    }
    */

    // Integrate kernel
    GVector values = integral.romberg(bounds, iter_delta);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    for (int i = 0; i < srcEngs.size(); ++i) {
        if (gammalib::is_notanumber(values[i]) ||
            gammalib::is_infinite(values[i])) {
            std::cout << "*** ERROR: GCTAResponseCube::psf_radial:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (value=" << values[i];
            std::cout << ", delta_min=" << delta_min;
            std::cout << ", delta_max=" << delta_max << ")";
            std::cout << std::endl;
        }
    }
    #endif

    // Return integrals
    return values;
}
