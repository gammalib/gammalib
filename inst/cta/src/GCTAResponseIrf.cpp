/***************************************************************************
 *       GCTAResponseIrf.cpp - CTA instrument response function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseIrf.cpp
 * @brief CTA response class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <vector>
#include <string>
#include "GFits.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFilename.hpp"
#include "GIntegral.hpp"
#include "GCaldb.hpp"
#include "GSource.hpp"
#include "GRan.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialRadialRing.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialComposite.hpp"
#include "GArf.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAPointing.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAEventList.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"
#include "GCTAAeff.hpp"
#include "GCTAAeff2D.hpp"
#include "GCTAAeffArf.hpp"
#include "GCTAAeffPerfTable.hpp"
#include "GCTAPsf.hpp"
#include "GCTAPsf2D.hpp"
#include "GCTAPsfVector.hpp"
#include "GCTAPsfPerfTable.hpp"
#include "GCTAPsfKing.hpp"
#include "GCTAPsfTable.hpp"
#include "GCTAEdisp.hpp"
#include "GCTAEdisp2D.hpp"
#include "GCTAEdispRmf.hpp"
#include "GCTAEdispPerfTable.hpp"
#include "GCTABackground.hpp"
#include "GCTABackgroundPerfTable.hpp"
#include "GCTABackground3D.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF   "GCTAResponseIrf::irf(GInstDir&, GEnergy&, GTime&, GSkyDir&,"\
                                          " GEnergy&, GTime&, GObservation&)"
#define G_MC   "GCTAResponseIrf::mc(double&, GPhoton&, GObservation&, GRan&)"
#define G_READ                          "GCTAResponseIrf::read(GXmlElement&)"
#define G_WRITE                        "GCTAResponseIrf::write(GXmlElement&)"
#define G_LOAD_AEFF                  "GCTAResponseIrf::load_aeff(GFilename&)"
#define G_LOAD_PSF                    "GCTAResponseIrf::load_psf(GFilename&)"
#define G_LOAD_EDISP                "GCTAResponseIrf::load_edisp(GFilename&)"
#define G_LOAD_BACKGROUND      "GCTAResponseIrf::load_background(GFilename&)"
#define G_NIRF            "GCTAResponseIrf::nirf(GPhoton&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_IRF_RADIAL         "GCTAResponseIrf::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_ELLIPTICAL "GCTAResponseIrf::irf_elliptical(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE       "GCTAResponseIrf::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_NROI_RADIAL    "GCTAResponseIrf::nroi_radial(GModelSky&, GEnergy&,"\
                                  " GTime&, GEnergy&, GTime&, GObservation&)"
#define G_NROI_ELLIPTICAL      "GCTAResponseIrf::nroi_elliptical(GModelSky&,"\
                        " GEnergy&, GTime&, GEnergy&, GTime&, GObservation&)"
#define G_NROI_DIFFUSE  "GCTAResponseIrf::nroi_diffuse(GModelSky&, GEnergy&,"\
                                  " GTime&, GEnergy&, GTime&, GObservation&)"
#define G_AEFF    "GCTAResponseIrf::aeff(double&, double&, double&, double&,"\
                                                                  " double&)"
#define G_PSF      "GCTAResponseIrf::psf(double&, double&, double&, double&,"\
                                                                  " double&)"
#define G_PSF_DELTA_MAX    "GCTAResponseIrf::psf_delta_max(double&, double&,"\
                                                " double&, double&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_USE_PSF_SYSTEM      //!< Do radial Irf integrations in Psf system

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_IRF_RADIAL                     //!< Debug irf_radial method
//#define G_DEBUG_IRF_DIFFUSE                   //!< Debug irf_diffuse method
//#define G_DEBUG_IRF_ELLIPTICAL             //!< Debug irf_elliptical method
//#define G_DEBUG_NROI_RADIAL                  //!< Debug npred_radial method
//#define G_DEBUG_NROI_ELLIPTICAL          //!< Debug npred_elliptical method
//#define G_DEBUG_NROI_DIFFUSE                //!< Debug npred_diffuse method
//#define G_DEBUG_PRINT_AEFF                   //!< Debug print() Aeff method
//#define G_DEBUG_PRINT_PSF                     //!< Debug print() Psf method
//#define G_DEBUG_PSF_DUMMY_SIGMA           //!< Debug psf_dummy_sigma method

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
GCTAResponseIrf::GCTAResponseIrf(void) : GCTAResponse()
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
 * Constructs CTA response by making a deep copy of an existing object.
 **************************************************************************/
GCTAResponseIrf::GCTAResponseIrf(const GCTAResponseIrf& rsp) : GCTAResponse(rsp)
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
GCTAResponseIrf::GCTAResponseIrf(const GXmlElement& xml) : GCTAResponse()
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
 * @param[in] rspname Response file name.
 * @param[in] caldb Calibration database.
 *
 * Create instance of CTA response by specifying the response name and the
 * calibration database. The response name can be either a response identifier
 * or a filename (see GCTAResponseIrf::load for more information).
 ***************************************************************************/
GCTAResponseIrf::GCTAResponseIrf(const std::string& rspname,
                                 const GCaldb& caldb) : GCTAResponse()
{
    // Initialise members
    init_members();

    // Set calibration database
    m_caldb = caldb;

    // Load IRF
    load(rspname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of CTA response object.
 ***************************************************************************/
GCTAResponseIrf::~GCTAResponseIrf(void)
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
GCTAResponseIrf& GCTAResponseIrf::operator=(const GCTAResponseIrf& rsp)
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
void GCTAResponseIrf::clear(void)
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
GCTAResponseIrf* GCTAResponseIrf::clone(void) const
{
    return new GCTAResponseIrf(*this);
}


/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 *
 * @todo Set polar angle phi of photon in camera system
 ***************************************************************************/
double GCTAResponseIrf::irf(const GEvent&       event,
                            const GPhoton&      photon,
                            const GObservation& obs) const
{
    // Retrieve CTA pointing and instrument direction
    const GCTAPointing& pnt = gammalib::cta_pnt(G_IRF, obs);
    const GCTAInstDir&  dir = gammalib::cta_dir(G_IRF, event);

    // Get event attributes
    const GSkyDir& obsDir = dir.dir();
    const GEnergy& obsEng = event.energy();

    // Get photon attributes
    const GSkyDir& srcDir  = photon.dir();
    const GEnergy& srcEng  = photon.energy();
    const GTime&   srcTime = photon.time();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta = pnt.dir().dist(srcDir);
    double phi   = 0.0; //TODO: Implement Phi dependence

    // Get log10(E/TeV) of true photon energy.
    double srcLogEng = srcEng.log10TeV();

    // Determine angular separation between true and measured photon
    // direction in radians
    double delta = obsDir.dist(srcDir);

    // Get maximum angular separation for which PSF is significant
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

    // Initialise IRF value
    double irf = 0.0;

    // Compute only if we're sufficiently close to PSF
    if (delta <= delta_max) {

        // Get effective area component
        irf = aeff(theta, phi, zenith, azimuth, srcLogEng);

        // Multiply-in PSF
        if (irf > 0.0) {

            // Get PSF component
            irf *= psf(delta, theta, phi, zenith, azimuth, srcLogEng);

            // Multiply-in energy dispersion
            if (use_edisp() && irf > 0.0) {

                // Multiply-in energy dispersion
                irf *= edisp(obsEng, srcEng, theta, phi, zenith, azimuth);

            } // endif: energy dispersion was available and PSF was non-zero

            // Apply deadtime correction
            irf *= obs.deadc(srcTime);

        } // endif: Aeff was non-zero

    } // endif: we were sufficiently close to PSF

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCTAResponseIrf::irf:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ", theta=" << theta;
        std::cout << ", phi=" << phi << ")";
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
 * @return Event probability.
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
 ***************************************************************************/
double GCTAResponseIrf::nroi(const GModelSky&    model,
                             const GEnergy&      obsEng,
                             const GTime&        obsTime,
                             const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    static const int iter = 6;

    // Initialise Nroi value
    double nroi = 0.0;

    // No time dispersion supported
    const GTime& srcTime = obsTime;

    // If energy dispersion is requested then integrate over the relevant
    // true photon energies ...
    if (use_edisp()) {

        // Retrieve true energy boundaries
        GEbounds ebounds = edisp()->etrue_bounds(obsEng);

        // Loop over all boundaries
        for (int i = 0; i < ebounds.size(); ++i) {

            // Get true energy boundaries in MeV
            double etrue_min = ebounds.emin(i).MeV();
            double etrue_max = ebounds.emax(i).MeV();

            // Continue only if valid
            if (etrue_max > etrue_min) {

                // Setup integration function
                cta_nroi_kern integrand(this, &obs, &model, srcTime, obsEng, obsTime);
                GIntegral integral(&integrand);

                // Set fixed number of iterations
                integral.fixed_iter(iter);

                // Do Romberg integration
                nroi += integral.romberg(etrue_min, etrue_max);

            } // endif: interval was valid

        } // endfor: looped over energy intervals

    } // endif: energy dispersion requested

    // ... otherwise evaluate
    else {

        // No energy dispersion
        const GEnergy& srcEng = obsEng;

        // Compute response components
        double nroi_spatial  = this->nroi(model, srcEng, srcTime, obsEng, obsTime, obs);
        double nroi_spectral = model.spectral()->eval(srcEng, srcTime);
        double nroi_temporal = model.temporal()->eval(srcTime);

        // Compute response
        nroi = nroi_spatial * nroi_spectral * nroi_temporal;

    } // endelse: no energy dispersion requested

    // If required, apply instrument specific model scaling
    if (model.has_scales()) {
        nroi *= model.scale(obs.instrument()).value();
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
        std::cout << "*** ERROR: GCTAResponseIrf::nroi:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (nroi=" << nroi;
        std::cout << ", obsEng=" << obsEng;
        std::cout << ", obsTime=" << obsTime;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return response value
    return nroi;
}


/***********************************************************************//**
 * @brief Return true energy boundaries for a specific observed energy
 *
 * @param[in] obsEnergy Observed Energy.
 * @return True energy boundaries for given observed energy.
 ***************************************************************************/
GEbounds GCTAResponseIrf::ebounds(const GEnergy& obsEnergy) const
{
    // Initialise an empty boundary object
    GEbounds ebounds;

    // If energy dispersion is available then set the energy boundaries
    if (edisp() != NULL) {
        ebounds = edisp()->etrue_bounds(obsEnergy);
    }

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Simulate event from photon
 *
 * @param[in] area Simulation surface area.
 * @param[in] photon Photon.
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return Simulated event.
 *
 * Simulates a CTA event using the response function from an incident photon.
 * If the event is not detected a NULL pointer is returned.
 *
 * The method also applies a deadtime correction using a Monte Carlo process,
 * taking into account temporal deadtime variations. For this purpose, the
 * method makes use of the time dependent GObservation::deadc method.
 *
 * @todo Set polar angle phi of photon in camera system
 * @todo Implement energy dispersion
 ***************************************************************************/
GCTAEventAtom* GCTAResponseIrf::mc(const double&       area,
                                   const GPhoton&      photon,
                                   const GObservation& obs,
                                   GRan&               ran) const
{
    // Initialise event
    GCTAEventAtom* event = NULL;

    // Retrieve CTA pointing
    const GCTAPointing& pnt = gammalib::cta_pnt(G_MC, obs);

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta = pnt.dir().dist(photon.dir());
    double phi   = 0.0;  //TODO Implement Phi dependence

    // Compute effective area for photon
    double srcLogEng      = photon.energy().log10TeV();
    double effective_area = aeff(theta, phi, zenith, azimuth, srcLogEng);

    // Compute limiting value
    double ulimite = effective_area / area;

    // Warning if ulimite is larger than one
    if (ulimite > 1.0) {
        std::string msg = "Effective area "+
                          gammalib::str(effective_area)+
                          " cm2 is larger than simulation surface area "+
                          gammalib::str(area)+
                          " cm2 for photon energy "+
                          gammalib::str(photon.energy().TeV())+
                          " TeV. Simulations are inaccurate.";
        gammalib::warning(G_MC, msg);
    }

    // Continue only if event is detected
    if ((ulimite > 0.0) && (ran.uniform() <= ulimite)) {

        // Apply deadtime correction
        double deadc = obs.deadc(photon.time());
        if (deadc >= 1.0 || ran.uniform() <= deadc) {

            // Simulate PSF and energy dispersion. Skip event if exception
            // occurs
           try {
                // Simulate offset from photon arrival direction.
                double delta = psf()->mc(ran, srcLogEng, theta, phi, zenith, azimuth) *
                               gammalib::rad2deg;
                double alpha = 360.0 * ran.uniform();

                // Rotate sky direction by offset
                GSkyDir sky_dir = photon.dir();
                sky_dir.rotate_deg(alpha, delta);

                // Set measured photon arrival direction in instrument direction
                GCTAInstDir inst_dir = pnt.instdir(sky_dir);

                // Set measured photon energy
                GEnergy energy = photon.energy();
                if (use_edisp()) {
                    energy = edisp()->mc(ran, photon.energy(), theta, phi,
                                         zenith, azimuth);
                }

                // Allocate event
                event = new GCTAEventAtom;

                // Set event attributes
                event->dir(inst_dir);
                event->energy(energy);
                event->time(photon.time());

            } // endtry: PSF and energy dispersion simulation successful

            // ... otherwise catch exception
            catch (GException::invalid_return_value) {
                ;
            }

        } // endif: detector was alive

    } // endif: event was detected

    // Return event
    return event;
}


/***********************************************************************//**
 * @brief Read response from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads information for a CTA observation from an XML element. The
 * calibration database and response name can be specified using
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="Calibration" database="..." response="..."/>
 *     </observation>
 *
 * If even more control is required over the response, individual file names
 * can be specified using
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *       <parameter name="Background"          file="..."/>
 *     </observation>
 *
 ***************************************************************************/
void GCTAResponseIrf::read(const GXmlElement& xml)
{
    // First check for "Calibration" parameter
    if (gammalib::xml_has_par(xml, "Calibration")) {

        // Get parameter
        const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Calibration");

        // Read database and response
        std::string xml_caldb   = gammalib::strip_whitespace(par->attribute("database"));
        std::string xml_rspname = gammalib::strip_whitespace(par->attribute("response"));

        // Set calibration database
        GCaldb caldb("cta", xml_caldb);
        this->caldb(caldb);

        // Load response
        this->load(xml_rspname);

        // If there is a sigma attribute then set sigma values of Aeff and
        // Background
        if (par->has_attribute("sigma")) {

            // Get sigma value
            double sigma = gammalib::todouble(par->attribute("sigma"));

            // If we have an effective area performance table then set sigma
            // value
            GCTAAeffPerfTable* perf = const_cast<GCTAAeffPerfTable*>(dynamic_cast<const GCTAAeffPerfTable*>(aeff()));
            if (perf != NULL) {
                perf->sigma(sigma);
            }

            // If we have a background performance table then set sigma value
            GCTABackgroundPerfTable* bgm = const_cast<GCTABackgroundPerfTable*>(dynamic_cast<const GCTABackgroundPerfTable*>(background()));
            if (bgm != NULL) {
                bgm->sigma(sigma);
            }

        } // endif: sigma attribute specified

        // Store database and response names for later writing into the XML
        // file (we do this now since the load() method clears all
        // GCTAResponseIrf data members)
        //m_xml_caldb   = xml_caldb;
        //m_xml_rspname = xml_rspname;

    } // endif: "Calibration" parameter found

    // ... otherwise check for components
    else {

        // Handle effective area
        if (gammalib::xml_has_par(xml, "EffectiveArea")) {

            // Get parameter
            const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "EffectiveArea");

            // Get filename
            std::string filename = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load effective area
            if (!filename.empty()) {

                // Load effective area
                load_aeff(filename);

                // Optional attributes
                double thetacut = 0.0;
                double scale    = 1.0;
                double sigma    = 3.0;

                // Optionally extract thetacut (0.0 if no thetacut)
                if (par->has_attribute("thetacut")) {
                    thetacut = gammalib::todouble(par->attribute("thetacut"));
                }

                // Optionally extract scale factor (1.0 if no scale)
                if (par->has_attribute("scale")) {
                    scale = gammalib::todouble(par->attribute("scale"));
                }

                // Optionally extract sigma (3.0 if no sigma)
                if (par->has_attribute("sigma")) {
                    sigma = gammalib::todouble(par->attribute("sigma"));
                }

                // If we have an ARF then set attributes
                GCTAAeffArf* arf = const_cast<GCTAAeffArf*>(dynamic_cast<const GCTAAeffArf*>(aeff()));
                if (arf != NULL) {
                    arf->thetacut(thetacut);
                    arf->scale(scale);
                    arf->sigma(sigma);
                }

                // If we have a performance table then set attributes
                GCTAAeffPerfTable* perf = const_cast<GCTAAeffPerfTable*>(dynamic_cast<const GCTAAeffPerfTable*>(aeff()));
                if (perf != NULL) {
                    perf->sigma(sigma);
                }

            } // endif: effective area filename was valid

        } // endif: handled effective area

        // Handle PSF
        if (gammalib::xml_has_par(xml, "PointSpreadFunction")) {

            // Get parameter
            const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "PointSpreadFunction");

            // Get filename
            std::string filename = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load point spread function
            if (!filename.empty()) {
                load_psf(filename);
            }

        } // endif: handled PSF

        // Handle energy dispersion
        if (gammalib::xml_has_par(xml, "EnergyDispersion")) {

            // Get parameter
            const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "EnergyDispersion");

            // Get filename
            std::string filename = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load energy dispersion
            if (!filename.empty()) {
                load_edisp(filename);
            }

        } // endif: handled energy dispersion

        // Handle Background
        if (gammalib::xml_has_par(xml, "Background")) {

            // Get parameter
            const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Background");

            // Get filename
            std::string filename = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load background model
            if (!filename.empty()) {
                load_background(filename);
            }

            // Optional attributes
            double sigma = 3.0;

            // Optionally extract sigma (3.0 if no sigma)
            if (par->has_attribute("sigma")) {
                sigma = gammalib::todouble(par->attribute("sigma"));
            }

            // If we have a performance table then set attributes
            GCTABackgroundPerfTable* bgm = const_cast<GCTABackgroundPerfTable*>(dynamic_cast<const GCTABackgroundPerfTable*>(background()));
            if (bgm != NULL) {
                bgm->sigma(sigma);
            }

        }

    } // endelse: handled components

    // If we have an ARF then remove thetacut if necessary
    GCTAAeffArf* arf = const_cast<GCTAAeffArf*>(dynamic_cast<const GCTAAeffArf*>(aeff()));
    if (arf != NULL) {
        if (arf->thetacut() > 0.0) {
            arf->remove_thetacut(*this);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response information into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes information for a CTA response into an XML element. If the
 * calibration database and response name had been specified, the following
 * output is written
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="Calibration" database="..." response="..."/>
 *     </observation>
 *
 * If even more control was required over the response and individual file
 * names were specified, the following output is written
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *       <parameter name="Background"          file="..."/>
 *     </observation>
 ***************************************************************************/
void GCTAResponseIrf::write(GXmlElement& xml) const
{
    // Determine number of existing parameter nodes in XML element
    //int npars = xml.elements("parameter");

    // If we have a calibration database and response name, then set
    // the information ...
    if (!m_xml_caldb.empty() || !m_xml_rspname.empty()) {
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "Calibration");
        par->attribute("database", m_xml_caldb);
        par->attribute("response", m_xml_rspname);
    }

    // ... otherwise add response components if they exist
    else {

        // Add effective area if it exists
        if (aeff() != NULL) {

            // Get effective area filename
            GFilename filename = aeff()->filename();

            // Continue only if filename is not empty
            if (!(filename.is_empty())) {

                // Get pointer to effective area
                GXmlElement* par =
                    gammalib::xml_need_par(G_WRITE, xml, "EffectiveArea");

                // Initialise attributes
                double thetacut = 0.0;
                double scale    = 1.0;
                double sigma    = 0.0;

                // Get optional ARF attributes
                const GCTAAeffArf* arf =
                      dynamic_cast<const GCTAAeffArf*>(aeff());
                if (arf != NULL) {
                    thetacut = arf->thetacut();
                    scale    = arf->scale();
                    sigma    = arf->sigma();
                }

                // Get optional performance table attributes
                const GCTAAeffPerfTable* perf =
                      dynamic_cast<const GCTAAeffPerfTable*>(aeff());
                if (perf != NULL) {
                    sigma = perf->sigma();
                }

                // Set attributes
                par->attribute("file", filename);
                if (thetacut > 0.0) {
                    par->attribute("thetacut", gammalib::str(thetacut));
                }
                if (scale != 1.0) {
                    par->attribute("scale", gammalib::str(scale));
                }
                if (sigma > 0.0) {
                    par->attribute("sigma", gammalib::str(sigma));
                }

            }
        }

        // Add PSF if it exists
        if (psf() != NULL) {
            if (!(psf()->filename().is_empty())) {
                GXmlElement* par =
                    gammalib::xml_need_par(G_WRITE, xml, "PointSpreadFunction");
                par->attribute("file", psf()->filename());
            }
        }

        // Add Edisp if it exists
        if (edisp() != NULL) {
            if (!(edisp()->filename().is_empty())) {
                GXmlElement* par =
                    gammalib::xml_need_par(G_WRITE, xml, "EnergyDispersion");
                par->attribute("file", edisp()->filename());
            }
        }

        // Add background if it exists
        if (background() != NULL) {
            if (!(background()->filename().is_empty())) {
                GXmlElement* par =
                    gammalib::xml_need_par(G_WRITE, xml, "Background");
                par->attribute("file", background()->filename());
            }
        }

    } // endelse: response components added

    // Return
    return;
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
 *
 * The method sets the calibration database and response names for writing
 * of the response information into the XML file.
 ***************************************************************************/
void GCTAResponseIrf::load(const std::string& rspname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // First attempt reading the response using the GCaldb interface
    std::string expr      = "NAME("+rspname+")";
    GFilename   aeffname  = m_caldb.filename("","","EFF_AREA","","",expr);
    GFilename   psfname   = m_caldb.filename("","","RPSF","","",expr);
    GFilename   edispname = m_caldb.filename("","","EDISP","","",expr);
    GFilename   bgdname   = m_caldb.filename("","","BKG","","",expr);

    // Kluge: if the background filename is empty it may be because we have
    // and old response file that used "BGD" as the name of the background
    // extension. So we try to get here the old name
    if (bgdname.is_empty()) {
        bgdname = m_caldb.filename("","","BGD","","",expr);
    }

    // Signal usage of calibration database
    bool use_caldb = true;

    // If filenames are empty then build filenames from CALDB root path and
    // response name
    if (aeffname.is_empty()) {
        aeffname  = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
        use_caldb = false;
    }
    if (psfname.is_empty()) {
        psfname   = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
        use_caldb = false;
    }
    if (edispname.is_empty()) {
        edispname = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
        use_caldb = false;
    }
    if (bgdname.is_empty()) {
        bgdname   = irf_filename(gammalib::filepath(m_caldb.rootdir(), rspname));
        use_caldb = false;
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

    // If the calibration database has been used then store database and
    // response names for later writing into the XML file
    if (use_caldb) {
        m_xml_caldb   = caldb.instrument();
        m_xml_rspname = rspname;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load effective area
 *
 * @param[in] filename Effective area filename.
 *
 * @exception GException::file_error
 *            File or extension not found.
 *
 * Loads the effective area from a response file.
 *
 * If the file is a FITS file, the method will either use the extension name
 * specified with the filename, or if no extension name is given, search for
 * an `EFFECTIVE AREA` or `SPECRESP` extension in the file and open the
 * corresponding FITS table. If a column named `SPECRESP` is found in the
 * table, a GCTAAeffArf object will be allocated. Otherwise, a GCTAAeff2D
 * object will be allocated. In both cases, the method will extract the
 * optional `LO_THRES` and `HI_THRES` safe energy thresholds from the FITS
 * file header.
 *
 * If the file is not a FITS file, it will be interpreted as a
 * GCTAAeffPerfTable performance table.
 ***************************************************************************/
void GCTAResponseIrf::load_aeff(const GFilename& filename)
{
    // Free any existing effective area instance
    if (m_aeff != NULL) delete m_aeff;
    m_aeff = NULL;

    // Check for existence of file
    if (!filename.exists()) {
        std::string msg = "File \""+filename+"\" not found. Please specify "
                          "a valid effective area response file.";
        throw GException::file_error(G_LOAD_AEFF, msg);
    }

    // If file is a FITS file ...
    if (filename.is_fits()) {

        // Open FITS file
        GFits fits(filename);

        // Get the extension name. If an extension name has been specified
        // then use this name, otherwise preset the extension name according
        // to the GADF specifications. If no GADF compliant name was found,
        // search either for "EFFECTIVE AREA" or the "SPECRESP" extension.
        std::string extname = gammalib::gadf_hduclas4(fits, "AEFF_2D");
        if (filename.has_extname()) {
            extname = filename.extname();
        }
        if (extname.empty()) {
            if (fits.contains(gammalib::extname_cta_aeff2d)) {
                extname = gammalib::extname_cta_aeff2d;
            }
            else if (fits.contains(gammalib::extname_arf)) {
                extname = gammalib::extname_arf;
            }
        }

        // Continue only if extension name is not empty
        if (!extname.empty()) {

            // Get FITS table
            const GFitsTable& table = *fits.table(extname);

            // Read safe energy thresholds if available
            if (table.has_card("LO_THRES")) {
                m_lo_safe_thres = table.real("LO_THRES");
            }
            if (table.has_card("HI_THRES")) {
                m_hi_safe_thres = table.real("HI_THRES");
            }

            // Check for specific table column
            if (table.contains(gammalib::extname_arf)) {

                // Close file
                fits.close();

                // Allocate Aeff from file
                m_aeff = new GCTAAeffArf(filename);

            } // endif: load as GCTAAeffArf

            else {

                // Close file
                fits.close();

                // Allocate Aeff from file
                m_aeff = new GCTAAeff2D(filename);

            } // endelse: load as GCTAAeff2D

        } // endif: extension name is not empty

        // Signal that no extension was found
        else {
            std::string msg = "FITS file \""+filename+"\" does not "
                              "contain a valid effective area table. "
                              "Please specify a valid effective area "
                              "response file.";
            throw GException::file_error(G_LOAD_AEFF, msg);
        }

    } // endif: file was FITS file

    // ... else try handling an ASCII performance table
    else {

        // Allocate a performance table
        m_aeff = new GCTAAeffPerfTable(filename);

    } // endif: load as GCTAAeffPerfTable

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load CTA PSF vector
 *
 * @param[in] filename FITS file name.
 *
 * @exception GException::file_error
 *            File or extension not found.
 *
 * Loads the point spead function from a response file.
 *
 * If the file is a FITS file, the method will either use the extension name
 * specified with the filename, or if no extension name is given, search for
 * one of the following extension names
 *
 *       POINT SPREAD FUNCTION
 *       PSF
 *       PSF_2D_TABLE
 *
 * in the file and open the corresponding FITS table. If columns named `GAMMA`
 * and `SIGMA` are found in the table, a GCTAPsfKing object will be allocated.
 *
 * If columns named `SCALE`, `SIGMA_1`, `AMPL_2`, `SIGMA_2`, `AMPL_3` and
 * `SIGMA_3` are found in the table, a GCTAPsf2D object will be allocated.
 *
 * If columns named `RAD_LO`, `RAD_HI`, and `RPSF` are found in the table, a
 * GCTAPsfTable object will be allocated.
 *
 * Otherwise, a CTAPsfVector object will be allocated.
 *
 * If the file is not a FITS file, it will be interpreted as a
 * GCTAPsfPerfTable performance table.
 ***************************************************************************/
void GCTAResponseIrf::load_psf(const GFilename& filename)
{
    // Free any existing point spread function instance
    if (m_psf != NULL) delete m_psf;
    m_psf = NULL;

    // Check for existence of file
    if (!filename.exists()) {
        std::string msg = "File \""+filename+"\" not found. Please specify "
                          "a valid point spread function response file.";
        throw GException::file_error(G_LOAD_PSF, msg);
    }

    // If file is a FITS file ...
    if (filename.is_fits()) {

        // Open FITS file
        GFits fits(filename);

        // Get the extension name. If an extension name has been specified
        // then use this name, otherwise preset the extension name according
        // to the GADF specifications. If no GADF compliant name was found,
        // search either for "POINT SPREAD FUNCTION" or the "PSF" extension.
        std::string extname = gammalib::gadf_hduclas4(fits, "PSF_TABLE");
        if (filename.has_extname()) {
            extname = filename.extname();
        }
        if (extname.empty()) {
            if (fits.contains("POINT SPREAD FUNCTION")) {
                extname = "POINT SPREAD FUNCTION";
            }
            else if (fits.contains("PSF")) {
                extname = "PSF";
            }
        }

        // Continue only if extension name is not empty
        if (!extname.empty()) {

            // Get FITS table
            const GFitsTable& table = *fits.table(extname);

            // Check for King profile specific table columns
            if (table.contains("GAMMA") && table.contains("SIGMA")) {

                // Close FITS file
                fits.close();

                // Allocate King profile PSF
                m_psf = new GCTAPsfKing(filename);

            }

            // ... otherwise check for Gaussian profile specific table
            // columns
            else if (table.contains("SCALE")  && table.contains("SIGMA_1") &&
                     table.contains("AMPL_2") && table.contains("SIGMA_2") &&
                     table.contains("AMPL_3") && table.contains("SIGMA_3")) {

                // Close FITS file
                fits.close();

                // Allocate Gaussian profile PSF
                m_psf = new GCTAPsf2D(filename);

            }

            // ... otherwise check for PSF table specific table columns
            else if (table.contains("RAD_LO") && table.contains("RAD_HI") &&
                     table.contains("RPSF")) {

                // Close FITS file
                fits.close();

                // Allocate PSF table
                m_psf = new GCTAPsfTable(filename);

            }

            // ... otherwise try opening as vector PSF
            else {

                // Close FITS file
                fits.close();

                // Allocate vector PSF
                m_psf = new GCTAPsfVector(filename);

            }

        } // endif: extension name is not empty

        // Signal that no extension was found
        else {
            std::string msg = "FITS file \""+filename+"\" does not "
                              "contain a valid point spread function table. "
                              "Please specify a valid point spread function "
                              "response file.";
            throw GException::file_error(G_LOAD_PSF, msg);
        }

    } // endif: file was FITS file

    // ... otherwise load file as a performance table
    else {

        // Allocate a performance table
        m_psf = new GCTAPsfPerfTable(filename);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy dispersion information
 *
 * @param[in] filename Energy dispersion file name.
 *
 * @exception GException::file_error
 *            File or extension not found.
 *
 * Loads the energy dispersion from a response file.
 *
 * If the file is a FITS file, the method will either use the extension name
 * specified with the filename, or if no extension name is given, search for
 * an `ENERGY DISPERSION` or `MATRIX` extension in the file and open the
 * corresponding FITS table. If columns named "MIGRA_LO" and "MIGRA_HI" are
 * found in the table, a GCTAEdisp2D object will be allocated. Otherwise, a
 * GCTAEdispRmf object will be allocated.
 *
 * If the file is not a FITS file, it will be interpreted as a
 * GCTAEdispPerfTable performance table.
 ***************************************************************************/
void GCTAResponseIrf::load_edisp(const GFilename& filename)
{
    // Free any existing energy dispersion instance
    if (m_edisp != NULL) delete m_edisp;
    m_edisp = NULL;

    // Check for existence of file
    if (!filename.exists()) {
        std::string msg = "File \""+filename+"\" not found. Please specify "
                          "a valid energy dispersion response file.";
        throw GException::file_error(G_LOAD_EDISP, msg);
    }

    // If file is a FITS file ...
    if (filename.is_fits()) {

        // Open FITS file
        GFits fits(filename);

        // Get the extension name. If an extension name has been specified
        // then use this name, otherwise preset the extension name according
        // to the GADF specifications. If no GADF compliant name was found,
        // search either for "ENERGY DISPERSION" or the "MATRIX" extension.
        std::string extname = gammalib::gadf_hduclas4(fits, "EDISP_2D");
        if (filename.has_extname()) {
            extname = filename.extname();
        }
        if (extname.empty()) {
            if (fits.contains(gammalib::extname_cta_edisp2d)) {
                extname = gammalib::extname_cta_edisp2d;
            }
            else if (fits.contains(gammalib::extname_rmf)) {
                extname = gammalib::extname_rmf;
            }
        }

        // Continue only if extension name is not empty
        if (!extname.empty()) {

            // Get FITS table
            const GFitsTable& table = *fits.table(extname);

            // Check for 2D migration matrix
            if (table.contains("MIGRA_LO") && table.contains("MIGRA_HI")) {

                // Close FITS file
                fits.close();

                // Allocate 2D migration matrix
                m_edisp = new GCTAEdisp2D(filename);

            }

            // ... otherwise allocate RMF
            else {

                // Close FITS file
                fits.close();

                // Allocate Gaussian profile PSF
                m_edisp = new GCTAEdispRmf(filename);

            }

        } // endif: extension name is not empty

        // Signal that no extension was found
        else {
            std::string msg = "FITS file \""+filename+"\" does not "
                              "contain a valid energy dispersion table. "
                              "Please specify a valid energy dispersion "
                              "response file.";
            throw GException::file_error(G_LOAD_EDISP, msg);
        }

    } // endif: file was FITS file

    // ... otherwise load file as a performance table
    else {

        // Allocate a performance table
        m_edisp = new GCTAEdispPerfTable(filename);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load background model
 *
 * @param[in] filename Background model file name.
 *
 * @exception GException::file_error
 *            File not found.
 *
 * Loads the background model from a response file.
 *
 * If the file is a FITS file the method allocates a GCTABackground3D
 * response and loads the information from the FITS file.
 *
 * If the file is not a FITS file it will be interpreted as a
 * GCTABackgroundPerfTable performance table.
 ***************************************************************************/
void GCTAResponseIrf::load_background(const GFilename& filename)
{
    // Free any existing background model instance
    if (m_background != NULL) delete m_background;
    m_background = NULL;

    // Check for existence of file
    if (!filename.exists()) {
        std::string msg = "File \""+filename+"\" not found. Please specify "
                          "a valid background response file.";
        throw GException::file_error(G_LOAD_BACKGROUND, msg);
    }

    // If file is a FITS file than load background as 3D background
    if (filename.is_fits()) {
        m_background = new GCTABackground3D(filename);
    }

    // ... otherwise load background as performance table
    else {
        m_background = new GCTABackgroundPerfTable(filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set offset angle dependence (degrees)
 *
 * @param[in] sigma Offset angle dependence value (degrees).
 *
 * Set the offset angle dependence for 1D effective area functions. The
 * method set the sigma value in case that the effective area function
 * is of type GCTAAeffArf or GCTAAeffPerfTable. Otherwise, nothing will
 * be done.
 ***************************************************************************/
void GCTAResponseIrf::offset_sigma(const double& sigma)
{
    // If effective area is an ARF then set offset angle
    GCTAAeffArf* arf = dynamic_cast<GCTAAeffArf*>(m_aeff);
    if (arf != NULL) {
        arf->sigma(sigma);
    }

    // If effective area is a performance table then set offset angle
    GCTAAeffPerfTable* prf = dynamic_cast<GCTAAeffPerfTable*>(m_aeff);
    if (prf != NULL) {
        prf->sigma(sigma);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return offset angle dependence (degrees)
 *
 * @return Offset angle dependence value (degrees).
 *
 * Return the offset angle dependence for 1D effective area functions. The
 * method returns the sigma value in case that the effective area function
 * is of type GCTAAeffArf or GCTAAeffPerfTable. Otherwise, 0.0 will be
 * returned.
 ***************************************************************************/
double GCTAResponseIrf::offset_sigma(void) const
{
    // Initialise value
    double sigma = 0.0;

    // If effective area is an ARF then get offset angle
    GCTAAeffArf* arf = dynamic_cast<GCTAAeffArf*>(m_aeff);
    if (arf != NULL) {
        sigma = arf->sigma();
    }

    // If effective area is a performance table then get offset angle
    GCTAAeffPerfTable* prf = dynamic_cast<GCTAAeffPerfTable*>(m_aeff);
    if (prf != NULL) {
        sigma = prf->sigma();
    }

    // Return sigma
    return sigma;
}


/***********************************************************************//**
 * @brief Print CTA response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing CTA response information.
 ***************************************************************************/
std::string GCTAResponseIrf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAResponseIrf ===");

        // Append response information
        result.append("\n"+gammalib::parformat("Caldb mission")+m_caldb.mission());
        result.append("\n"+gammalib::parformat("Caldb instrument")+m_caldb.instrument());
        result.append("\n"+gammalib::parformat("Response name")+m_rspname);
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

        // Append safe energy threshold information
        result.append("\n"+gammalib::parformat("Safe energy range"));
        if (m_lo_safe_thres > 0.0 && m_hi_safe_thres) {
            result.append(gammalib::str(m_lo_safe_thres));
            result.append(" - ");
            result.append(gammalib::str(m_hi_safe_thres));
            result.append(" TeV");
        }
        else if (m_lo_safe_thres > 0.0) {
            result.append("> ");
            result.append(gammalib::str(m_lo_safe_thres));
            result.append(" TeV");
        }
        else if (m_hi_safe_thres > 0.0) {
            result.append("< ");
            result.append(gammalib::str(m_hi_safe_thres));
            result.append(" TeV");
        }
        else {
            result.append("undefined");
        }

        // Get reduced chatter level
        GChatter reduced_chatter = gammalib::reduce(chatter);

        // Append detailed information
        if (chatter >= NORMAL) {

            // Append calibration database
            result.append("\n"+m_caldb.print(reduced_chatter));

            // Append effective area information
            if (m_aeff != NULL) {
                result.append("\n"+m_aeff->print(reduced_chatter));
            }

            // Append point spread function information
            if (m_psf != NULL) {
                result.append("\n"+m_psf->print(reduced_chatter));
            }

            // Append energy dispersion information
            if (m_edisp != NULL) {
                result.append("\n"+m_edisp->print(reduced_chatter));
            }

            // Append background information
            if (m_background != NULL) {
                result.append("\n"+m_background->print(reduced_chatter));
            }

            // Append cache information
            result.append("\n"+m_irf_cache.print(reduced_chatter));
            result.append("\n"+m_nroi_cache.print(reduced_chatter));

        } // endif: appended detailed information

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                    Low-level CTA response methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return effective area (in units of cm2)
 *
 * @param[in] theta Radial offset angle of photon in camera (radians).
 * @param[in] phi Polar angle of photon in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 * @return Effective area in units fo cm2.
 *
 * @exception GException::invalid_value
 *            No effective area information found.
 *
 * Returns the effective area as function of the true photon position in the
 * camera system and the telescope pointing direction in the Earth system.
 ***************************************************************************/
double GCTAResponseIrf::aeff(const double& theta,
                             const double& phi,
                             const double& zenith,
                             const double& azimuth,
                             const double& srcLogEng) const
{
    // Throw an exception if instrument response is not defined
    if (m_aeff == NULL) {
        std::string msg = "No effective area information found in response.\n"
                          "Please make sure that the instrument response is"
                          " properly defined.";
        throw GException::invalid_value(G_AEFF, msg);
    }

    // Get effective area
    double aeff = (*m_aeff)(srcLogEng, theta, phi, zenith, azimuth);

    // Return effective area
    return aeff;
}


/***********************************************************************//**
 * @brief Return point spread function (in units of sr^-1)
 *
 * @param[in] delta Angular separation between true and measured photon
 *            directions (radians).
 * @param[in] theta Radial offset angle of photon in camera (radians).
 * @param[in] phi Polar angle of photon in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * @exception GException::invalid_value
 *            No point spread function information found.
 *
 * Returns the point spread function for a given offset angle as function
 * of the true photon position in the camera system and the telescope
 * pointing direction in the Earth system.
 ***************************************************************************/
double GCTAResponseIrf::psf(const double& delta,
                            const double& theta,
                            const double& phi,
                            const double& zenith,
                            const double& azimuth,
                            const double& srcLogEng) const
{
    // Throw an exception if instrument response is not defined
    if (m_psf == NULL) {
        std::string msg = "No point spread function information found in"
                          " response.\n"
                          "Please make sure that the instrument response is"
                          " properly defined.";
        throw GException::invalid_value(G_PSF, msg);
    }

    // Compute PSF
    double psf = (*m_psf)(delta, srcLogEng, theta, phi, zenith, azimuth);

    // Return PSF
    return psf;
}


/***********************************************************************//**
 * @brief Return maximum angular separation (in radians)
 *
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * @exception GException::invalid_value
 *            No point spread function information found.
 *
 * This method returns the maximum angular separation between true and
 * measured photon directions for which the PSF is non zero as function
 * of the true photon position in the camera system and the telescope
 * pointing direction in the Earth system.
 ***************************************************************************/
double GCTAResponseIrf::psf_delta_max(const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth,
                                      const double& srcLogEng) const
{
    // Throw an exception if instrument response is not defined
    if (m_psf == NULL) {
        std::string msg = "No point spread function information found in"
                          " response.\n"
                          "Please make sure that the instrument response is"
                          " properly defined.";
        throw GException::invalid_value(G_PSF_DELTA_MAX, msg);
    }

    // Compute PSF
    double delta_max = m_psf->delta_max(srcLogEng, theta, phi, zenith, azimuth);

    // Return PSF
    return delta_max;
}


/***********************************************************************//**
 * @brief Return energy dispersion (in units of MeV\f$^{-1}\f$)
 *
 * @param[in] ereco Reconstructed event energy.
 * @param[in] etrue True photon energy.
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @return Energy dispersion (MeV\f$^{-1}\f$).
 ***************************************************************************/
double GCTAResponseIrf::edisp(const GEnergy& ereco,
                              const GEnergy& etrue,
                              const double&  theta,
                              const double&  phi,
                              const double&  zenith,
                              const double&  azimuth) const
{
    // Compute energy dispersion
    double edisp = (*m_edisp)(ereco, etrue, theta, phi, zenith, azimuth);

    // Return energy dispersion
    return edisp;
}


/***********************************************************************//**
 * @brief Return spatial integral of sky model
 *
 * @param[in] model Sky Model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event arrival time.
 * @param[in] obs Observation.
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t'|E,t) = \int_{\rm ROI} P(p',E',t'|E,t) dp'
 * \f]
 *
 * of
 *
 * \f[
 *    P(p',E',t'|E,t) = \int
 *                      S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * over the Region of Interest (ROI) for a sky model \f$S(p,E,t)\f$ and the
 * response function \f$R(p',E',t'|p,E,t)\f$.
 ***************************************************************************/
double GCTAResponseIrf::nroi(const GModelSky&    model,
                             const GEnergy&      srcEng,
                             const GTime&        srcTime,
                             const GEnergy&      obsEng,
                             const GTime&        obsTime,
                             const GObservation& obs) const
{
    // Initialise response value
    double nroi = 0.0;

    // Set response value attributes
    std::string name = obs.id() + "::" + model.name();

    // Signal if spatial model has free parameters
    bool has_free_pars = model.spatial()->has_free_pars();

    // If the spatial model component has free parameters, or the response
    // cache should not be used, or the cache does not contain the requested
    // IRF value then compute the IRF value for the spatial model.
    if (has_free_pars     ||
        !m_use_nroi_cache ||
        !m_nroi_cache.contains(name, obsEng, srcEng, &nroi)) {

        // Select method depending on the spatial model type
        switch (model.spatial()->code()) {
            case GMODEL_SPATIAL_POINT_SOURCE:
                nroi = nroi_ptsrc(model, srcEng, srcTime, obsEng, obsTime, obs);
                break;
            case GMODEL_SPATIAL_RADIAL:
                nroi = nroi_radial(model, srcEng, srcTime, obsEng, obsTime, obs);
                break;
            case GMODEL_SPATIAL_ELLIPTICAL:
                nroi = nroi_elliptical(model, srcEng, srcTime, obsEng, obsTime, obs);
                break;
            case GMODEL_SPATIAL_DIFFUSE:
                nroi = nroi_diffuse(model, srcEng, srcTime, obsEng, obsTime, obs);
                break;
            case GMODEL_SPATIAL_COMPOSITE:
                nroi = nroi_composite(model, srcEng, srcTime, obsEng, obsTime, obs);
                break;
            default:
                break;
        }

    } // endif: computed spatial model

    // If the spatial model has no free parameters and the response cache
    // should be used then put the IRF value in the response cache.
    if (!has_free_pars && m_use_nroi_cache) {
        m_nroi_cache.set(name, obsEng, srcEng, nroi);
    }

    // Return response value
    return nroi;
}


/***********************************************************************//**
 * @brief Return spatial integral of Instrument Response Function
 *
 * @param[in] photon Photon.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event time.
 * @param[in] obs Observation.
 *
 * Computes the integral of the instrument response function over the Region
 * of Interest
 *
 * \f[
 *    R(E',t'|p,E,t) = \int_{\rm ROI} R(p',E',t'|p,E,t) dp'
 * \f]
 ***************************************************************************/
double GCTAResponseIrf::nirf(const GPhoton&      photon,
                             const GEnergy&      obsEng,
                             const GTime&        obsTime,
                             const GObservation& obs) const
{
    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = gammalib::cta_obs(G_NIRF, obs);
    const GCTARoi&         roi = gammalib::cta_event_list(G_NIRF, obs).roi();
    const GCTAPointing&    pnt = cta.pointing();

    // Get photon attributes
    const GSkyDir& srcDir  = photon.dir();
    const GEnergy& srcEng  = photon.energy();
    const GTime&   srcTime = photon.time();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta = pnt.dir().dist(srcDir);
    double phi   = 0.0; //TODO: Implement Phi dependence

    // Get log10(E/TeV) of true photon energy.
    double srcLogEng = srcEng.log10TeV();

    // Get effectve area components
    double nroi = aeff(theta, phi, zenith, azimuth, srcLogEng);

    // Multiply-in PSF
    if (nroi > 0.0) {

        // Get PSF
        nroi *= npsf(srcDir, srcLogEng, srcTime, pnt, roi);

        // Multiply-in energy dispersion
        if (use_edisp() && nroi > 0.0) {

            // Multiply-in energy dispersion
            nroi *= edisp(obsEng, srcEng, theta, phi, zenith, azimuth);

        } // endif: had energy dispersion

        // Apply deadtime correction
        nroi *= obs.deadc(srcTime);

    } // endif: had non-zero effective area

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
        std::cout << "*** ERROR: GCTAResponseIrf::nroi:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (nroi=" << nroi;
        std::cout << ", theta=" << theta;
        std::cout << ", phi=" << phi << ")";
        std::cout << std::endl;
    }
    #endif

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Return result of PSF integration over ROI.
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] pnt CTA pointing.
 * @param[in] roi CTA region of interest.
 *
 * This method integrates the PSF over the circular region of interest (ROI).
 * Integration is done in a polar coordinate system centred on the PSF since
 * the PSF is assumed to be azimuthally symmetric. The polar integration is
 * done using the method npsf_kern_rad_azsym() that computes analytically
 * the arclength that is comprised within the ROI.
 * 
 * Note that the integration is only performed when the PSF is spilling out
 * of the ROI border, otherwise the integral is simply 1. Numerical
 * integration is done using the standard Romberg method. The integration
 * boundaries are computed so that only the PSF section that falls in the ROI
 * is considered.
 *
 * @todo Enhance romberg() integration method for small integration regions
 *       (see comment about kluge below)
 * @todo Implement phi dependence in camera system
 ***************************************************************************/
double GCTAResponseIrf::npsf(const GSkyDir&      srcDir,
                             const double&       srcLogEng,
                             const GTime&        srcTime,
                             const GCTAPointing& pnt,
                             const GCTARoi&      roi) const
{
    // Set number of iterations for Romberg integration.
    static const int iter = 6;

    // Declare result
    double value = 0.0;

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Compute offset angle of source direction in camera system
    double theta = pnt.dir().dist(srcDir);

    // Compute azimuth angle of source direction in camera system
    double phi = 0.0; //TODO: Implement phi dependence

    // Extract relevant parameters from arguments
    double roi_radius       = roi.radius() * gammalib::deg2rad;
    double roi_psf_distance = roi.centre().dir().dist(srcDir);
    double rmax             = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

    // If PSF is fully enclosed by the ROI then skip the numerical
    // integration and assume that the integral is 1.0
    if (roi_psf_distance + rmax <= roi_radius) {
        value = 1.0;
    }

    // ... otherwise perform numerical integration
    else {

        // Compute minimum PSF integration radius
        double rmin = (roi_psf_distance > roi_radius) 
                      ? roi_psf_distance - roi_radius : 0.0;

        // Continue only if integration range is valid
        if (rmax > rmin) {

            // Setup integration kernel
            cta_npsf_kern_rad_azsym integrand(this,
                                              roi_radius,
                                              roi_psf_distance,
                                              srcLogEng,
                                              theta,
                                              phi,
                                              zenith,
                                              azimuth);

            // Setup integration
            GIntegral integral(&integrand);

            // Set fixed number of iterations
            integral.fixed_iter(iter);

            // Radially integrate PSF. In case that the radial integration
            // region is small, we do the integration using a simple
            // trapezoidal rule. This is a kluge to prevent convergence
            // problems in the romberg() method for small integration intervals.
            // Ideally, the romberg() method should be enhanced to handle this
            // case automatically. The kluge threshold was fixed manually!
            if (rmax-rmin < 1.0e-12) {
                value = integral.trapzd(rmin, rmax);
            }
            else {
                value = integral.romberg(rmin, rmax);
            }

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
                std::cout << "*** ERROR: GCTAResponseIrf::npsf:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (value=" << value;
                std::cout << ", roi_radius=" << roi_radius * gammalib::rad2deg;
                std::cout << ", roi_psf_distance=";
                std::cout << roi_psf_distance * gammalib::rad2deg;
                std::cout << ", r=[" << rmin * gammalib::rad2deg;
                std::cout << "," << rmax * gammalib::rad2deg << "])";
                std::cout << std::endl;
            }
            #endif

        } // endif: integration range was valid

    } // endelse: numerical integration required

    // Return integrated PSF
    return value;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAResponseIrf::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();
    m_aeff           = NULL;
    m_psf            = NULL;
    m_edisp          = NULL;
    m_background     = NULL;
    m_apply_edisp    = false;  //!< Switched off by default
    m_lo_safe_thres  = 0.0;
    m_hi_safe_thres  = 0.0;

    // XML response filenames
    m_xml_caldb.clear();
    m_xml_rspname.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GCTAResponseIrf::copy_members(const GCTAResponseIrf& rsp)
{
    // Copy members
    m_caldb          = rsp.m_caldb;
    m_rspname        = rsp.m_rspname;
    m_apply_edisp    = rsp.m_apply_edisp;
    m_lo_safe_thres  = rsp.m_lo_safe_thres;
    m_hi_safe_thres  = rsp.m_hi_safe_thres;

    // Copy response filenames
    m_xml_caldb      = rsp.m_xml_caldb;
    m_xml_rspname    = rsp.m_xml_rspname;

    // Copy nroi cache

    // Clone members
    m_aeff       = (rsp.m_aeff       != NULL) ? rsp.m_aeff->clone()  : NULL;
    m_psf        = (rsp.m_psf        != NULL) ? rsp.m_psf->clone()   : NULL;
    m_edisp      = (rsp.m_edisp      != NULL) ? rsp.m_edisp->clone() : NULL;
    m_background = (rsp.m_background != NULL) ? rsp.m_background->clone() : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponseIrf::free_members(void)
{
    // Free memory
    if (m_aeff       != NULL) delete m_aeff;
    if (m_psf        != NULL) delete m_psf;
    if (m_edisp      != NULL) delete m_edisp;
    if (m_background != NULL) delete m_background;

    // Initialise pointers
    m_aeff       = NULL;
    m_psf        = NULL;
    m_edisp      = NULL;
    m_background = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return filename with appropriate extension
 *
 * @param[in] filename File name.
 * @return File name.
 *
 * Checks if the specified @p filename exists, and if not, checks whether a
 * file with the added suffix .dat exists. Returns the file name with the
 * appropriate extension.
 ***************************************************************************/
GFilename GCTAResponseIrf::irf_filename(const std::string& filename) const
{
    // Set input filename as result filename
    GFilename result = filename;

    // If file does not exist then try a variant with extension .dat
    if (!result.exists()) {
        GFilename testname = filename + ".dat";
        if (testname.exists()) {
            result = testname;
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return value of point source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Value of instrument response function for a point source.
 *
 * This method returns the value of the instrument response function for a
 * point source. The method assumes that source.model() is of type
 * GModelSpatialPointSource.
 ***************************************************************************/
double GCTAResponseIrf::irf_ptsrc(const GEvent&       event,
                                  const GSource&      source,
                                  const GObservation& obs) const
{
    // Get point source spatial model
    const GModelSpatialPointSource* src =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Set Photon
    GPhoton photon(src->dir(), source.energy(), source.time());

    // Compute IRF
    double irf = this->irf(event, photon, obs);

    // Return IRF
    return irf;
}


/***********************************************************************//**
 * @brief Return IRF value for radial source model
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_model_type
 *            Model is not a radial model.
 *
 * Integrates the product of the spatial model and the instrument response
 * function over the true photon arrival direction using
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *    \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *    {\rm Irf}(\rho, \omega) d\omega d\rho
 * \f]
 *
 * where
 * \f$S_{\rm p}(\rho | E, t)\f$ is the radial spatial model,
 * \f${\rm Irf}(\rho, \omega)\f$ is the instrument response function,
 * \f$\rho\f$ is the radial distance from the model centre, and
 * \f$\omega\f$ is the position angle with respect to the connecting line
 * between the model centre and the observed photon arrival direction.
 *
 * The integration is performed in the coordinate system of the source
 * model spanned by \f$\rho\f$ and \f$\omega\f$ which allows to benefit
 * from the symmetry of the source model.
 *
 * The source centre is located at \f$\vec{m}\f$, and a spherical system
 * is defined around this location with \f$(\omega,\rho)\f$ being the
 * azimuth and zenith angles, respectively. \f$\omega=0\f$ is defined
 * by the direction that connects the source centre \f$\vec{m}\f$ to the
 * measured photon direction \f$\vec{p'}\f$, and \f$\omega\f$ increases
 * counterclockwise.
 *
 * Note that this method approximates the true theta angle (angle between
 * incident photon and pointing direction) by the measured theta angle
 * (angle between the measured photon arrival direction and the pointing
 * direction). Given the slow variation of the Psf shape over the field of
 * view, this approximation should be fine. It helps in fact a lot in
 * speeding up the computations.
 ***************************************************************************/
double GCTAResponseIrf::irf_radial(const GEvent&       event,
                                   const GSource&      source,
                                   const GObservation& obs) const
{
    // Retrieve CTA pointing
    const GCTAPointing& pnt = gammalib::cta_pnt(G_IRF_RADIAL, obs);
    const GCTAInstDir&  dir = gammalib::cta_dir(G_IRF_RADIAL, event);

    // Get pointer on radial model
    const GModelSpatialRadial* model =
          dynamic_cast<const GModelSpatialRadial*>(source.model());
    if (model == NULL) {
        throw GCTAException::bad_model_type(G_IRF_RADIAL);
    }

    // Get pointer on shell model (will be NULL for other models)
    const GModelSpatialRadialShell* shell =
          dynamic_cast<const GModelSpatialRadialShell*>(model);

    // Get pointer on ring model (will be NULL for other models)
    const GModelSpatialRadialRing* ring =
          dynamic_cast<const GModelSpatialRadialRing*>(model);

    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1299
    // https://cta-redmine.irap.omp.eu/issues/1521
    // https://cta-redmine.irap.omp.eu/issues/2860
    static const int iter_rho = 6;
    static const int iter_phi = 6;

    // Get event attributes
    const GSkyDir& obsDir = dir.dir();
    const GEnergy& obsEng = event.energy();

    // Get source attributes
    const GSkyDir& centre  = model->dir();
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Determine angular distance between measured photon direction and model
    // centre [radians]
    double zeta = centre.dist(obsDir);

    // Determine angular distance between measured photon direction and
    // pointing direction [radians]
    double eta = pnt.dir().dist(obsDir);

    // Determine angular distance between model centre and pointing direction
    // [radians]
    double lambda = centre.dist(pnt.dir());

    // Compute azimuth angle of pointing in model system [radians]
    // Will be comprised in interval [0,pi]
    double omega0 = 0.0;
    double denom  = std::sin(lambda) * std::sin(zeta);
    if (denom != 0.0) {
        double arg = (std::cos(eta) - std::cos(lambda) * std::cos(zeta))/denom;
        omega0     = gammalib::acos(arg);
    }

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Assign the observed theta angle (eta) as the true theta angle
    // between the source and the pointing directions. This is a (not
    // too bad) approximation which helps to speed up computations.
    // If we want to do this correctly, however, we would need to move
    // the psf_dummy_sigma down to the integration kernel, and we would
    // need to make sure that psf_delta_max really gives the absolute
    // maximum (this is certainly less critical)
    double theta = eta;
    double phi   = 0.0; //TODO: Implement IRF Phi dependence

    // Get maximum PSF and source radius in radians.
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);
    double src_max   = model->theta_max();

    // Set radial model zenith angle range
    double rho_min = (zeta > delta_max) ? zeta - delta_max : 0.0;
    double rho_max = zeta + delta_max;
    if (rho_max > src_max) {
        rho_max = src_max;
    }

    // Initialise IRF value
    double irf = 0.0;

    // Perform zenith angle integration if interval is valid
    if (rho_max > rho_min) {

        // Setup integration kernel
        cta_irf_radial_kern_rho integrand(this,
                                          model,
                                          zenith,
                                          azimuth,
                                          srcEng,
                                          srcTime,
                                          obsEng,
                                          zeta,
                                          lambda,
                                          omega0,
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
        double transition_point = delta_max - zeta;
        if (transition_point > rho_min && transition_point < rho_max) {
            bounds.push_back(transition_point);
        }

        // If we have a shell model then add an integration boundary for the
        // shell radius as a function discontinuity will occur at this
        // location
        if (shell != NULL) {
            double shell_radius = shell->radius() * gammalib::deg2rad;
            if (shell_radius > rho_min && shell_radius < rho_max) {
                bounds.push_back(shell_radius);
            }
        }

        // If we have a ring model then add an integration boundary for the
        // ring radius as a function discontinuity will occur at this
        // location
        if (ring != NULL) {
            double ring_radius = ring->radius() * gammalib::deg2rad;
            if (ring_radius > rho_min && ring_radius < rho_max) {
                bounds.push_back(ring_radius);
            }
        }

        // Integrate kernel
        irf = integral.romberg(bounds, iter_rho);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GCTAResponseIrf::irf_radial:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", rho_min=" << rho_min;
            std::cout << ", rho_max=" << rho_max;
            std::cout << ", omega0=" << omega0 << ")";
            std::cout << std::endl;
        }
        #endif

        // Apply deadtime correction
        irf *= obs.deadc(srcTime);

    } // endif: integration interval is valid

    // Compile option: Show integration results
    #if defined(G_DEBUG_IRF_RADIAL)
    std::cout << "GCTAResponseIrf::irf_radial:";
    std::cout << " rho=[" << rho_min*gammalib::rad2deg;
    std::cout << ", " << rho_max*gammalib::rad2deg << " deg;";
    std::cout << " zeta=" << zeta*gammalib::rad2deg << " deg;";
    std::cout << " eta=" << eta*gammalib::rad2deg << " deg;";
    std::cout << " lambda=" << lambda*gammalib::rad2deg << " deg;";
    std::cout << " omega0=" << omega0*gammalib::rad2deg << " deg;";
    std::cout << " theta_max=" << src_max*gammalib::rad2deg << " deg;";
    std::cout << " delta_max=" << delta_max*gammalib::rad2deg << " deg;";
    std::cout << " obsDir=" << obsDir << ";";
    std::cout << " modelDir=" << model->dir() << ";";
    std::cout << " irf=" << irf << std::endl;
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return Irf value for elliptical source model
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_instdir_type
 *            Instrument direction is not a valid CTA instrument direction.
 * @exception GCTAException::bad_model_type
 *            Model is not an elliptical model.
 *
 * Integrates the product of the model and the IRF over the true photon
 * arrival direction using
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times
 *    \int_{\omega} 
 *    S_{\rm p}(\rho, \omega | E, t) \, IRF(\rho, \omega) d\omega d\rho
 * \f]
 *
 * where
 * \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical model,
 * \f$IRF(\rho, \omega)\f$ is the instrument response function,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the position angle with respect to the connecting line
 * between the model centre and the observed photon arrival direction.
 *
 * The source model centre is located at \f$\vec{m}\f$, and a spherical
 * coordinate system is defined around this location with \f$(\rho,\omega)\f$
 * being the zenith and azimuth angles, respectively. The azimuth angle
 * \f$(\omega)\f$ is counted counterclockwise from the vector that runs from
 * the model centre \f$\vec{m}\f$ to the measured photon direction
 * \f$\vec{p'}\f$.
 ***************************************************************************/
double GCTAResponseIrf::irf_elliptical(const GEvent&       event,
                                       const GSource&      source,
                                       const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1299
    // https://cta-redmine.irap.omp.eu/issues/1521
    // https://cta-redmine.irap.omp.eu/issues/2860
    static const int iter_rho = 6;
    static const int iter_phi = 6;

    // Retrieve CTA pointing
    const GCTAPointing& pnt = gammalib::cta_pnt(G_IRF_ELLIPTICAL, obs);
    const GCTAInstDir&  dir = gammalib::cta_dir(G_IRF_ELLIPTICAL, event);

    // Get pointer on elliptical model
    const GModelSpatialElliptical* model =
          dynamic_cast<const GModelSpatialElliptical*>(source.model());
    if (model == NULL) {
        throw GCTAException::bad_model_type(G_IRF_ELLIPTICAL);
    }

    // Get event attributes (measured photon)
    const GSkyDir& obsDir = dir.dir();
    const GEnergy& obsEng = event.energy();

    // Get source attributes
    const GSkyDir& centre  = model->dir();
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Determine angular distance between observed photon direction and model
    // centre and position angle of observed photon direction seen from the
    // model centre [radians]
    double rho_obs      = centre.dist(obsDir);
    double posangle_obs = centre.posang(obsDir); // Celestial system

    // Determine angular distance between model centre and pointing direction
    // [radians]
    double rho_pnt      = centre.dist(pnt.dir());
    double posangle_pnt = centre.posang(pnt.dir()); // Celestial system

    // Compute azimuth angle of pointing in model coordinate system [radians]
    double omega_pnt = posangle_pnt - posangle_obs;

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius [radians]. We assign here the measured theta
    // angle (eta) as the true theta angle between the source and the pointing
    // directions. As we only use the angle to determine the maximum PSF size,
    // this should be sufficient.
    double theta     = pnt.dir().dist(obsDir);
    double phi       = 0.0; //TODO: Implement IRF Phi dependence
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

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

    // Initialise IRF value
    double irf = 0.0;

    // Perform zenith angle integration if interval is valid
    if (rho_max > rho_min) {

        // Setup integration kernel
        cta_irf_elliptical_kern_rho integrand(this,
                                              model,
                                              semimajor,
                                              semiminor,
                                              posangle,
                                              zenith,
                                              azimuth,
                                              srcEng,
                                              srcTime,
                                              obsEng,
                                              rho_obs,
                                              posangle_obs,
                                              rho_pnt,
                                              omega_pnt,
                                              delta_max,
                                              iter_phi);

        // Integrate over model's zenith angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter_rho);

        // Setup integration boundaries
        std::vector<double> bounds;
        bounds.push_back(rho_min);
        bounds.push_back(rho_max);

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
            std::cout << "*** ERROR: GCTAResponseIrf::irf_elliptical:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", rho_min=" << rho_min;
            std::cout << ", rho_max=" << rho_max << ")";
            std::cout << std::endl;
        }
        #endif

        // Apply deadtime correction
        irf *= obs.deadc(srcTime);

    } // endif: integration interval is valid

    // Compile option: Show integration results
    #if defined(G_DEBUG_IRF_ELLIPTICAL)
    std::cout << "GCTAResponseIrf::irf_elliptical:";
    std::cout << " rho_min=" << rho_min;
    std::cout << " rho_max=" << rho_max;
    std::cout << " irf=" << irf << std::endl;
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of diffuse source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * Integrates the product of the model and the IRF over the true photon
 * arrival direction using
 *
 * \f[
 *    \int_{0}^{\theta_{\rm max}}
 *    \sin \theta \times PSF(\theta)
 *    \int_{0}^{2\pi}
 *    S_{\rm p}(\theta, \phi | E, t) \,
 *    Aeff(\theta, \phi) \,
 *    Edisp(\theta, \phi) d\phi
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\theta, \phi | E, t)\f$ is the diffuse model,
 * - \f$PSF(\theta)\f$ is the azimuthally symmetric Point Spread Function,
 * - \f$Aeff(\theta, \phi)\f$ is the effective area,
 * - \f$Edisp(\theta, \phi)\f$ is the energy dispersion,
 * - \f$\theta\f$ is the distance from the PSF centre, and
 * - \f$\phi\f$ is the azimuth angle.
 *
 * The integration is performed in the reference of the observed arrival
 * direction. Integration is done first over the azimuth angle \f$\phi\f$ and
 * then over the offset angle \f$\theta\f$.
 *
 * The integration kernels for this method are implemented by the response
 * helper classes cta_irf_diffuse_kern_theta and cta_irf_diffuse_kern_phi.
 ***************************************************************************/
double GCTAResponseIrf::irf_diffuse(const GEvent&       event,
                                    const GSource&      source,
                                    const GObservation& obs) const
{
    // Set minimum and maximum number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1248
    // https://cta-redmine.irap.omp.eu/issues/2458
    // https://cta-redmine.irap.omp.eu/issues/2860
    static const int min_iter_rho = 6;
    static const int min_iter_phi = 6;
    static const int max_iter_rho = 8;
    static const int max_iter_phi = 8;

    // Initialise IRF value
    double irf = 0.0;

    // Retrieve CTA pointing
    const GCTAPointing& pnt = gammalib::cta_pnt(G_IRF_DIFFUSE, obs);

    // Get CTA instrument direction
    const GCTAInstDir& dir = gammalib::cta_dir(G_IRF_DIFFUSE, event);

    // Get pointer on spatial model
    const GModelSpatial* model =
        dynamic_cast<const GModelSpatial*>(source.model());
    if (model == NULL) {
        throw GCTAException::bad_model_type(G_IRF_DIFFUSE);
    }

    // Get resolution of spatial model
    double resolution = gammalib::resolution(model);

    // Get event attributes
    const GEnergy& obsEng = event.energy();

    // Get source attributes
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Determine angular distance between measured photon direction and
    // pointing direction [radians]
    double eta = pnt.dir().dist(dir.dir());

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Assign the observed theta angle (eta) as the true theta angle
    // between the source and the pointing directions. This is a (not
    // too bad) approximation which helps to speed up computations.
    // If we want to do this correctly, however, we would need to move
    // the psf_dummy_sigma down to the integration kernel, and we would
    // need to make sure that psf_delta_max really gives the absolute
    // maximum (this is certainly less critical)
    double theta = eta;
    double phi   = 0.0; //TODO: Implement Phi dependence

    // Get maximum PSF radius in radians
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

    // Perform zenith angle integration if interval is valid
    if (delta_max > 0.0) {

        // Compute rotation matrix to convert from coordinates (theta,phi)
        // in the reference frame of the observed arrival direction into
        // celestial coordinates
        GMatrix ry;
        GMatrix rz;
        ry.eulery(dir.dir().dec_deg() - 90.0);
        rz.eulerz(-dir.dir().ra_deg());
        GMatrix rot = (ry * rz).transpose();

        // Setup integration kernel
        cta_irf_diffuse_kern_theta integrand(this,
                                             model,
                                             &rot,
                                             theta,
                                             phi,
                                             zenith,
                                             azimuth,
                                             srcEng,
                                             srcTime,
                                             srcLogEng,
                                             obsEng,
                                             eta,
                                             min_iter_phi,
                                             max_iter_phi,
                                             resolution);

        // Set number of radial integration iterations
        int iter  = gammalib::iter_rho(delta_max, resolution,
                                       min_iter_rho, max_iter_rho);

        // Integrate over Psf delta angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter);
        irf = integral.romberg(0.0, delta_max);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GCTAResponseIrf::irf_diffuse:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", delta_max=" << delta_max << ")";
            std::cout << std::endl;
        }
        #endif
    }

    // Apply deadtime correction
    irf *= obs.deadc(srcTime);

    // Compile option: Show integration results
    #if defined(G_DEBUG_IRF_DIFFUSE)
    std::cout << "GCTAResponseIrf::irf_diffuse:";
    std::cout << " srcLogEng=" << srcLogEng;
    std::cout << " obsLogEng=" << obsLogEng;
    std::cout << " eta=" << eta;
    std::cout << " delta_max=" << delta_max;
    std::cout << " irf=" << irf << std::endl;
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return spatial integral of point source model
 *
 * @param[in] model Sky Model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event arrival time.
 * @param[in] obs Observation.
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t'|E,t) = \int_{\rm ROI} P(p',E',t'|E,t) dp'
 * \f]
 *
 * of
 *
 * \f[
 *    P(p',E',t'|E,t) = \int
 *                      S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * over the Region of Interest (ROI) for a point source model \f$S(p,E,t)\f$
 * and the response function \f$R(p',E',t'|p,E,t)\f$.
 ***************************************************************************/
double GCTAResponseIrf::nroi_ptsrc(const GModelSky&    model,
                                   const GEnergy&      srcEng,
                                   const GTime&        srcTime,
                                   const GEnergy&      obsEng,
                                   const GTime&        obsTime,
                                   const GObservation& obs) const
{
    // Get point source spatial model
    const GModelSpatialPointSource* src =
          static_cast<const GModelSpatialPointSource*>(model.spatial());

    // Set Photon
    GPhoton photon(src->dir(), srcEng, srcTime);

    // Compute Nroi
    double nroi = nirf(photon, obsEng, obsTime, obs);

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Return spatial integral of radial source model
 *
 * @param[in] model Sky Model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event arrival time.
 * @param[in] obs Observation.
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t'|E,t) = \int_{\rm ROI} P(p',E',t'|E,t) dp'
 * \f]
 *
 * of
 *
 * \f[
 *    P(p',E',t'|E,t) = \int
 *                      S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * over the Region of Interest (ROI) for a radial source model \f$S(p,E,t)\f$
 * and the response function \f$R(p',E',t'|p,E,t)\f$.
 ***************************************************************************/
double GCTAResponseIrf::nroi_radial(const GModelSky&    model,
                                    const GEnergy&      srcEng,
                                    const GTime&        srcTime,
                                    const GEnergy&      obsEng,
                                    const GTime&        obsTime,
                                    const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1299
    static const int iter_rho = 6;
    static const int iter_phi = 6;

    // Initialise Nroi value
    double nroi = 0.0;

    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = gammalib::cta_obs(G_NROI_RADIAL, obs);
    const GCTARoi&         roi = gammalib::cta_event_list(G_NROI_RADIAL, obs).roi();
    const GCTAPointing&    pnt = cta.pointing();

    // Get pointer on radial model
    const GModelSpatialRadial* spatial =
          dynamic_cast<const GModelSpatialRadial*>(model.spatial());
    if (spatial == NULL) {
        throw GCTAException::bad_model_type(G_NROI_RADIAL);
    }

    // Get source attributes
    const GSkyDir& centre  = spatial->dir();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius (radians). We do this for the onaxis PSF only,
    // as this allows us doing this computation in the outer loop. This
    // should be sufficient here, unless the offaxis PSF becomes much worse
    // than the onaxis PSF. In this case, we may add a safety factor here
    // to make sure we encompass the entire PSF.
    double psf_max_radius = psf_delta_max(0.0, 0.0, zenith, azimuth, srcLogEng);

    // Extract ROI radius (radians)
    double roi_radius = roi.radius() * gammalib::deg2rad;

    // Compute distance between ROI and model centre (radians)
    double roi_model_distance = roi.centre().dir().dist(centre);

    // Compute the ROI radius plus maximum PSF radius (radians). Any photon
    // coming from beyond this radius will not make it in the dataspace and
    // thus can be neglected.
    double roi_psf_radius = roi_radius + psf_max_radius;

    // Set offset angle integration range. We take here the ROI+PSF into
    // account to make no integrations beyond the point where the
    // contribution drops to zero.
    double rho_min = (roi_model_distance > roi_psf_radius)
                     ? roi_model_distance - roi_psf_radius: 0.0;
    double rho_max = spatial->theta_max();

    // Perform offset angle integration only if interval is valid
    if (rho_max > rho_min) {

        // Compute rotation matrix to convert from native model coordinates,
        // given by (rho,omega), into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        ry.eulery(spatial->dec() - 90.0);
        rz.eulerz(-spatial->ra());
        GMatrix rot = (ry * rz).transpose();

        // Compute position angle of ROI centre with respect to model
        // centre (radians)
        double omega0 = centre.posang(roi.centre().dir()); // Celestial system

        // Setup integration kernel
        cta_nroi_radial_kern_rho integrand(this,
                                           &cta,
                                           spatial,
                                           &rot,
                                           srcEng,
                                           srcTime,
                                           obsEng,
                                           obsTime,
                                           roi_model_distance,
                                           roi_psf_radius,
                                           omega0,
                                           iter_phi);

        // Integrate over model's zenith angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter_rho);

        // Setup integration boundaries
        std::vector<double> bounds;
        bounds.push_back(rho_min);
        bounds.push_back(rho_max);

        // If the integration range includes a transition between full
        // containment of model within ROI and partial containment, then
        // add a boundary at this location
        double transition_point = roi_psf_radius - roi_model_distance;
        if (transition_point > rho_min && transition_point < rho_max) {
            bounds.push_back(transition_point);
        }

        // If the integration range includes a transition between full
        // containment of ROI within model and partial containment, then
        // add a boundary at this location
        transition_point = roi_psf_radius + roi_model_distance;
        if (transition_point > rho_min && transition_point < rho_max) {
            bounds.push_back(transition_point);
        }

        // If we have a shell model then add an integration boundary for the
        // shell radius as a function discontinuity will occur at this
        // location
        const GModelSpatialRadialShell* shell = dynamic_cast<const GModelSpatialRadialShell*>(spatial);
        if (shell != NULL) {
            double shell_radius = shell->radius() * gammalib::deg2rad;
            if (shell_radius > rho_min && shell_radius < rho_max) {
                bounds.push_back(shell_radius);
            }
        }

        // If we have a ring model then add an integration boundary for the
        // ring radius as a function discontinuity will occur at this
        // location
        const GModelSpatialRadialRing* ring = dynamic_cast<const GModelSpatialRadialRing*>(spatial);
        if (ring != NULL) {
            double ring_radius = ring->radius() * gammalib::deg2rad;
            if (ring_radius > rho_min && ring_radius < rho_max) {
                bounds.push_back(ring_radius);
            }
        }

        // Integrate kernel
        nroi = integral.romberg(bounds, iter_rho);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NROI_RADIAL)
        std::cout << "GCTAResponseIrf::nroi_radial:";
        std::cout << " rho_min=" << rho_min;
        std::cout << " rho_max=" << rho_max;
        std::cout << " nroi=" << nroi << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
        std::cout << "*** ERROR: GCTAResponseIrf::nroi_radial:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (nroi=" << nroi;
        std::cout << ", rho_min=" << rho_min;
        std::cout << ", rho_max=" << rho_max;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Return spatial integral of elliptical source model
 *
 * @param[in] model Sky Model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event arrival time.
 * @param[in] obs Observation.
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t'|E,t) = \int_{\rm ROI} P(p',E',t'|E,t) dp'
 * \f]
 *
 * of
 *
 * \f[
 *    P(p',E',t'|E,t) = \int
 *                      S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * over the Region of Interest (ROI) for an elliptical source model
 * \f$S(p,E,t)\f$ and the response function \f$R(p',E',t'|p,E,t)\f$.
 ***************************************************************************/
double GCTAResponseIrf::nroi_elliptical(const GModelSky&    model,
                                        const GEnergy&      srcEng,
                                        const GTime&        srcTime,
                                        const GEnergy&      obsEng,
                                        const GTime&        obsTime,
                                        const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1299
    static const int iter_rho = 6;
    static const int iter_phi = 6;

    // Initialise Nroi value
    double nroi = 0.0;

    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = gammalib::cta_obs(G_NROI_ELLIPTICAL, obs);
    const GCTARoi&         roi = gammalib::cta_event_list(G_NROI_ELLIPTICAL, obs).roi();
    const GCTAPointing&    pnt = cta.pointing();

    // Get pointer on elliptical model
    const GModelSpatialElliptical* spatial =
          dynamic_cast<const GModelSpatialElliptical*>(model.spatial());
    if (spatial == NULL) {
        throw GCTAException::bad_model_type(G_NROI_ELLIPTICAL);
    }

    // Get source attributes
    const GSkyDir& centre  = spatial->dir();

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius (radians). We do this for the onaxis PSF only,
    // as this allows us doing this computation in the outer loop. This
    // should be sufficient here, unless the offaxis PSF becomes much worse
    // than the onaxis PSF. In this case, we may add a safety factor here
    // to make sure we encompass the entire PSF.
    double psf_max_radius = psf_delta_max(0.0, 0.0, zenith, azimuth, srcLogEng);

    // Extract ROI radius plus maximum PSF radius (radians). Any photon
    // coming from beyond this radius will not make it in the dataspace and
    // thus can be neglected.
    double radius_roi = roi.radius() * gammalib::deg2rad + psf_max_radius;

    // Compute distance between ROI and model centre (radians)
    double rho_roi = roi.centre().dir().dist(centre);

    // Get the ellipse boundary (radians). Note that these are NOT the
    // parameters of the ellipse but of a boundary ellipse that is used
    // for computing the relevant omega angle intervals for a given angle
    // rho. The boundary ellipse takes care of the possibility that the
    // semiminor axis is larger than the semimajor axis
    double semimajor;    // Will be the larger axis
    double semiminor;    // Will be the smaller axis
    double posangle;     // Will be the corrected position angle
    double aspect_ratio; // Ratio between smaller/larger axis of model
    if (spatial->semimajor() >= spatial->semiminor()) {
        aspect_ratio = (spatial->semimajor() > 0.0) ?
                        spatial->semiminor() / spatial->semimajor() : 0.0;
        posangle     = spatial->posangle() * gammalib::deg2rad;
    }
    else {
        aspect_ratio = (spatial->semiminor() > 0.0) ?
                        spatial->semimajor() / spatial->semiminor() : 0.0;
        posangle     = spatial->posangle() * gammalib::deg2rad + gammalib::pihalf;
    }
    semimajor = spatial->theta_max();
    semiminor = semimajor * aspect_ratio;

    // Set offset angle integration range. We take here the ROI+PSF into
    // account to make no integrations beyond the point where the
    // contribution drops to zero.
    double rho_min = (rho_roi > radius_roi) ? rho_roi - radius_roi: 0.0;
    double rho_max = rho_roi + radius_roi;
    if (rho_max > semimajor) {
        rho_max = semimajor;
    }

    // Perform offset angle integration only if interval is valid
    if (rho_max > rho_min) {

        // Compute rotation matrix to convert from native model coordinates,
        // given by (rho,omega), into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        ry.eulery(spatial->dec() - 90.0);
        rz.eulerz(-spatial->ra());
        GMatrix rot = (ry * rz).transpose();

        // Compute position angle of ROI centre with respect to model
        // centre (radians)
        double posangle_roi = centre.posang(roi.centre().dir()); // Celestial system

        // Setup integration kernel
        cta_nroi_elliptical_kern_rho integrand(this,
                                               &cta,
                                               spatial,
                                               &rot,
                                               semimajor,
                                               semiminor,
                                               posangle,
                                               srcEng,
                                               srcTime,
                                               obsEng,
                                               obsTime,
                                               rho_roi,
                                               posangle_roi,
                                               radius_roi,
                                               iter_phi);

        // Integrate over model's zenith angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter_rho);

        // Setup integration boundaries
        std::vector<double> bounds;
        bounds.push_back(rho_min);
        bounds.push_back(rho_max);

        // If the integration range includes the semiminor boundary, then
        // add an integration boundary at that location
        if (semiminor > rho_min && semiminor < rho_max) {
            bounds.push_back(semiminor);
        }

        // Integrate kernel
        nroi = integral.romberg(bounds, iter_rho);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NROI_ELLIPTICAL)
        std::cout << "GCTAResponseIrf::nroi_elliptical:";
        std::cout << " rho_min=" << rho_min;
        std::cout << " rho_max=" << rho_max;
        std::cout << " nroi=" << nroi << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
        std::cout << "*** ERROR: GCTAResponseIrf::nroi_elliptical:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (nroi=" << nroi;
        std::cout << ", rho_min=" << rho_min;
        std::cout << ", rho_max=" << rho_max;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Return spatial integral of diffuse source model
 *
 * @param[in] model Sky Model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event arrival time.
 * @param[in] obs Observation.
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t'|E,t) = \int_{\rm ROI} P(p',E',t'|E,t) dp'
 * \f]
 *
 * of
 *
 * \f[
 *    P(p',E',t'|E,t) = \int
 *                      S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * over the Region of Interest (ROI) for a diffuse source model
 * \f$S(p,E,t)\f$ and the response function \f$R(p',E',t'|p,E,t)\f$.
 ***************************************************************************/
double GCTAResponseIrf::nroi_diffuse(const GModelSky&    model,
                                     const GEnergy&      srcEng,
                                     const GTime&        srcTime,
                                     const GEnergy&      obsEng,
                                     const GTime&        obsTime,
                                     const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    // These values have been determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1248
    static const int iter_rho = 9;
    static const int iter_phi = 9;

    // Initialise Nroi value
    double nroi = 0.0;

    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = gammalib::cta_obs(G_NROI_DIFFUSE, obs);
    const GCTARoi&         roi = gammalib::cta_event_list(G_NROI_DIFFUSE, obs).roi();
    const GCTAPointing&    pnt = cta.pointing();

    // Get pointer on spatial model
    const GModelSpatial* spatial =
        dynamic_cast<const GModelSpatial*>(model.spatial());
    if (spatial == NULL) {
        throw GCTAException::bad_model_type(G_NROI_DIFFUSE);
    }

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt.zenith();
    double azimuth = pnt.azimuth();

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius (radians). We do this for the onaxis PSF only,
    // as this allows us doing this computation in the outer loop. This
    // should be sufficient here, unless the offaxis PSF becomes much worse
    // than the onaxis PSF. In this case, we may add a safety factor here
    // to make sure we encompass the entire PSF.
    double psf_max_radius = psf_delta_max(0.0, 0.0, zenith, azimuth, srcLogEng);

    // Extract ROI radius (radians)
    double roi_radius = roi.radius() * gammalib::deg2rad;

    // Compute the ROI radius plus maximum PSF radius (radians). Any photon
    // coming from beyond this radius will not make it in the dataspace and
    // thus can be neglected.
    double roi_psf_radius = roi_radius + psf_max_radius;

    // Perform offset angle integration only if interval is valid
    if (roi_psf_radius > 0.0) {

        // Compute rotation matrix to convert from native ROI coordinates,
        // given by (theta,phi), into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        ry.eulery(roi.centre().dir().dec_deg() - 90.0);
        rz.eulerz(-roi.centre().dir().ra_deg());
        GMatrix rot = (ry * rz).transpose();

        // Setup integration kernel
        cta_nroi_diffuse_kern_theta integrand(this,
                                              &cta,
                                              spatial,
                                              &rot,
                                              srcEng,
                                              srcTime,
                                              obsEng,
                                              obsTime,
                                              iter_phi);

        // Integrate over model's zenith angle
        GIntegral integral(&integrand);
        integral.fixed_iter(iter_rho);
        nroi = integral.romberg(0.0, roi_psf_radius);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NROI_DIFFUSE)
        std::cout << "GCTAResponseIrf::nroi_diffuse:";
        std::cout << " roi_psf_radius=" << roi_psf_radius;
        std::cout << " nroi=" << nroi;
        std::cout << " Etrue=" << srcEng;
        std::cout << " Ereco=" << obsEng;
        std::cout << " id=" << id << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
        std::cout << "*** ERROR: GCTAResponseIrf::nroi_diffuse:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (nroi=" << nroi;
        std::cout << ", roi_psf_radius=" << roi_psf_radius;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Return spatial integral of composite source model
 *
 * @param[in] model Sky Model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obsEng Observed event energy.
 * @param[in] obsTime Observed event arrival time.
 * @param[in] obs Observation.
 *
 * Computes the integral
 *
 * \f[
 *    N_{\rm ROI}(E',t'|E,t) = \int_{\rm ROI} P(p',E',t'|E,t) dp'
 * \f]
 *
 * of
 *
 * \f[
 *    P(p',E',t'|E,t) = \int
 *                      S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * over the Region of Interest (ROI) for a composite source model
 * \f$S(p,E,t)\f$ and the response function \f$R(p',E',t'|p,E,t)\f$.
 ***************************************************************************/
double GCTAResponseIrf::nroi_composite(const GModelSky&    model,
                                       const GEnergy&      srcEng,
                                       const GTime&        srcTime,
                                       const GEnergy&      obsEng,
                                       const GTime&        obsTime,
                                       const GObservation& obs) const
{
    // Initialise nroi
    double nroi = 0.0;

    // Get composite model
    GModelSpatialComposite* comp =
        dynamic_cast<GModelSpatialComposite*>(model.spatial());

    // Loop over model components
    for (int i = 0; i < comp->components(); ++i) {

        // Create new sky model
        GModelSky sky = GModelSky(model);

        // Assign spatial part
        sky.spatial(comp->component(i));

        // Compute nroi
        nroi += this->nroi(sky, srcEng, srcTime, obsEng, obsTime, obs) *
                comp->scale(i);

    }

    // Divide by number of model components
    double sum = comp->sum_of_scales();
    if (sum > 0.0) {
        nroi /= sum;
    }

    // Return nroi
    return nroi;

}
