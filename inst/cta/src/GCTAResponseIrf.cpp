/***************************************************************************
 *       GCTAResponseIrf.cpp - CTA instrument response function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#include "GIntegral.hpp"
#include "GCaldb.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAPointing.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAEventList.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"
#include "GCTAAeff2D.hpp"
#include "GCTAAeffArf.hpp"
#include "GCTAAeffPerfTable.hpp"
#include "GCTAPsf2D.hpp"
#include "GCTAPsfVector.hpp"
#include "GCTAPsfPerfTable.hpp"
#include "GCTAPsfKing.hpp"
#include "GCTAAeff.hpp"
#include "GCTAPsf.hpp"
#include "GCTAEdisp.hpp"
#include "GCTAEdispRmf.hpp"
#include "GCTAEdispPerfTable.hpp"
#include "GCTABackground.hpp"
#include "GCTABackgroundPerfTable.hpp"
#include "GCTABackground3D.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                             "GCTAResponseIrf::caldb(GCaldb&)"
#define G_IRF   "GCTAResponseIrf::irf(GInstDir&, GEnergy&, GTime&, GSkyDir&,"\
                                          " GEnergy&, GTime&, GObservation&)"
#define G_NPRED          "GCTAResponseIrf::npred(GSkyDir&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC   "GCTAResponseIrf::mc(double&, GPhoton&, GObservation&, GRan&)"
#define G_IRF_RADIAL         "GCTAResponseIrf::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_ELLIPTICAL "GCTAResponseIrf::irf_elliptical(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE       "GCTAResponseIrf::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_NPRED_RADIAL              "GCTAResponseIrf::npred_radial(GSource&,"\
                                                            " GObservation&)"
#define G_NPRED_ELLIPTICAL      "GCTAResponseIrf::npred_elliptical(GSource&,"\
                                                            " GObservation&)"
#define G_NPRED_DIFFUSE            "GCTAResponseIrf::npred_diffuse(GSource&,"\
                                                            " GObservation&)"
#define G_AEFF    "GCTAResponseIrf::aeff(double&, double&, double&, double&,"\
                                                                  " double&)"
#define G_PSF      "GCTAResponseIrf::psf(double&, double&, double&, double&,"\
                                                                  " double&)"
#define G_PSF_DELTA_MAX    "GCTAResponseIrf::psf_delta_max(double&, double&,"\
                                                " double&, double&, double&)"
#define G_EDISP  "GCTAResponseIrf::edisp(double&, double&, double&, double&,"\
                                                                  " double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_IRF_CACHE            //!< Use IRF cache in irf_diffuse method
#define G_USE_NPRED_CACHE      //!< Use Npred cache in npred_diffuse method

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_IRF_RADIAL                     //!< Debug irf_radial method
//#define G_DEBUG_IRF_DIFFUSE                   //!< Debug irf_diffuse method
//#define G_DEBUG_IRF_ELLIPTICAL             //!< Debug irf_elliptical method
//#define G_DEBUG_NPRED_RADIAL                 //!< Debug npred_radial method
//#define G_DEBUG_NPRED_DIFFUSE               //!< Debug npred_diffuse method
//#define G_DEBUG_NPRED_ELLIPTICAL         //!< Debug npred_elliptical method
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
GCTAResponseIrf::GCTAResponseIrf(void) : GResponse()
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
GCTAResponseIrf::GCTAResponseIrf(const GCTAResponseIrf& rsp) : GResponse(rsp)
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
 * @param[in] rspname Response file name.
 * @param[in] caldb Calibration database.
 *
 * Create instance of CTA response by specifying the response name and the
 * calibration database. The response name can be either a response identifier
 * or a filename (see GCTAResponseIrf::load for more information).
 ***************************************************************************/
GCTAResponseIrf::GCTAResponseIrf(const std::string& rspname,
                                 const GCaldb& caldb) : GResponse()
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
void GCTAResponseIrf::clear(void)
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
    const GCTAPointing& pnt = retrieve_pnt(G_IRF, obs);
    const GCTAInstDir&  dir = retrieve_dir(G_IRF, event);

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

                // Get log10(E/TeV) of measured photon energy.
                //double obsLogEng = obsEng.log10TeV();

                // Multiply-in energy dispersion
                irf *= edisp(obsEng, theta, phi, zenith, azimuth, srcLogEng);

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
 * @brief Return spatial integral of point spread function
 *
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 *
 * @todo Set polar angle phi of photon in camera system
 * @todo Write method documentation
 ***************************************************************************/
double GCTAResponseIrf::npred(const GPhoton&      photon,
                              const GObservation& obs) const
{
    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = retrieve_obs(G_NPRED, obs);
    const GCTARoi&         roi = retrieve_roi(G_NPRED, obs);
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
    double npred = aeff(theta, phi, zenith, azimuth, srcLogEng);

    // Multiply-in PSF
    if (npred > 0.0) {

        // Get PSF
        npred *= npsf(srcDir, srcLogEng, srcTime, pnt, roi);

        // Multiply-in energy dispersion
        if (use_edisp() && npred > 0.0) {

            // Get energy dispersion
            npred *= nedisp(srcDir, srcEng, srcTime, pnt, cta.events()->ebounds());

        } // endif: had energy dispersion

        // Apply deadtime correction
        npred *= obs.deadc(srcTime);

    } // endif: had non-zero effective area

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GCTAResponseIrf::npred:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", theta=" << theta;
        std::cout << ", phi=" << phi << ")";
        std::cout << std::endl;
    }
    #endif

    // Return Npred
    return npred;
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
GCTAEventAtom* GCTAResponseIrf::mc(const double& area, const GPhoton& photon,
                                   const GObservation& obs, GRan& ran) const
{
    // Initialise event
    GCTAEventAtom* event = NULL;

    // Retrieve CTA pointing
    const GCTAPointing& pnt = retrieve_pnt(G_MC, obs);

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
    if (ran.uniform() <= ulimite) {

        // Apply deadtime correction
        double deadc = obs.deadc(photon.time());
        if (deadc >= 1.0 || ran.uniform() <= deadc) {

            // Simulate offset from photon arrival direction
            double delta = psf()->mc(ran, srcLogEng, theta, phi, zenith, azimuth) *
                           gammalib::rad2deg;
            double alpha = 360.0 * ran.uniform();

            // Rotate sky direction by offset
            GSkyDir sky_dir = photon.dir();
            sky_dir.rotate_deg(alpha, delta);

            // Set measured photon arrival direction
            GCTAInstDir inst_dir;
            inst_dir.dir(sky_dir);

            // Set measured photon energy
            GEnergy energy = photon.energy();
            if (use_edisp()) {
                energy = edisp()->mc(ran, srcLogEng, theta, phi, zenith, azimuth);
            }

            // Allocate event
            event = new GCTAEventAtom;

            // Set event attributes
            event->dir(inst_dir);
            event->energy(energy);
            event->time(photon.time());

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
 * @todo Still supports old ARF, PSF and RMF parameter names.
 ***************************************************************************/
void GCTAResponseIrf::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Extract parameters
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle Calibration
        if (par->attribute("name") == "Calibration") {

            // Read database and response
            m_xml_caldb   = gammalib::strip_whitespace(par->attribute("database"));
            m_xml_rspname = gammalib::strip_whitespace(par->attribute("response"));

            // Set response
            GCaldb caldb;
            if (gammalib::dir_exists(m_xml_caldb)) {
                caldb.rootdir(m_xml_caldb);
            }
            else {
                caldb.open("cta", m_xml_caldb);
            }
            this->caldb(caldb);
            load(m_xml_rspname);
        }

        // Handle effective area
        else if ((par->attribute("name") == "EffectiveArea") ||
                 (par->attribute("name") == "ARF")) {

            // Get filename
            m_xml_aeff = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load effective area
            if (!m_xml_aeff.empty()) {

                // Load effective area
                load_aeff(m_xml_aeff);

                // Optional attributes
                double thetacut = 0.0;
                double scale    = 1.0;
                double sigma    = 0.0;

                // Optionally extract thetacut (0.0 if no thetacut)
                std::string s_thetacut = par->attribute("thetacut");
                if (s_thetacut.length() > 0) {
                    thetacut = gammalib::todouble(s_thetacut);
                }

                // Optionally extract scale factor (1.0 if no scale)
                std::string s_scale = par->attribute("scale");
                if (s_scale.length() > 0) {
                    scale = gammalib::todouble(s_scale);
                }

                // Optionally extract sigma (0.0 if no sigma)
                std::string s_sigma = par->attribute("sigma");
                if (s_sigma.length() > 0) {
                    sigma = gammalib::todouble(s_sigma);
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

        }

        // Handle PSF
        else if ((par->attribute("name") == "PointSpreadFunction") ||
                 (par->attribute("name") == "PSF")) {

            // Get filename
            m_xml_psf = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load point spread function
            if (!m_xml_psf.empty()) {
                load_psf(m_xml_psf);
            }

        }


        // Handle RMF
        else if ((par->attribute("name") == "EnergyDispersion") ||
                 (par->attribute("name") == "RMF")) {

            // Get filename
            m_xml_edisp = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load energy dispersion
            if (!m_xml_edisp.empty()) {
                load_edisp(m_xml_edisp);
            }

        }

        // Handle background model
        else if (par->attribute("name") == "Background") {

            // Get filename
            m_xml_background = gammalib::strip_whitespace(par->attribute("file"));

            // If filename is not empty then load background model
            if (!m_xml_background.empty()) {
                load_background(m_xml_background);
            }

        }

    } // endfor: looped over all parameters

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
    int npars = xml.elements("parameter");

    // If we have a calibration database and response name, then set
    // the information ...
    if (!m_xml_caldb.empty() || !m_xml_rspname.empty()) {
        GXmlElement* par = gammalib::parameter(xml, "Calibration");
        par->attribute("database", m_xml_caldb);
        par->attribute("response", m_xml_rspname);
    }

    // ... otherwise add response components if they exist
    else {

        // Add effective area if it exists
        if (aeff() != NULL) {
            if (!(m_xml_aeff.empty())) {

                // Get pointer to effective area
                GXmlElement* par = gammalib::parameter(xml, "EffectiveArea");

                // Initialise attributes
                double thetacut = 0.0;
                double scale    = 1.0;
                double sigma    = 0.0;

                // Get optional ARF attributes
                const GCTAAeffArf* arf = dynamic_cast<const GCTAAeffArf*>(aeff());
                if (arf != NULL) {
                    thetacut = arf->thetacut();
                    scale    = arf->scale();
                    sigma    = arf->sigma();
                }

                // Get optional performance table attributes
                const GCTAAeffPerfTable* perf = dynamic_cast<const GCTAAeffPerfTable*>(aeff());
                if (perf != NULL) {
                    sigma = perf->sigma();
                }

                // Set attributes
                par->attribute("file", m_xml_aeff);
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
            if (!(m_xml_psf.empty())) {

                // Get pointer to PSD
                GXmlElement* par = gammalib::parameter(xml, "PointSpreadFunction");

                // Write PSF filename
                par->attribute("file", m_xml_psf);

            }
        }

        // Add Edisp if it exists
        if (edisp() != NULL) {
            if (!(m_xml_edisp.empty())) {

                // Get pointer to energy dispersion
                GXmlElement* par = gammalib::parameter(xml, "EnergyDispersion");

                // Write Edisp filename
                par->attribute("file", m_xml_edisp);
                
            }
        }

        // Add background if it exists
        if (background() != NULL) {
            if (!(m_xml_background.empty())) {

                // Get pointer to energy dispersion
                GXmlElement* par = gammalib::parameter(xml, "Background");

                // Write background filename
                par->attribute("file", m_xml_background);

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
 ***************************************************************************/
void GCTAResponseIrf::load(const std::string& rspname)
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


/***********************************************************************//**
 * @brief Load effective area
 *
 * @param[in] filename Effective area filename.
 *
 * This method allocates an effective area instance and load the effective
 * area information from a response file. The following response file formats
 * are supported:
 *
 * (1) A CTA performance table. This is an ASCII file which specifies the
 *     on-axis effective area as function of energy.
 *
 * (2) A ARF FITS file. This is a FITS file which stores the effective area
 *     in a vector.
 *
 * (3) A CTA response table. This is a FITS file which specifies the
 *     effective area as function of energy and offset angle.
 *
 * This method examines the file, and depending on the detected format,
 * allocates the appropriate effective area class and loads the data.
 *
 * First, the method checks whether the file is a FITS file or not. If the
 * file is not a FITS file, it is assumed that the file is an ASCII
 * performance table. If the file is a FITS file, the number of rows found
 * in the table is used to distinguish between an ARF (multiple rows) and
 * a CTA response table (single row).
 *
 * @todo Implement a method that checks if a file is a FITS file instead
 *       of using try-catch.
 ***************************************************************************/
void GCTAResponseIrf::load_aeff(const std::string& filename)
{
    // Free any existing effective area instance
    if (m_aeff != NULL) delete m_aeff;
    m_aeff = NULL;

    // Try opening the file as a FITS file
    try {

        // Open FITS file
        GFits file(filename);

        // If file contains an "EFFECTIVE AREA" extension then load it
        // as CTA response table
        if (file.contains("EFFECTIVE AREA")) {
            file.close();
            m_aeff = new GCTAAeff2D(filename);
        }

        // ... else if file contains a "SPECRESP" extension then load it
        // as ARF
        else if (file.contains("SPECRESP")) {
            file.close();
            m_aeff = new GCTAAeffArf(filename);
        }

    }

    // If FITS file opening failed then assume that we have a performance
    // table
    catch (GException::fits_open_error &e) {
        m_aeff = new GCTAAeffPerfTable(filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load CTA PSF vector
 *
 * @param[in] filename FITS file name.
 *
 * This method loads CTA PSF information from a FITS table. Two FITS file
 * formats are supported by the method:
 *
 * (1) A PSF vector, stored in a format similar to an ARF vector. It is
 * expected that this format is only a preliminary format that will
 * disappear in the future (m_psf_version=-9).
 *
 * (2) A PSF response table, where PSF parameters are given as function of
 * energy, offset angle, and eventually some other parameters. This format
 * is expected to be the definitive response format for CTA
 * (m_psf_version=-8).
 *
 * This method examines the FITS file, and depending on the detected format,
 * calls the relevant methods. Detection is done by the number of rows that
 * are found in the table. A single row means that we deal with a response
 * table, while multiple rows mean that we deal with a response vector.
 *
 * @todo Implement a method that checks if a file is a FITS file instead
 *       of using try-catch.
 ***************************************************************************/
void GCTAResponseIrf::load_psf(const std::string& filename)
{
    // Free any existing point spread function instance
    if (m_psf != NULL) delete m_psf;
    m_psf = NULL;

    // Try opening the file as a FITS file
    try {

        // Open FITS file
        GFits file(filename);

        // If file contains a "POINT SPREAD FUNCTION" extension then load it
        // as either a King profile PSF or a 2D PSF
        if (file.contains("POINT SPREAD FUNCTION")) {
            const GFitsTable& table = *file.table("POINT SPREAD FUNCTION");
            if (table.contains("GAMMA") && table.contains("SIGMA")) {
                file.close();
                m_psf = new GCTAPsfKing(filename);
            }
            else if (table.contains("SCALE") && table.contains("SIGMA_1") &&
                     table.contains("AMPL_2") && table.contains("SIGMA_2") &&
                     table.contains("AMPL_3") && table.contains("SIGMA_3")) {
                file.close();
                m_psf = new GCTAPsf2D(filename);
            }
            else {
                file.close();
            }
        }

        // ... else load it has PSF vector 
        else {
            file.close();
            m_psf = new GCTAPsfVector(filename);
        }

    }

    // If FITS file opening failed then assume that we have a performance
    // table
    catch (GException::fits_open_error &e) {
        m_psf = new GCTAPsfPerfTable(filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy dispersion information
 *
 * @param[in] filename Energy dispersion file name.
 ***************************************************************************/
void GCTAResponseIrf::load_edisp(const std::string& filename)
{
    // Free any existing energy dispersion instance
    if (m_edisp != NULL) delete m_edisp;
    m_edisp = NULL;

    // Try opening the file as a FITS file
    try {

        // Open FITS file
        GFits file(filename);

        // If file contains an "ENERGY DISPERSION" extension then load it
        // as CTA response table
        if (file.contains("ENERGY DISPERSION")) {
            file.close();
            //m_edisp = new GCTAEdisp2D(filename);
        }

        // ... else load it as RMF
        else {
            file.close();
            m_edisp = new GCTAEdispRmf(filename);
        }

    }

    // If FITS file opening failed then assume that we have a performance
    // table
    catch (GException::fits_open_error &e) {
        m_edisp = new GCTAEdispPerfTable(filename);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load background model
 *
 * @param[in] filename Background model file name.
 ***************************************************************************/
void GCTAResponseIrf::load_background(const std::string& filename)
{
    // Free any existing background model instance
    if (m_background != NULL) delete m_background;
    m_background = NULL;

    // Try opening the file as a FITS file
    try {
        // Load background as 3D background
        m_background = new GCTABackground3D(filename);
    }
    catch (GException::fits_open_error &e) {
        // Load background as performance table background
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
 * @param[in] chatter Chattiness (defaults to NORMAL).
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

        // Append calibration database
        result.append("\n"+m_caldb.print(chatter));

        // Append effective area information
        if (m_aeff != NULL) {
            result.append("\n"+m_aeff->print(chatter));
        }

        // Append point spread function information
        if (m_psf != NULL) {
            result.append("\n"+m_psf->print(chatter));
        }

        // Append energy dispersion information
        if (m_edisp != NULL) {
            result.append("\n"+m_edisp->print(chatter));
        }

        // Append background information
        if (m_background != NULL) {
            result.append("\n"+m_background->print(chatter));
        }

        // EXPLICIT: Append Npred cache information
        if (chatter >= EXPLICIT) {
            if (!m_npred_names.empty()) {
                for (int i = 0; i < m_npred_names.size(); ++i) {
                    result.append("\n"+gammalib::parformat("Npred cache " +
                                  gammalib::str(i)));
                    result.append(m_npred_names[i]+", ");
                    result.append(m_npred_energies[i].print()+", ");
                    result.append(m_npred_times[i].print()+" = ");
                    result.append(gammalib::str(m_npred_values[i]));
                }
            }
        } // endif: chatter was explicit

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =              Model type dependent CTA response methods                  =
 =                                                                         =
 ==========================================================================*/

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
 * Integrates the product of the model and the IRF over the true photon
 * arrival direction using
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *    \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *    IRF(\rho, \omega) d\omega d\rho
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho | E, t)\f$ is the radial model,
 * - \f$IRF(\rho, \omega)\f$ is the IRF
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
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
 * direction). Given the slow variation of the PSF shape over the field of
 * view, this approximation should be fine. It helps in fact a lot in
 * speeding up the computations.
 ***************************************************************************/
double GCTAResponseIrf::irf_radial(const GEvent&       event,
                                   const GSource&      source,
                                   const GObservation& obs) const
{
    // Retrieve CTA pointing
    const GCTAPointing& pnt = retrieve_pnt(G_IRF_RADIAL, obs);
    const GCTAInstDir&  dir = retrieve_dir(G_IRF_RADIAL, event);

    // Get pointer on radial model
    const GModelSpatialRadial* model =
          dynamic_cast<const GModelSpatialRadial*>(source.model());
    if (model == NULL) {
        throw GCTAException::bad_model_type(G_IRF_RADIAL);
    }

    // Get event attributes
    //const GSkyDir& obsDir = dir->dir();
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
    double zeta = centre.dist(dir.dir());

    // Determine angular distance between measured photon direction and
    // pointing direction [radians]
    double eta = pnt.dir().dist(dir.dir());

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
        cta_irf_radial_kern_rho integrand(*this,
                                          *model,
                                          zenith,
                                          azimuth,
                                          srcEng,
                                          srcTime,
                                          srcLogEng,
                                          obsEng,
                                          zeta,
                                          lambda,
                                          omega0,
                                          delta_max);

        // Integrate over zenith angle
        GIntegral integral(&integrand);
        integral.eps(1.0e-5);
        irf = integral.romb(rho_min, rho_max);

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

    }

    // Compile option: Show integration results
    #if defined(G_DEBUG_IRF_RADIAL)
    std::cout << "GCTAResponseIrf::irf_radial:";
    std::cout << " rho_min=" << rho_min;
    std::cout << " rho_max=" << rho_max;
    std::cout << " irf=" << irf << std::endl;
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return IRF value for elliptical source model
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
 *    \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *    S_{\rm p}(\rho, \omega | E, t) \, IRF(\rho, \omega) d\omega d\rho
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the radial model,
 * - \f$IRF(\rho, \omega)\f$ is the IRF
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
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
    // Retrieve CTA pointing
    const GCTAPointing& pnt = retrieve_pnt(G_IRF_ELLIPTICAL, obs);
    const GCTAInstDir&  dir = retrieve_dir(G_IRF_ELLIPTICAL, event);

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
    double zeta     = centre.dist(obsDir);
    double obsOmega = centre.posang(obsDir);

    // Determine angular distance between measured photon direction and
    // pointing direction [radians]
    double eta = pnt.dir().dist(obsDir);

    // Determine angular distance between model centre and pointing direction
    // [radians]
    double lambda = centre.dist(pnt.dir());

    // Compute azimuth angle of pointing in model coordinate system [radians]
    // This azimuth angle is comprised in the interval [0,pi], and defines
    // the zero point of the model coordinate system.
    double omega0 = 0.0;
    double denom  = std::sin(lambda) * std::sin(zeta);
    if (denom != 0.0) {
        double arg = (std::cos(eta) - std::cos(lambda) * std::cos(zeta))/denom;
        omega0     = gammalib::acos(arg);
    }

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius [radians]. We assign here the measured theta
    // angle (eta) as the true theta angle between the source and the pointing
    // directions. As we only use the angle to determine the maximum PSF size,
    // this should be sufficient.
    double theta     = eta;
    double phi       = 0.0; //TODO: Implement IRF Phi dependence
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

    // Get maximum source model radius [radians]
    double src_max = model->theta_max();

    // Set zenith angle integration range for elliptical model
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
        cta_irf_elliptical_kern_rho integrand(*this,
                                              *model,
                                              zenith,
                                              azimuth,
                                              srcEng,
                                              srcTime,
                                              srcLogEng,
                                              obsEng,
                                              zeta,
                                              lambda,
                                              obsOmega,
                                              omega0,
                                              delta_max);

        // Integrate over zenith angle
        GIntegral integral(&integrand);
        integral.eps(1.0e-5);
        irf = integral.romb(rho_min, rho_max);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GCTAResponseIrf::irf_elliptical:";
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
    }

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
    // Initialise IRF value
    bool   has_irf = false;
    double irf     = 0.0;

    // Retrieve CTA pointing
    const GCTAPointing& pnt = retrieve_pnt(G_IRF_DIFFUSE, obs);

    // Try getting the IRF value from cache
    #if defined(G_USE_IRF_CACHE)
    const GCTAEventList* list = dynamic_cast<const GCTAEventList*>(obs.events());
    const GCTAEventAtom* atom = dynamic_cast<const GCTAEventAtom*>(&event);
    if (list != NULL && atom != NULL) {
        irf = list->irf_cache(source.name(), atom->index());
        if (irf >= 0.0) {
            has_irf = true;
            #if defined(G_DEBUG_IRF_DIFFUSE)
            std::cout << "GCTAResponseIrf::irf_diffuse:";
            std::cout << " cached irf=" << irf << std::endl;
            #endif
        }
        else {
            irf = 0.0;
        }
    }
    #endif

    // Continue only if we have no IRF value
    if (!has_irf) {

        // Get CTA instrument direction
        const GCTAInstDir&  dir = retrieve_dir(G_IRF_ELLIPTICAL, event);

        // Get pointer on spatial model
        const GModelSpatial* model =
            dynamic_cast<const GModelSpatial*>(source.model());
        if (model == NULL) {
            throw GCTAException::bad_model_type(G_IRF_DIFFUSE);
        }

        // Get event attributes
        //const GSkyDir& obsDir = dir.dir();
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
            cta_irf_diffuse_kern_theta integrand(*this,
                                                 *model,
                                                 theta,
                                                 phi,
                                                 zenith,
                                                 azimuth,
                                                 srcEng,
                                                 srcTime,
                                                 srcLogEng,
                                                 obsEng,
                                                 rot,
                                                 eta);

            // Integrate over zenith angle
            GIntegral integral(&integrand);
            integral.eps(1.0e-4);
            irf = integral.romb(0.0, delta_max);

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

        // Put IRF value in cache
        #if defined(G_USE_IRF_CACHE)
        if (list != NULL && atom != NULL) {
            list->irf_cache(source.name(), atom->index(), irf);
        }
        #endif

        // Compile option: Show integration results
        #if defined(G_DEBUG_IRF_DIFFUSE)
        std::cout << "GCTAResponseIrf::irf_diffuse:";
        std::cout << " srcLogEng=" << srcLogEng;
        std::cout << " obsLogEng=" << obsLogEng;
        std::cout << " eta=" << eta;
        std::cout << " delta_max=" << delta_max;
        std::cout << " irf=" << irf << std::endl;
        #endif

    } // endif: has no IRF

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return spatial integral of radial source model over ROI
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_model_type
 *            Model is not a radial model.
 *
 * Integrates the product of the radial model and Npred over the Region Of
 * Interest using
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *    \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *    N_{\rm pred}(\rho, \omega) d\omega
 *    d\rho
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho | E, t)\f$ is the radial model,
 * - \f$N_{\rm pred}(\rho,\omega)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest,
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 *
 * The integration is performed in a spherical coordinate system that is
 * centred on the source model centre \f$\vec{m}\f$, with \f$(\rho,\omega)\f$
 * being the zenith and azimuth angles, respectively.
 *
 * The zenith angle integration range \f$[\rho_{\rm min}, \rho_{\rm max}\f$
 * and azimuth angle integration range 
 * \f$[\omega_{\rm min}, \omega_{\rm max}\f$
 * are adjusted so that only coordinates within the circular region of
 * interest will be considered.
 *
 * Note that we estimate the integration radius based on the size of the
 * onaxis PSF. This should be fine as long as the offaxis PSF is not
 * considerably larger than the onaxis PSF. We should verify this, however.
 *
 * @todo Verify that offaxis PSF is not considerably larger than onaxis
 *       PSF. 
 ***************************************************************************/
double GCTAResponseIrf::npred_radial(const GSource& source,
                                     const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = retrieve_obs(G_NPRED_RADIAL, obs);
    const GCTARoi&         roi = retrieve_roi(G_NPRED_RADIAL, obs);
    const GCTAPointing&    pnt = cta.pointing();

    // Get pointer on radial model
    const GModelSpatialRadial* model =
          dynamic_cast<const GModelSpatialRadial*>(source.model());
    if (model == NULL) {
        throw GCTAException::bad_model_type(G_NPRED_RADIAL);
    }

    // Get source attributes
    const GSkyDir& centre  = model->dir();
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

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
    double rho_max = model->theta_max();

    // Perform offset angle integration only if interval is valid
    if (rho_max > rho_min) {

        // Compute rotation matrix to convert from native model coordinates,
        // given by (rho,omega), into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        ry.eulery(model->dec() - 90.0);
        rz.eulerz(-model->ra());
        GMatrix rot = (ry * rz).transpose();

        // Compute position angle of ROI centre with respect to model
        // centre (radians)
        double omega0 = centre.posang(roi.centre().dir());

        // Setup integration kernel
        cta_npred_radial_kern_rho integrand(*this,
                                            *model,
                                            source.energy(),
                                            source.time(),
                                            cta,
                                            rot,
                                            roi_model_distance,
                                            roi_psf_radius,
                                            omega0);

        // Integrate over theta
        GIntegral integral(&integrand);
        integral.eps(1.0e-5);
        npred = integral.romb(rho_min, rho_max);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NPRED_RADIAL)
        std::cout << "GCTAResponseIrf::npred_radial:";
        std::cout << " rho_min=" << rho_min;
        std::cout << " rho_max=" << rho_max;
        std::cout << " npred=" << npred << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GCTAResponseIrf::npred_radial:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", rho_min=" << rho_min;
        std::cout << ", rho_max=" << rho_max;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return spatial integral of elliptical source model over ROI
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_model_type
 *            Model is not an elliptical model.
 *
 * Integrates the product of the elliptical model and Npred over the Region
 * Of Interest using
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}}
 *    \sin \rho \times
 *    \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *    S_{\rm p}(\rho,\omega | E, t) \,
 *    N_{\rm pred}(\rho,\omega) d\omega
 *    d\rho
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho,\omega | E, t)\f$ is the elliptical model,
 * - \f$N_{\rm pred}(\rho,\omega)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest,
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 *
 * The integration is performed in a spherical coordinate system that is
 * centred on the source model centre \f$\vec{m}\f$, with \f$(\rho,\omega)\f$
 * being the zenith and azimuth angles, respectively.
 *
 * The zenith angle integration range \f$[\rho_{\rm min}, \rho_{\rm max}\f$
 * and azimuth angle integration range 
 * \f$[\omega_{\rm min}, \omega_{\rm max}\f$
 * are adjusted so that only coordinates within the circular region of
 * interest will be considered.
 *
 * Note that we estimate the integration radius based on the size of the
 * onaxis PSF. This should be fine as long as the offaxis PSF is not
 * considerably larger than the onaxis PSF. We should verify this, however.
 *
 * @todo Verify that offaxis PSF is not considerably larger than onaxis
 *       PSF.
 ***************************************************************************/
double GCTAResponseIrf::npred_elliptical(const GSource& source,
                                         const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Retrieve CTA observation, ROI and pointing
    const GCTAObservation& cta = retrieve_obs(G_NPRED_ELLIPTICAL, obs);
    const GCTARoi&         roi = retrieve_roi(G_NPRED_ELLIPTICAL, obs);
    const GCTAPointing&    pnt = cta.pointing();

    // Get pointer on elliptical model
    const GModelSpatialElliptical* model =
          dynamic_cast<const GModelSpatialElliptical*>(source.model());
    if (model == NULL) {
        throw GCTAException::bad_model_type(G_NPRED_ELLIPTICAL);
    }

    // Get source attributes
    const GSkyDir& centre  = model->dir();
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

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
    double rho_max = model->theta_max();

    // Perform offset angle integration only if interval is valid
    if (rho_max > rho_min) {

        // Compute rotation matrix to convert from native model coordinates,
        // given by (rho,omega), into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        ry.eulery(model->dec() - 90.0);
        rz.eulerz(-model->ra());
        GMatrix rot = (ry * rz).transpose();

        // Compute position angle of ROI centre with respect to model
        // centre (radians)
        double omega0 = centre.posang(roi.centre().dir());

        // Setup integration kernel
        cta_npred_elliptical_kern_rho integrand(*this,
                                                *model,
                                                source.energy(),
                                                source.time(),
                                                cta,
                                                rot,
                                                roi_model_distance,
                                                roi_psf_radius,
                                                omega0);

        // Integrate over theta
        GIntegral integral(&integrand);
        integral.eps(1.0e-5);
        npred = integral.romb(rho_min, rho_max);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NPRED_ELLIPTICAL)
        std::cout << "GCTAResponseIrf::npred_elliptical:";
        std::cout << " rho_min=" << rho_min;
        std::cout << " rho_max=" << rho_max;
        std::cout << " npred=" << npred << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GCTAResponseIrf::npred_elliptical:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", rho_min=" << rho_min;
        std::cout << ", rho_max=" << rho_max;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return spatial integral of diffuse source model over ROI
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_model_type
 *            Model is not a radial model.
 *
 * Integrates the product of the diffuse model and Npred over the Region
 * Of Interest using
 *
 * \f[
 *    \int_{0}^{\theta_{\rm max}}
 *    \sin \theta \times
 *    \int_{0}^{2\pi}
 *    S_{\rm p}(\theta, \phi | E, t) \,
 *    N_{\rm pred}(\theta, \phi) d\phi
 *    d\theta
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\theta, \phi | E, t)\f$ is the diffuse model,
 * - \f$N_{\rm pred}(\theta, \phi)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest in the reference frame of the diffuse source
 *   model
 * - \f$\theta\f$ is the distance from the ROI centre, and
 * - \f$\phi\f$ is the azimuth angle.
 *
 * Note that the integration precision was adjusted trading-off between
 * computation time and computation precision. A value of 1e-4 was judged
 * appropriate.
 *
 * Note that we estimate the integration radius based on the size of the
 * onaxis PSF in this method. This should be fine as long as the offaxis
 * PSF is not considerably larger than the onaxis PSF. We should verify
 * this, however.
 *
 * @todo Verify that offaxis PSF is not considerably larger than onaxis PSF.
 ***************************************************************************/
double GCTAResponseIrf::npred_diffuse(const GSource& source,
                                      const GObservation& obs) const
{
    // Initialise Npred value
    double npred     = 0.0;
    bool   has_npred = false;

    // Build unique identifier
    std::string id = source.name() + "::" + obs.id();

    // Check if Npred value is already in cache
    #if defined(G_USE_NPRED_CACHE)
    if (!m_npred_names.empty()) {

         // Search for unique identifier, and if found, recover Npred value
         // and break
         for (int i = 0; i < m_npred_names.size(); ++i) {
             if (m_npred_names[i]    == id &&
                 m_npred_energies[i] == source.energy() &&
                 m_npred_times[i]    == source.time()) {
                 npred = m_npred_values[i];
                 has_npred = true;
                 #if defined(G_DEBUG_NPRED_DIFFUSE)
                 std::cout << "GCTAResponseIrf::npred_diffuse:";
                 std::cout << " cache=" << i;
                 std::cout << " npred=" << npred << std::endl;
                 #endif
                 break;
             }
         }

    } // endif: there were values in the Npred cache
    #endif

    // Continue only if no Npred cache value was found
    if (!has_npred) {

        // Retrieve CTA observation, ROI and pointing
        const GCTAObservation& cta = retrieve_obs(G_NPRED_DIFFUSE, obs);
        const GCTARoi&         roi = retrieve_roi(G_NPRED_DIFFUSE, obs);
        const GCTAPointing&    pnt = cta.pointing();

        // Get pointer on spatial model
        const GModelSpatial* model =
            dynamic_cast<const GModelSpatial*>(source.model());
        if (model == NULL) {
            throw GCTAException::bad_model_type(G_NPRED_DIFFUSE);
        }

        // Get source attributes
        const GEnergy& srcEng  = source.energy();
        const GTime&   srcTime = source.time();

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
            cta_npred_diffuse_kern_theta integrand(*this,
                                                   *model,
                                                   source.energy(),
                                                   source.time(),
                                                   cta,
                                                   rot);

            // Integrate over theta
            GIntegral integral(&integrand);
            integral.eps(1.0e-5);
            npred = integral.romb(0.0, roi_psf_radius);

            // Compile option: Show integration results
            #if defined(G_DEBUG_NPRED_DIFFUSE)
            std::cout << "GCTAResponseIrf::npred_diffuse:";
            std::cout << " roi_psf_radius=" << roi_psf_radius;
            std::cout << " npred=" << npred;
            std::cout << " id=" << id << std::endl;
            #endif

        } // endif: offset angle range was valid

        // Store result in Npred cache
        #if defined(G_USE_NPRED_CACHE)
        m_npred_names.push_back(id);
        m_npred_energies.push_back(source.energy());
        m_npred_times.push_back(source.time());
        m_npred_values.push_back(npred);
        #endif

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
            std::cout << "*** ERROR: GCTAResponseIrf::npred_diffuse:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", roi_psf_radius=" << roi_psf_radius;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: Npred computation required

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return true energy boundaries for a specific observed energy
 *
 * @param[in] obsEnergy Observed Energy.
 * @return True energy boundaries for given observed energy.
 *
 * @todo So far we have no means to pass additional parameters to the
 * GCTAEdisp::ebounds() method.
 ***************************************************************************/
GEbounds GCTAResponseIrf::ebounds_src(const GEnergy& obsEnergy) const
{
    // Initialise an empty boundary object
    GEbounds ebounds;

    // If energy dispersion is available then set the energy boundaries
    if (edisp() != NULL) {
        double obsLogEng = obsEnergy.log10TeV();
        ebounds          = edisp()->ebounds_src(obsLogEng); // Requires TeV
    }

    // Return energy boundaries
    return ebounds;
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
 * @brief Return energy dispersion (in units or MeV^-1)
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 ***************************************************************************/
double GCTAResponseIrf::edisp(const GEnergy& obsEng,
                              const double&  theta,
                              const double&  phi,
                              const double&  zenith,
                              const double&  azimuth,
                              const double&  srcLogEng) const
{
    // Throw an exception if instrument response is not defined
    if (m_edisp == NULL) {
        std::string msg = "No energy dispersion information found in"
                          " response.\n"
                          "Please make sure that the instrument response is"
                          " properly defined.";
        throw GException::invalid_value(G_EDISP, msg);
    }

    // Compute log10 energy in TeV and linear energy in MeV
    double obsLogEng = obsEng.log10TeV();
    double energy    = obsEng.MeV();

    // Compute energy dispersion
    double edisp = (*m_edisp)(obsLogEng, srcLogEng, theta, phi, zenith, azimuth) /
                   (gammalib::ln10 * energy);

    // Return energy dispersion
    return edisp;
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
 * @todo Enhance romb() integration method for small integration regions
 *       (see comment about kluge below)
 * @todo Implement phi dependence in camera system
 ***************************************************************************/
double GCTAResponseIrf::npsf(const GSkyDir&      srcDir,
                             const double&       srcLogEng,
                             const GTime&        srcTime,
                             const GCTAPointing& pnt,
                             const GCTARoi&      roi) const
{
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
            cta_npsf_kern_rad_azsym integrand(*this,
                                              roi_radius,
                                              roi_psf_distance,
                                              srcLogEng,
                                              theta,
                                              phi,
                                              zenith,
                                              azimuth);

            // Setup integration
            GIntegral integral(&integrand);
            integral.eps(1.0e-5);

            // Radially integrate PSF. In case that the radial integration
            // region is small, we do the integration using a simple
            // trapezoidal rule. This is a kluge to prevent convergence
            // problems in the romb() method for small integration intervals.
            // Ideally, the romb() method should be enhanced to handle this
            // case automatically. The kluge threshold was fixed manually!
            if (rmax-rmin < 1.0e-12) {
                value = integral.trapzd(rmin, rmax);
            }
            else {
                value = integral.romb(rmin, rmax);
            }

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
                std::cout << "*** ERROR: GCTAResponseIrf::npsf:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (value=" << value;
                std::cout << ", roi_radius=" << roi_radius;
                std::cout << ", roi_psf_distance=" << roi_psf_distance;
                //std::cout << ", sigma=" << sigma;
                std::cout << ", r=[" << rmin << "," << rmax << "])";
                std::cout << std::endl;
            }
            #endif

        } // endif: integration range was valid

    } // endelse: numerical integration required

    // Return integrated PSF
    return value;
}


/***********************************************************************//**
 * @brief Return result of energy dispersion integral over energy range
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt CTA pointing.
 * @param[in] ebds Energy boundaries of data selection.
 *
 * @todo Implement phi dependence in camera system
 ***************************************************************************/
double GCTAResponseIrf::nedisp(const GSkyDir&      srcDir,
                               const GEnergy&      srcEng,
                               const GTime&        srcTime,
                               const GCTAPointing& pnt,
                               const GEbounds&     ebds) const
{
    // Initialise energy dispersion integral
    double nedisp = 1.0;

    // Continue only if energy dispersion information is available
    if (edisp() != NULL) {

        // Get the observed energy boundaries for specified true energy
        GEbounds ebounds = edisp()->ebounds_obs(srcEng.log10TeV());

        // Check if at least one of the energy boundaries covered by the
        // energy dispersion for the specified true energy lies outside
        // any of the energy boundaries of the data selection (or lies
        // within different energy boundaries of the data selection)
        bool outside = false;
        for (int k = 0; k < ebounds.size(); ++k) {
            int imin = ebds.index(ebounds.emin(k));
            int imax = ebds.index(ebounds.emax(k));
            if (imin != imax || imin == -1 || imax == -1) {
                outside = true;
                break;
            }
        }

        // If energy boundaries are not fully covered then integrate
        // numerically
        if (outside) {

            // Initialise energy dispersion integral
            nedisp = 0.0;

            // Get pointing direction zenith angle and azimuth [radians]
            double zenith  = pnt.zenith();
            double azimuth = pnt.azimuth();

            // Compute offset angle of source direction in camera system
            double theta = pnt.dir().dist(srcDir);

            // Compute azimuth angle of source direction in camera system
            double phi = 0.0; //TODO: Implement phi dependence

            // Loop over energy boundaries in observed energy
            for (int i = 0; i < ebds.size(); ++i) {

                // Get boundaries in observed energy
                GEnergy emin_obs = ebds.emin(i);
                GEnergy emax_obs = ebds.emax(i);

                // Loop over energy boundaries of energy dispersion
                for (int k = 0; k < ebounds.size(); ++k) {

                    // Get boundaries of energy dispersion
                    GEnergy emin_edisp = ebounds.emin(k);
                    GEnergy emax_edisp = ebounds.emax(k);

                    // Get energy dispersion interval that overlaps with
                    // the observed energy interval
                    GEnergy emin = (emin_edisp < emin_obs) ? emin_obs : emin_edisp;
                    GEnergy emax = (emax_edisp > emax_obs) ? emax_obs : emax_edisp;

                    // If interval has positive length then integrate over
                    // the energy dispersion
                    if (emin < emax) {

                        // Get log10 of energy boundaries in TeV
                        double e_log_min = emin.log10TeV();
                        double e_log_max = emax.log10TeV();

                        // Setup integration function
                        cta_nedisp_kern integrand(*this,
                                                  srcEng.log10TeV(),
                                                  theta,
                                                  phi,
                                                  zenith,
                                                  azimuth);
                        GIntegral integral(&integrand);

                        // Set integration precision
                        integral.eps(1.0e-3);

                        // Do Romberg integration
                        nedisp += integral.romb(e_log_min, e_log_max);

                    } // endif: integration range was valid
                
                } // endfor: looped over energy boundaries of energy dispersion
            } // endfor: looped over energy boundaries in observed energy

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(nedisp) || gammalib::is_infinite(nedisp)) {
                std::cout << "*** ERROR: GCTAResponseIrf::nedisp:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (nedisp=" << nedisp;
                std::cout << ", srcEng=" << srcEng;
                std::cout << ", srcTime=" << srcTime;
                std::cout << ", pnt=" << pnt;
                std::cout << ", ebds=" << ebds;
                std::cout << ")" << std::endl;
            }
            #endif

        } // endif: numerical integration was needed

    } // endif: there is an energy dispersion response

    // Return integral
    return nedisp;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * We set the relative integration precision to 1e-5 as test images are
 * pretty smooth with this precision and computations are still reasonably
 * fast.
 *
 * The following timing was obtained on a 64 Bit machine (fermi) using the
 * script ./test_model for a disk, a Gaussian, and a shell model:
 *
 *                           Disk      Gauss      Shell
 *      m_eps = 1e-3 : user 0m03.80s  0m13.41s   0m03.83s
 *      m_eps = 1e-4 : user 0m03.85s  0m13.71s   0m04.68s
 *      m_eps = 1e-5 : user 0m06.29s  0m23.22s   0m16.94s
 *      m_eps = 1e-6 : user 0m12.65s  0m55.08s   1m32.52s
 ***************************************************************************/
void GCTAResponseIrf::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();
    m_eps         = 1.0e-5; //!< Precision for Romberg integration
    m_aeff        = NULL;
    m_psf         = NULL;
    m_edisp       = NULL;
    m_background  = NULL;
    m_apply_edisp = false;  //!< Switched off by default

    // XML response filenames
    m_xml_caldb.clear();
    m_xml_rspname.clear();
    m_xml_aeff.clear();
    m_xml_psf.clear();
    m_xml_edisp.clear();
    m_xml_background.clear();

    // Initialise Npred cache
    m_npred_names.clear();
    m_npred_energies.clear();
    m_npred_times.clear();
    m_npred_values.clear();

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
    m_caldb       = rsp.m_caldb;
    m_rspname     = rsp.m_rspname;
    m_eps         = rsp.m_eps;
    m_apply_edisp = rsp.m_apply_edisp;

    // Copy response filenames
    m_xml_caldb      = rsp.m_xml_caldb;
    m_xml_rspname    = rsp.m_xml_rspname;
    m_xml_aeff       = rsp.m_xml_aeff;
    m_xml_psf        = rsp.m_xml_psf;
    m_xml_edisp      = rsp.m_xml_edisp;
    m_xml_background = rsp.m_xml_background;

    // Copy cache
    m_npred_names    = rsp.m_npred_names;
    m_npred_energies = rsp.m_npred_energies;
    m_npred_times    = rsp.m_npred_times;
    m_npred_values   = rsp.m_npred_values;

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
const GCTAObservation& GCTAResponseIrf::retrieve_obs(const std::string& origin,
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
 * @brief Retrieve CTA pointing from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 *
 * Extract CTA pointing from a CTA observation.
 ***************************************************************************/
const GCTAPointing& GCTAResponseIrf::retrieve_pnt(const std::string& origin,
                                                  const GObservation& obs) const
{
    // Retrieve CTA observation and pointing
    const GCTAObservation& cta = retrieve_obs(origin, obs);

    // Return CTA pointing
    return (cta.pointing());
}


/***********************************************************************//**
 * @brief Retrieve CTA ROI from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA event list.
 *
 * Extract CTA Region of Interest from a CTA observation.
 ***************************************************************************/
const GCTARoi& GCTAResponseIrf::retrieve_roi(const std::string& origin,
                                             const GObservation& obs) const
{
    // Retrieve CTA observation
    const GCTAObservation& cta = retrieve_obs(origin, obs);

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(cta.events());

    // If pointer is not valid then throw an exception
    if (events == NULL) {
        std::string msg = "Specified observation does not contain a CTA event"
                          " list.\n"
                          "Please specify a CTA observation containing a CTA"
                          " event list when calling this method.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA ROI
    return (events->roi());
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
const GCTAInstDir& GCTAResponseIrf::retrieve_dir(const std::string& origin,
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
std::string GCTAResponseIrf::irf_filename(const std::string& filename) const
{
    // Set input filename as result filename
    std::string result = filename;

    // If file does not exist then try a variant with extension .dat
    if (!gammalib::file_exists(result)) {
        std::string testname = result + ".dat";
        if (gammalib::file_exists(testname)) {
            result = testname;
        }
    }

    // Return result
    return result;
}
