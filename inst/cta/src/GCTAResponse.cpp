/***************************************************************************
 *                  GCTAResponse.cpp  -  CTA Response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Juergen Knoedlseder                         *
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
 * @file GCTAResponse.cpp
 * @brief CTA response class implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>           // access() function
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include <math.h>
#include <vector>
#include <string>
#include "GModelSpatialPtsrc.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"
#include "GVector.hpp"
#include "GCaldb.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GCTAPointing.hpp"
#include "GCTAEventList.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"
#include "GCTADir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                           "GCTAResponse::caldb(std::string&)"
#define G_IRF      "GCTAResponse::irf(GInstDir&, GEnergy&, GTime&, GSkyDir&,"\
                                          " GEnergy&, GTime&, GObservation&)"
#define G_NPRED             "GCTAResponse::npred(GSkyDir&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC            "GCTAResponse::mc(double&,GPhoton&,GPointing&,GRan&)"
#define G_IRF_EXTENDED      "GCTAResponse::irf_extended(GInstDir&, GEnergy&,"\
           " GTime&, GModelExtendedSource&, GEnergy&, GTime&, GObservation&)"
#define G_IRF_DIFFUSE     "GCTAResponse::irf_diffuse(GCTAInstDir&, GEnergy&,"\
         " GTime&, GModelDiffuseSource&, GEnergy&, GTime&, GCTAObservation&)"
#define G_NPRED_EXTENDED                      "GCTAResponse::npred_extended(" \
                       "GModelExtendedSource&,GEnergy&,GTime&,GObservation&)"
#define G_READ           "GCTAResponse::read_performance_table(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_IRF_EXTENDED                 //!< Debug irf_extended method
//#define G_DEBUG_NPRED_EXTENDED                    //!< Debug npred_extended

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAResponse::GCTAResponse(void) : GResponse()
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
 **************************************************************************/
GCTAResponse::GCTAResponse(const GCTAResponse& rsp) : GResponse(rsp)
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
 * @param[in] caldb Calibration database path (defaults to "").
 *
 * Create instance by specifying the response file name and the calibration
 * database path. If an empty string is passed as calibration database path
 * the method will use the CALDB environment variable to determine
 * calibration database path.
 ***************************************************************************/
GCTAResponse::GCTAResponse(const std::string& rspname,
                           const std::string& caldb) : GResponse()
{
    // Initialise members
    init_members();

    // Set calibration database
    this->caldb(caldb);

    // Load IRF
    this->load(rspname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAResponse::~GCTAResponse(void)
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
 ***************************************************************************/
GCTAResponse& GCTAResponse::operator=(const GCTAResponse& rsp)
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
***************************************************************************/
void GCTAResponse::clear(void)
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
***************************************************************************/
GCTAResponse* GCTAResponse::clone(void) const
{
    return new GCTAResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon arrival direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_observation_type
 *            Observation is not a CTA observations.
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found.
 * @exception GCTAException::bad_instdir_type
 *            Instrument direction is not a valid CTA instrument direction.
 *
 * @todo Implement Phi dependence in CTA IRF
 ***************************************************************************/
double GCTAResponse::irf(const GInstDir&     obsDir,
                         const GEnergy&      obsEng,
                         const GTime&        obsTime,
                         const GSkyDir&      srcDir,
                         const GEnergy&      srcEng,
                         const GTime&        srcTime,
                         const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL)
        throw GCTAException::bad_observation_type(G_IRF);

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing(srcTime);
    if (pnt == NULL)
        throw GCTAException::no_pointing(G_IRF);

    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&obsDir);
    if (dir == NULL)
        throw GCTAException::bad_instdir_type(G_IRF);

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt->zenith();
    double azimuth = pnt->azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta = pnt->dir().dist(srcDir);
    double phi   = 0.0; //TODO: Implement Phi dependence

    // Get log10(E/TeV) of true photon energy.
    double srcLogEng = srcEng.log10TeV();

    // Determine angular separation between true and measured photon
    // direction in radians
    double delta = dir->dist(srcDir);

    // Get maximum angular separation for which PSF is significant
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

    // Initialise IRF value
    double irf = 0.0;

    // Compute only if we're sufficiently close to PSF
    if (delta <= delta_max) {

        // Get effective area component
        irf = aeff(theta, phi, zenith, azimuth, srcLogEng);
        
        // Multiply-in PSF
        if (irf > 0) {
        
            // Get PSF component
            irf *= psf(delta, theta, phi, zenith, azimuth, srcLogEng);
            
            // Multiply-in energy dispersion
            if (hasedisp() && irf > 0) {

                // Get log10(E/TeV) of measured photon energy.
                double obsLogEng = obsEng.log10TeV();

                // Multiply-in energy dispersion
                irf *= edisp(obsLogEng, theta, phi, zenith, azimuth, srcLogEng);

            } // endif: energy dispersion was available and PSF was non-zero
            
        } // endif: Aeff was non-zero

    } // endif: we were sufficiently close to PSF

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(irf) || isinfinite(irf)) {
        std::cout << "*** ERROR: GCTAResponse::irf:";
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
 * @param[in] srcDir True photon arrival direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_observation_type
 *            Observation is not a CTA observations.
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found.
 * @exception GException::no_list
 *            Observation does not contain a valid CTA event list.
 *
 * @todo Set polar angle of photon in camera system
 * @todo Set telescope zenith and azimuth angles
 *
 * @todo Implement Phi dependence in CTA IRF
 ***************************************************************************/
double GCTAResponse::npred(const GSkyDir&      srcDir,
                           const GEnergy&      srcEng,
                           const GTime&        srcTime,
                           const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL)
        throw GCTAException::bad_observation_type(G_NPRED);

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing(srcTime);
    if (pnt == NULL)
        throw GCTAException::no_pointing(G_NPRED);

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(ctaobs->events());
    if (events == NULL)
        throw GException::no_list(G_NPRED);

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt->zenith();
    double azimuth = pnt->azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta = pnt->dir().dist(srcDir);
    double phi   = 0.0; //TODO: Implement Phi dependence

    // Get log10(E/TeV) of true photon energy.
    double srcLogEng = srcEng.log10TeV();

    // Get effectve area components
    double npred = aeff(theta, phi, zenith, azimuth, srcLogEng);
    
    // Multiply-in PSF
    if (npred > 0.0) {
    
        // Get PSF
        npred *= npsf(srcDir, srcLogEng, srcTime, *pnt, events->roi());

        // Multiply-in energy dispersion
        if (hasedisp() && npred > 0.0) {

            // Get energy dispersion
            npred *= nedisp(srcDir, srcEng, srcTime, *pnt, events->ebounds());

        } // endif: had energy dispersion
    
    } // endif: had non-zero effective area

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GCTAResponse::npred:";
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
 * @param[in] pnt Pointing.
 * @param[in] ran Random number generator.
 *
 * @exception GCTAException::bad_pointing_type
 *            Specified pointing is not a CTA pointing.
 *
 * Simulates a CTA event using the response function from an incident photon.
 * If the event is not detected a NULL pointer is returned.
 *
 * @todo Implement Phi dependence in CTA IRF
 * @todo Implement energy dispersion.
 ***************************************************************************/
GCTAEventAtom* GCTAResponse::mc(const double& area, const GPhoton& photon,
                                const GPointing& pnt, GRan& ran) const
{
    // Initialise event
    GCTAEventAtom* event = NULL;

    // Get pointer on CTA pointing
    const GCTAPointing *ctapnt = dynamic_cast<const GCTAPointing*>(&pnt);
    if (ctapnt == NULL)
        throw GCTAException::bad_pointing_type(G_MC);

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = ctapnt->zenith();
    double azimuth = ctapnt->azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta_p = ctapnt->dir().dist(photon.dir());
    double phi_p   = 0.0;  //TODO Implement Phi dependence

    // Compute effective area for photon
    double srcLogEng      = photon.energy().log10TeV();
    double effective_area = aeff(theta_p, phi_p, zenith, azimuth, srcLogEng);

    // Compute limiting value
    double ulimite = effective_area / area;

    // Continue only if event is detected
    if (ran.uniform() <= ulimite) {

        // Simulate offset from photon arrival direction
        double theta = psf_dummy_sigma(srcLogEng) * ran.chisq2() * rad2deg;
        double phi   = 360.0 * ran.uniform();

        // Rotate sky direction by offset
        GSkyDir sky_dir = photon.dir();
        sky_dir.rotate_deg(phi, theta);

        // Set measured photon arrival direction
        GCTAInstDir inst_dir;
        inst_dir.skydir(sky_dir);

        // Allocate event
        event = new GCTAEventAtom;

        // Set event attributes
        event->dir(inst_dir);
        event->energy(photon.energy());
        event->time(photon.time());

    } // endif: event was detected

    // Return event
    return event;
}


/***********************************************************************//**
 * @brief Set path to the calibration database
 *
 * @param[in] caldb Path to calibration database
 *
 * This method stores the CALDB root directory as the path to the CTA
 * calibration database. Once the final calibration format has been decided
 * on, this method should be adjusted.
 ***************************************************************************/
void GCTAResponse::caldb(const std::string& caldb)
{
    // Allocate calibration database
    GCaldb db(caldb);

    // Store the path to the calibration database
    m_caldb = db.dir();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load CTA response.
 *
 * @param[in] irfname Name of CTA response (without any file extension).
 *
 * The actually dummy version of the CTA response loads a CTA performance
 * table given in ASCII format into memory.
 ***************************************************************************/
void GCTAResponse::load(const std::string& irfname)
{
    // Save calibration database name
    std::string caldb = m_caldb;

    // Clear instance
    clear();

    // Restore calibration database name
    m_caldb = caldb;

    // Build filename
    std::string filename = m_caldb + "/" + irfname + ".dat";

    // Read performance table
    read_performance_table(filename);

    // Store response name
    m_rspname = irfname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA response information
 ***************************************************************************/
std::string GCTAResponse::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAResponse ===");
    result.append("\n"+parformat("Calibration database")+m_caldb);
    result.append("\n"+parformat("Response name")+m_rspname);
    if (m_offset_sigma == 0) {
        result.append("\n"+parformat("Offset angle dependence")+"none");
    }
    else {
        std::string txt = "Fixed sigma="+str(m_offset_sigma);
        result.append("\n"+parformat("Offset angle dependence")+txt);
    }
    /*
    result.append("\n"+parformat("Response definiton"));
    for (int i = 0; i < m_logE.size(); ++i) {
        result.append("\n"+parformat("logE="+str(m_logE.at(i))));
        result.append("Aeff="+str(m_aeff.at(i))+" m2");
        result.append(", r68="+str(m_r68.at(i))+" deg");
        result.append(", r80="+str(m_r68.at(i))+" deg");
    }
    */

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =              Model type dependent CTA response methods                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of extended source instrument response function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] model Extended source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_observation_type
 *            Specified observation is not a CTA observations.
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found.
 * @exception GCTAException::bad_instdir_type
 *            Instrument direction is not a valid CTA instrument direction.
 *
 * Performs integration of the model times IRF over the true photon arrival
 * direction in the coordinate system of the source model for azimuthally
 * independent models \f$M(\rho)\f$:
 * \f[\int_{\omega_{\rm min}}^{\omega_{\rm max}}
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} M(\rho) IRF(\rho, \omega)
 *    d\rho d\omega\f],
 *
 * The source centre is located at \f$\vec{m}\f$, and a spherical system
 * is defined around this location with \f$(\omega,\rho)\f$ being the
 * azimuth and zenith angles, respectively. \f$\omega=0\f$ is defined
 * by the direction that connects the source centre \f$\vec{m}\f$ to the
 * measured photon direction \f$\vec{p'}\f$, and \f$\omega\f$ increases
 * counterclockwise.
 ***************************************************************************/
double GCTAResponse::irf_extended(const GInstDir&             obsDir,
                                  const GEnergy&              obsEng,
                                  const GTime&                obsTime,
                                  const GModelExtendedSource& model,
                                  const GEnergy&              srcEng,
                                  const GTime&                srcTime,
                                  const GObservation&         obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL)
        throw GCTAException::bad_observation_type(G_IRF_EXTENDED);

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing(srcTime);
    if (pnt == NULL)
        throw GCTAException::no_pointing(G_IRF_EXTENDED);

    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&obsDir);
    if (dir == NULL)
        throw GCTAException::bad_instdir_type(G_IRF_EXTENDED);

    // Determine angular distance between measured photon direction and model
    // centre [radians]
    double zeta = model.dir().dist(dir->skydir());

    // Determine angular distance between measured photon direction and
    // pointing direction [radians]
    double eta = pnt->dir().dist(dir->skydir());

    // Determine angular distance between model centre and pointing direction
    // [radians]
    double lambda = model.dir().dist(pnt->dir());

    // Compute azimuth angle of pointing in model system [radians]
    // Will be comprised in interval [0,pi]
    double omega0 = 0.0;
    double denom  = std::sin(lambda) * std::sin(zeta);
    if (denom != 0.0) {
        double arg = (std::cos(eta) - std::cos(lambda) * std::cos(zeta))/denom;
        omega0     = arccos(arg);
    }

    // Get log10(E/TeV) of true and measured photon energies
    double srcLogEng = srcEng.log10TeV();
    double obsLogEng = obsEng.log10TeV();

    // Get PSF sigma
    double sigma = psf_dummy_sigma(srcLogEng);

    // Get maximum PSF and source radius in radians.
    double delta_max = psf_delta_max(0.0, 0.0, 0.0, 0.0, srcLogEng);
    double src_max   = model.radial()->theta_max();

    // Set radial model zenith angle range
    double rho_min = (zeta > delta_max) ? zeta - delta_max : 0.0;
    double rho_max = zeta + delta_max;
    if (rho_max > src_max)
        rho_max = src_max;

    // Initialise IRF value
    double irf = 0.0;

    // Perform zenith angle integration if interval is valid
    if (rho_max > rho_min) {

        // Setup integration kernel
        GCTAResponse::irf_kern_rho integrand(this,
                                             model.radial(),
                                             pnt->zenith(),
                                             pnt->azimuth(),
                                             srcLogEng,
                                             obsLogEng,
                                             sigma,
                                             zeta,
                                             lambda,
                                             omega0,
                                             delta_max);

        // Integrate over zenith angle
        GIntegral integral(&integrand);
        integral.eps(m_eps);
        irf = integral.romb(rho_min, rho_max);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(irf) || isinfinite(irf)) {
            std::cout << "*** ERROR: GCTAResponse::irf_extended:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", rho_min=" << rho_min;
            std::cout << ", rho_max=" << rho_max;
            std::cout << ", omega0=" << omega0 << ")";
            std::cout << std::endl;
        }
        #endif
    }

    // Compile option: Show integration results
    #if defined(G_DEBUG_IRF_EXTENDED)
    std::cout << "irf_extended";
    std::cout << " sigma=" << sigma;
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
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] model Diffuse source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs CTA observation.
 *
 * @exception GException::feature_not_implemented
 *            Diffuse source method is not yet implemented.
 *
 * @todo Implement method.
 ***************************************************************************/
double GCTAResponse::irf_diffuse(const GInstDir&            obsDir,
                                 const GEnergy&             obsEng,
                                 const GTime&               obsTime,
                                 const GModelDiffuseSource& model,
                                 const GEnergy&             srcEng,
                                 const GTime&               srcTime,
                                 const GObservation&        obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_DIFFUSE,
          "Diffuse IRF not yet implemented.");

    // Return IRF value
    return 0.0;
}


/***********************************************************************//**
 * @brief Return spatial integral of extended source model
 *
 * @param[in] model Extended source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 ***************************************************************************/
double GCTAResponse::npred_extended(const GModelExtendedSource& model,
                                    const GEnergy&              srcEng,
                                    const GTime&                srcTime,
                                    const GObservation&         obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL)
        throw GCTAException::bad_observation_type(G_NPRED_EXTENDED);

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(ctaobs->events());
    if (events == NULL)
        throw GException::no_list(G_NPRED_EXTENDED);

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius (radians)
    double sigma          = psf_dummy_sigma(srcLogEng);
    double psf_max_radius = psf_dummy_max(sigma);

    // Extract ROI radius (radians)
    double roi_radius = events->roi().radius() * deg2rad;

    // Compute distance between ROI and model centre (radians)
    double roi_model_distance = events->roi().centre().dist(model.radial()->dir());

    // Compute the ROI radius plus maximum PSF radius (radians). Any photon
    // coming from beyond this radius will not make it in the dataspace and
    // thus can be neglected.
    double roi_psf_radius = roi_radius + psf_max_radius;

    // Set offset angle integration range. We take here the ROI+PSF into
    // account to make no integrations beyond the point where the
    // contribution drops to zero.
    double theta_min = (roi_model_distance > roi_psf_radius)
                       ? roi_model_distance - roi_psf_radius: 0.0;
    double theta_max = model.radial()->theta_max();

    // Perform offset angle integration only if interval is valid
    if (theta_max > theta_min) {

        // Compute rotation matrix to convert from native model coordinates,
        // given by (theta,phi), into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        GMatrix rot;
        ry.eulery(model.radial()->dec() - 90.0);
        rz.eulerz(-model.radial()->ra());
        rot = transpose(ry * rz);
        
        // Compute position angle of ROI centre with respect to model
        // centre (radians)
        double phi = model.radial()->dir().posang(events->roi().centre().skydir());

        // Setup integration kernel
        GCTAResponse::npred_kern_theta integrand(this,
                                                 model.radial(),
                                                 &srcEng,
                                                 &srcTime,
                                                 ctaobs,
                                                 &rot,
                                                 roi_model_distance,
                                                 roi_psf_radius,
                                                 phi);

        // Integrate over theta
        GIntegral integral(&integrand);
        npred = integral.romb(theta_min, theta_max);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NPRED_EXTENDED)
        std::cout << "npred_extended:";
        std::cout << " theta_min=" << theta_min;
        std::cout << " theta_max=" << theta_max;
        std::cout << " npred=" << npred << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GCTAResponse::npred_extended:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", theta_min=" << theta_min;
        std::cout << ", theta_max=" << theta_max;
        std::cout << ")" << std::endl;
    }
    #endif
    
    // Return Npred
    return npred;
}


/*==========================================================================
 =                                                                         =
 =                    Low-level CTA response methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return effective area (in units of cm2)
 *
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * A simple offset angle dependence has been implemented using a Gaussian
 * if theta^2 (similar to the radial acceptance model implemented for the
 * background). The width of the Gaussian is assumed independent of energy
 * and given by the fixed parameter m_offset_sigma.
 *
 * @todo Implement offset angle dependence from response file.
 * @todo So far the parameters phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::aeff(const double& theta,
                          const double& phi,
                          const double& zenith,
                          const double& azimuth,
                          const double& srcLogEng) const
{
    // Interpolate effective area using node array
    double aeff = m_aeff_logE.interpolate(srcLogEng, m_aeff);
    if (aeff < 0)
        aeff = 0.0;
    
    // Optionally add in offset angle dependence
    if (m_offset_sigma != 0.0) {
        double offset = theta * rad2deg;
        double arg    = offset * offset / m_offset_sigma;
        double scale  = exp(-0.5 * arg * arg);
        aeff         *= scale;
    }

    // Return effective area
    return aeff;
}


/***********************************************************************//**
 * @brief Return point spread function (in units of sr^-1)
 *
 * @param[in] delta Angular separation between true and measured photon
 *            directions (radians).
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * @todo So far the parameters theta, phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::psf(const double& delta,
                         const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const double& srcLogEng) const
{
    // Determine energy dependent width of PSF
    double sigma = psf_dummy_sigma(srcLogEng);

    // Compute PSF
    double psf = psf_dummy(delta, sigma);

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
 * This method returns the maximum angular separation between true and
 * measured photon directions for which the PSF is non zero. The maximum
 * separation is actually fixed to 5 sigma, which corresponds to less than
 * 1e-5 of the central IRF value.
 *
 * @todo So far the parameters theta, phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::psf_delta_max(const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth,
                                   const double& srcLogEng) const
{
    // Determine energy dependent width of PSF
    double sigma = psf_dummy_sigma(srcLogEng);

    // Set maximum angular separation
    double delta_max = psf_dummy_max(sigma);

    // Return PSF
    return delta_max;
}


/***********************************************************************//**
 * @brief Return energy dispersion (in units or MeV^-1)
 *
 * @param[in] obsLogEng Log10 of measured photon energy (E/TeV).
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * @todo So far the parameters theta, phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::edisp(const double& obsLogEng,
                           const double& theta,
                           const double& phi,
                           const double& zenith,
                           const double& azimuth,
                           const double& srcLogEng) const
{
    // Dirac energy dispersion
    double edisp = (obsLogEng == srcLogEng) ? 1.0 : 0.0;

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
 *       (see comment abobr kluge below)
 ***************************************************************************/
double GCTAResponse::npsf(const GSkyDir&      srcDir,
                          const double&       srcLogEng,
                          const GTime&        srcTime,
                          const GCTAPointing& pnt,
                          const GCTARoi&      roi) const
{
    // Declare result
    double value = 0.0;

    // Extract relevant parameters from arguments
    double roi_radius       = roi.radius() * deg2rad;
    double roi_psf_distance = roi.centre().dist(srcDir);
    double sigma            = psf_dummy_sigma(srcLogEng);
    double rmax             = psf_dummy_max(sigma);

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
            GCTAResponse::npsf_kern_rad_azsym integrand(this,
                                                        roi_radius,
                                                        roi_psf_distance,
                                                        sigma);

            // Setup integration
            GIntegral integral(&integrand);
            integral.eps(m_eps);

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
            if (isnotanumber(value) || isinfinite(value)) {
                std::cout << "*** ERROR: GCTAResponse::npsf:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (value=" << value;
                std::cout << ", roi_radius=" << roi_radius;
                std::cout << ", roi_psf_distance=" << roi_psf_distance;
                std::cout << ", sigma=" << sigma;
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
 * @todo Implement integration over energy range.
 ***************************************************************************/
double GCTAResponse::nedisp(const GSkyDir&      srcDir,
                            const GEnergy&      srcEng,
                            const GTime&        srcTime,
                            const GCTAPointing& pnt,
                            const GEbounds&     ebds) const
{
    // Dummy
    double nedisp = 1.0;

    // Return integral
    return nedisp;
}


/*==========================================================================
 =                                                                         =
 =                   Analytical CTA PSF implementation                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return dummy point spread function (in units of sr^-1)
 *
 * @param[in] delta Angular separation between true and measured photon
 *            directions (radians).
 * @param[in] sigma Width of point spread function (radians).
 *
 * The Point Spread Function defines the probability density 
 * \f$d^2P/d\theta d\phi\f$
 * that a photon coming from direction "srcDir" is measured towards direction
 * "obsDir". The actual method implements a simple 2D Gaussian in small
 * angle approximation for the PSF.
 * The performance table quotes the size of the PSF as the 68%
 * containment radius \f$r_{68}\f$ in degrees. 
 * The containment radius \f$r\f$ is related to the 2D Gaussian 
 * \f$\sigma\f$ by the relation \f$r=\sigma \sqrt{-2 \ln (1-P)}\f$, where
 * \f$P\f$ is the containment fraction. For 68% one obtains
 * \f$\sigma=0.6624 \times r_{68}\f$.
 *
 * @todo The actual PSF is only valid in the small angle approximation.
 ***************************************************************************/
double GCTAResponse::psf_dummy(const double& delta, const double& sigma) const
{
    // Compute Psf value
    double sigma2 = sigma * sigma;
    double value  = exp(-0.5 * delta * delta / sigma2) / (twopi * sigma2);

    // Return PSF value
    return value;
}


/***********************************************************************//**
 * @brief Return width parameter of point spread function (in radians)
 *
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * This method returns the Gaussian sigma of the CTA PSF as function of
 * incident photon energy.
 ***************************************************************************/
double GCTAResponse::psf_dummy_sigma(const double& srcLogEng) const
{
    // Set conversion factor
    const double conv = 0.6624 * deg2rad;

    // Determine Gaussian sigma in radians
    double sigma = m_psf_logE.interpolate(srcLogEng, m_r68) * conv;

    // Return result
    return sigma;
}


/***********************************************************************//**
 * @brief Returns radius beyond which PSF is negligible (in radians)
 *
 * @param[in] sigma Width of point spread function (radians).
 ***************************************************************************/
double GCTAResponse::psf_dummy_max(const double& sigma) const
{
    // Return radius
    return (5.0*sigma);
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
 *                      Disk      Gauss      Shell
 * m_eps = 1e-3 : user 0m03.80s  0m13.41s   0m03.83s
 * m_eps = 1e-4 : user 0m03.85s  0m13.71s   0m04.68s
 * m_eps = 1e-5 : user 0m06.29s  0m23.22s   0m16.94s
 * m_eps = 1e-6 : user 0m12.65s  0m55.08s   1m32.52s
 ***************************************************************************/
void GCTAResponse::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();
    m_arffile.clear();
    m_rmffile.clear();
    m_psffile.clear();
    m_psf_logE.clear();
    m_aeff_logE.clear();
    m_logE.clear();
    m_aeff.clear();
    m_r68.clear();
    m_r80.clear();
    m_eps          = 1.0e-5; // Precision for Romberg integration
    m_offset_sigma = 3.0;    // Default value for now ...

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GCTAResponse::copy_members(const GCTAResponse& rsp)
{
    // Copy attributes
    m_caldb        = rsp.m_caldb;
    m_rspname      = rsp.m_rspname;
    m_arffile      = rsp.m_arffile;
    m_rmffile      = rsp.m_rmffile;
    m_psffile      = rsp.m_psffile;
    m_psf_logE     = rsp.m_psf_logE;
    m_aeff_logE    = rsp.m_aeff_logE;
    m_logE         = rsp.m_logE;
    m_aeff         = rsp.m_aeff;
    m_r68          = rsp.m_r68;
    m_r80          = rsp.m_r80;
    m_eps          = rsp.m_eps;
    m_offset_sigma = rsp.m_offset_sigma;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponse::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA performance table
 *
 * @param[in] filename Filename of CTA performance table.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method reads a CTA performance table given in the format that is
 * distributed within the CTA collaboration. Note that the effective area
 * is converted from m2 to cm2 and stored in units of cm2.
 ***************************************************************************/
void GCTAResponse::read_performance_table(const std::string& filename)
{
    // Clear arrays
    m_psf_logE.clear();
    m_aeff_logE.clear();
    m_logE.clear();
    m_aeff.clear();
    m_r68.clear();
    m_r80.clear();

    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Open performance table readonly
    FILE* fptr = std::fopen(filename.c_str(), "r");
    if (fptr == NULL)
        throw GCTAException::file_open_error(G_READ, filename);

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements
        std::vector<std::string> elements = split(line, " ");
        for (std::vector<std::string>::iterator it = elements.begin();
             it != elements.end(); ++it) {
            if (strip_whitespace(*it).length() == 0)
                elements.erase(it);
        }

        // Skip header
        if (elements[0].find("log(E)") != std::string::npos)
            continue;

        // Break loop if end of data table has been reached
        if (elements[0].find("----------") != std::string::npos)
            break;

        // Push elements in vectors
        m_logE.push_back(todouble(elements[0]));
        m_aeff.push_back(todouble(elements[1])*10000.0);
        m_r68.push_back(todouble(elements[2]));
        m_r80.push_back(todouble(elements[3]));

    } // endwhile: looped over lines

    // If we have nodes then setup node array
    int num = m_logE.size();
    if (num > 0) {
        for (int i = 0; i < num; ++i) {
            m_psf_logE.append(m_logE.at(i));
            m_aeff_logE.append(m_logE.at(i));
        }
    }

    // Close file
    std::fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA ARF matrix
 *
 * @param[in] hdu FITS table pointer.
 *
 * This method reads a CTA 1DC ARF matrix from the FITS HDU. Note that the
 * effective area is converted from m2 to cm2 and stored in units of cm2.
 *
 * @todo So far we expect energies in units of TeV
 * @todo So far we expect the effective area in units of m2
 ***************************************************************************/
void GCTAResponse::read_arf(const GFitsTable* hdu)
{
    // Clear arrays
    m_aeff_logE.clear();
    m_aeff.clear();

    // Get pointers to table columns
    const GFitsTableCol* energy_lo = &(*hdu)["ENERG_LO"];
    const GFitsTableCol* energy_hi = &(*hdu)["ENERG_HI"];
    const GFitsTableCol* specresp  = &(*hdu)["SPECRESP"];

    // Extract number of energy bins
    int num = energy_lo->number();

    // Set nodes
    for (int i = 0; i < num; ++i) {
    
        // Compute log10 mean energy
        double e_min = energy_lo->real(i);
        double e_max = energy_hi->real(i);
        double logE  = 0.5 * (log10(e_min) + log10(e_max));
        
        // Compute effective area in cm2
        double aeff = specresp->real(i) * 10000.0;
        
        // Store log10 mean energy and effective area value
        m_aeff_logE.append(logE);
        m_aeff.push_back(aeff);

    } // endfor: looped over nodes
    
    // Disable offset angle dependence
    m_offset_sigma = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Kernel for model zenith angle integration of IRF
 *
 * @param[in] rho Zenith angle with respect to model centre [radians].
 *
 * This method evaluates the kernel \f$K(\rho)\f$ for the zenith angle
 * integration
 * \f[\int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho) d\rho\f]
 * of the product between model and IRF, where
 * \f[K(\rho) = \int_{\omega_{\rm min}}^{\omega_{\rm max}} M(\rho)
 *              IRF(\rho, \omega) d\omega\f],
 * \f$M(\rho)\f$ is azimuthally symmetric the source model, and
 * \f$IRF(\rho, \omega)\f$ is the instrument response function.
 ***************************************************************************/
double GCTAResponse::irf_kern_rho::eval(double rho)
{
    // Compute half length of arc that lies within PSF validity circle
    // (in radians)
    double domega = 0.5 * cta_roi_arclength(rho,
                                            m_zeta,
                                            m_cos_zeta,
                                            m_sin_zeta,
                                            m_delta_max,
                                            m_cos_delta_max);

    // Initialise result
    double irf = 0.0;

    // Continue only if arc length is positive
    if (domega > 0.0) {

        // Compute omega integration range
        double omega_min = -domega;
        double omega_max = +domega;

        // Evaluate sky model M(rho)
        double model = m_radial->eval(rho);

        // Precompute cosine and sine terms for azimuthal integration
        double cos_rho = std::cos(rho);
        double sin_rho = std::sin(rho);
        double cos_psf = cos_rho*m_cos_zeta;
        double sin_psf = sin_rho*m_sin_zeta;
        double cos_ph  = cos_rho*m_cos_lambda;
        double sin_ph  = sin_rho*m_sin_lambda;

        // Setup integration kernel
        GCTAResponse::irf_kern_omega integrand(m_rsp,
                                               m_zenith,
                                               m_azimuth,
                                               m_srcLogEng,
                                               m_obsLogEng,
                                               m_sigma,
                                               m_zeta,
                                               m_lambda,
                                               m_omega0,
                                               rho,
                                               cos_psf,
                                               sin_psf,
                                               cos_ph,
                                               sin_ph);

        // Integrate over phi
        GIntegral integral(&integrand);
        integral.eps(m_rsp->m_eps);
        irf = integral.romb(omega_min, omega_max) * model * sin_rho;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(irf) || isinfinite(irf)) {
            std::cout << "*** ERROR: GCTAResponse::irf_kern_rho::eval";
            std::cout << "(rho=" << rho << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", domega=" << domega;
            std::cout << ", model=" << model;
            std::cout << ", sin_rho=" << sin_rho << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: arc length was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle IRF integration
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * This method evaluates the instrument response function
 * \f$IRF(\rho,\omega)\f$ for the azimuth angle integration of the IRF.
 *
 * From the model coordinates \f$(\rho,\omega)\f$ it computes the PSF
 * offset angle \f$\delta\f$, defined as the angle between true 
 * (\f$\vec{p}\f$) and observed (\f$\vec{p'}\f$) photon arrival direction,
 * using
 * \f[\delta = \arccos(\cos \rho \cos \zeta + 
 *                     \sin \rho \sin \zeta \cos \omega)\f]
 * where
 * \f$\zeta\f$ is the angular distance between the observed photon direction
 * \f$\vec{p}\f$ and the model centre \f$\vec{m}\f$.
 *
 * Furthermore, it computes the observed photon offset angle \f$\theta\f$,
 * defined as the angle between observed photon direction and camera pointing,
 * using
 * \f[\theta = \arccos(\cos \rho \cos \lambda + 
 *                     \sin \rho \sin \lambda \cos \omega_0 - \omega)\f]
 * where
 * \f$\lambda\f$ is the angular distance between the model centre and the
 * camera pointing direction.
 ***************************************************************************/
double GCTAResponse::irf_kern_omega::eval(double omega)
{
    // Compute PSF offset angle [radians]
    double delta = arccos(m_cos_psf + m_sin_psf * std::cos(omega));
    
    // Compute observed photon offset angle in camera system [radians]
    double theta = arccos(m_cos_ph + m_sin_ph * std::cos(m_omega0 - omega));
    
    //TODO: Compute true photon offset angle in camera system [radians]
    double phi = 0.0;

    // Evaluate IRF
    double irf = m_rsp->aeff(theta, phi, m_zenith, m_azimuth, m_srcLogEng) *
                 m_rsp->psf_dummy(delta, m_sigma);

    // Optionally take energy dispersion into account
    if (m_rsp->hasedisp() && irf > 0.0) {
        irf *= m_rsp->edisp(m_obsLogEng, theta, phi, m_zenith, m_azimuth, m_srcLogEng);
    }
    
    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(irf) || isinfinite(irf)) {
        std::cout << "*** ERROR: GCTAResponse::irf_kern_omega::eval";
        std::cout << "(omega=" << omega << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ", delta=" << delta;
        std::cout << ", theta=" << theta;
        std::cout << ", phi=" << phi << ")";
        std::cout << std::endl;
    }
    #endif

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for zenith angle Npred integration or radial model
 *
 * @param[in] theta Radial model zenith angle (radians).
 *
 * This method integrates a radial model for a given zenith angle theta over
 * all azimuth angles that fall within the ROI+PSF radius. The limitation to
 * an arc assures that the integration converges properly.
 *
 * @todo Rename method to npred_radial_kern_theta.
 ***************************************************************************/
double GCTAResponse::npred_kern_theta::eval(double theta)
{
    // Initialise Npred value
    double npred = 0.0;

    // Compute half length of arc that lies within ROI+PSF radius (radians)
    double dphi = 0.5 * cta_roi_arclength(theta,
                                          m_dist,
                                          m_cos_dist,
                                          m_sin_dist,
                                          m_radius,
                                          m_cos_radius);


    // Continue only if arc length is positive
    if (dphi > 0.0) {

        // Compute phi integration range
        double phi_min = m_phi - dphi;
        double phi_max = m_phi + dphi;

        // Get radial model value
        double model = m_radial->eval(theta);

        // Compute sine of offset angle
        double sin_theta = std::sin(theta);

        // Setup phi integration kernel
        GCTAResponse::npred_kern_phi integrand(m_rsp,
                                               m_srcEng,
                                               m_srcTime,
                                               m_obs,
                                               m_rot,
                                               theta,
                                               sin_theta);

        // Integrate over phi
        GIntegral integral(&integrand);
        npred = integral.romb(phi_min, phi_max) * sin_theta * model;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (isnotanumber(npred) || isinfinite(npred)) {
            std::cout << "*** ERROR: GCTAResponse::npred_kern_theta::eval";
            std::cout << "(theta=" << theta << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", model=" << model;
            std::cout << ", phi=[" << phi_min << "," << phi_max << "]";
            std::cout << ", sin_theta=" << sin_theta;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: arc length was positive

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Npred integration of radial model
 *
 * @param[in] phi Azimuth angle (radians).
 *
 * @todo Rename method to npred_radial_kern_phi.
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double GCTAResponse::npred_kern_phi::eval(double phi)
{
    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Compute Npred for this sky direction
    double npred = m_rsp->npred(srcDir, *m_srcEng, *m_srcTime, *m_obs);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GCTAResponse::npred_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", cos_phi=" << cos_phi;
        std::cout << ", sin_phi=" << sin_phi;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Integration kernel for npsf() method
 *
 * @param[in] theta Zenith angle with respect to PSF centre.
 *
 * This method implements the integration kernel needed for the npsf()
 * method.
 ***************************************************************************/
double GCTAResponse::npsf_kern_rad_azsym::eval(double theta)
{
    // Initialise PSF value
    double value = 0.0;
    
    // Get arclength for given radius in radians
    double phi = cta_roi_arclength(theta,
                                   m_psf,
                                   m_cospsf,
                                   m_sinpsf,
                                   m_roi,
                                   m_cosroi);

    // If arclength is positive then compute the PSF value
    if (phi > 0) {
    
        // Compute PSF value
        value = m_parent->psf_dummy(theta, m_sigma) * phi * std::sin(theta);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(value) || isinfinite(value)) {
            std::cout << "*** ERROR: GCTAResponse::npsf_kern_rad_azsym::eval";
            std::cout << "(theta=" << theta << ").";
            std::cout << " NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", theta=" << theta;
            std::cout << ", phi=" << phi << ")";
            std::cout << std::endl;
        }
        #endif
        
    } // endif: arclength was positive

    // Return
    return value;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
