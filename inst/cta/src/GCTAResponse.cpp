/***************************************************************************
 *                  GCTAResponse.cpp  -  CTA Response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAResponse.cpp
 * @brief CTA response class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>           // access() function
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include <vector>
#include <string>
#include "GModelSpatialPtsrc.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"
#include "GVector.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GCTAPointing.hpp"
#include "GCTAEventList.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                           "GCTAResponse::caldb(std::string&)"
#define G_IRF      "GCTAResponse::irf(GEvent&, GModelSky&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_NPRED1             "GCTAResponse::npred(GModel&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_NPRED2            "GCTAResponse::npred(GSkyDir&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE        "GCTAResponse::irf_diffuse(GInstDir&, GEnergy&,"\
                      " GTime&, GModelSky&, GEnergy&, GTime&, GObservation&)"
#define G_READ           "GCTAResponse::read_performance_table(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_IRF_DIFFUSE 1

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
 * @brief Response constructor
 *
 * @param[in] rspname Response file name.
 * @param[in] caldb Calibration database.
 *
 * Create instance by specifying the response file name and the calibration
 * database path.
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
GCTAResponse& GCTAResponse::operator= (const GCTAResponse& rsp)
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
 * @param[in] event Event.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GCTAException::bad_observation_type
 *            Observation is not a CTA observation.
 * @exception GCTAException::bad_instdir_type
 *            Instrument direction is not a CTA instrument direction.
 ***************************************************************************/
double GCTAResponse::irf(const GEvent& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL)
        throw GCTAException::bad_observation_type(G_IRF);

    // Get pointer on CTA instrument direction
    const GCTAInstDir* ctadir = dynamic_cast<const GCTAInstDir*>(&event.dir());
    if (ctadir == NULL)
        throw GCTAException::bad_instdir_type(G_IRF);

    // Get pointer on point source spatial component
    GModelSpatialPtsrc* ptsrc = dynamic_cast<GModelSpatialPtsrc*>(model.spatial());

    // If model is a point source then compute point source IRF
    if (ptsrc != NULL)
        irf = irf_ptsrc(*ctadir, event.energy(), event.time(),
                        ptsrc->dir(), srcEng, srcTime, *ctaobs);

    // ... otherwise compute diffuse IRF
    else
        irf = irf_diffuse(*ctadir, event.energy(), event.time(),
                          model, srcEng, srcTime, *ctaobs);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return spatial integral of instrument response function
 *
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @todo Implement diffuse computation.
 ***************************************************************************/
double GCTAResponse::npred(const GModelSky& model, const GEnergy& srcEng,
                           const GTime& srcTime, const GObservation& obs) const
{
    // Initialise response value
    double rsp = 0.0;

    // Get model pointers
    GModelSpatialPtsrc* ptsrc =
                        dynamic_cast<GModelSpatialPtsrc*>(model.spatial());

    // If model is a point source then compute the point source IRF
    if (ptsrc != NULL)
        rsp = npred(ptsrc->dir(), srcEng, srcTime, obs);

    // ... otherwise return diffuse IRF
    else {

        // Feature not yet implemented
        throw GException::feature_not_implemented(G_NPRED1,
              "Diffuse IRF not yet implemented.");

    } // endelse: model was not a point source

    // Return response value
    return rsp;
}


/***********************************************************************//**
 * @brief Simulate event from photon
 *
 * @param[in] area Simulation surface area.
 * @param[in] photon Photon.
 * @param[in] pnt Pointing.
 * @param[in] ran Random number generator.
 *
 * Simulates a CTA event using the response function from an incident photon.
 * If the event is not detected a NULL pointer is returned.
 *
 * @todo Implement computation of telescope pointing and radial offset and
 *       polar angles.
 * @todo Implement energy dispersion.
 ***************************************************************************/
GCTAEventAtom* GCTAResponse::mc(const double& area, const GPhoton& photon,
                                const GPointing& pnt, GRan& ran) const
{
    // Initialise event
    GCTAEventAtom* event = NULL;

    // Compute effective area for photon
    double theta          = 0.0;
    double phi            = 0.0;
    double zenith         = 0.0;
    double azimuth        = 0.0;
    double srcLogEng      = photon.energy().log10TeV();
    double effective_area = aeff(theta, phi, zenith, azimuth, srcLogEng);

    // Compute limiting value
    double ulimite = effective_area / area;

    // Continue only if event is detected
    if (ran.uniform() <= ulimite) {

        // Simulate offset from photon arrival direction
        double theta = psf_dummy_sigma(srcLogEng) * ran.chisq2() * rad2deg;
        double phi   = 360.0 * ran.uniform();

        // Rotate sky direction by offset
        GSkyDir sky_dir = photon.dir();
        sky_dir.rotate(phi, theta);

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
 * @exception GException::caldb_not_found
 *            Calibration database repository not found.
 *
 * This default method simply checks if the calibration database directory
 * exists. If the directory exists, the path will be stored. No checking is
 * implemented that checks for the consistency of the calibration database.
 *
 * @todo Implement a GCalDB class that handles any calibration database
 *       issues. GCalDB may be an abstract class for which instrument
 *       specific methods are implement to handle any instrument specific
 *       IRF database issues. 
 ***************************************************************************/
void GCTAResponse::caldb(const std::string& caldb)
{
    // Check if calibration database directory is accessible
    if (access(caldb.c_str(), R_OK) != 0)
        throw GException::caldb_not_found(G_CALDB, caldb);

    // Store the path to the calibration database
    m_caldb = caldb;

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
 =                         CTA specific IRF methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of point source instrument response function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observations.
 *
 * @todo The telescope zenith and azimuth angles as well as the radial offset
 *       and polar angles in the camera are not yet computed correctly (as
 *       the actual IRF implement does not need these values).
 ***************************************************************************/
double GCTAResponse::irf_ptsrc(const GCTAInstDir& obsDir, const GEnergy& obsEng,
                               const GTime& obsTime,
                               const GSkyDir& srcDir, const GEnergy& srcEng,
                               const GTime& srcTime, const GCTAObservation& obs) const
{
    // Get CTA pointing
    const GCTAPointing* pnt = obs.pointing(srcTime);

    // Get pointing direction zenith angle and azimuth
    double zenith  = 0.0;
    double azimuth = 0.0;

    // Get radial offset and polar angles in camera
    double theta = 0.0;
    double phi   = 0.0;

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

        // Get point source IRF
        irf  = aeff(theta, phi, zenith, azimuth, srcLogEng);
        irf *= psf(delta, theta, phi, zenith, azimuth, srcLogEng);

        // Multiply-in energy dispersion
        if (hasedisp()) {

            // Get log10(E/TeV) of measured photon energy.
            double obsLogEng = obsEng.log10TeV();

            // Multiply-in energy dispersion
            irf *= edisp(obsLogEng, theta, phi, zenith, azimuth, srcLogEng);

        } // endif: energy dispersion was available

    } // endif: we were sufficiently close to PSF

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of diffuse source instrument response function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] model Diffuse model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observations.
 *
 * @todo The telescope zenith and azimuth angles as well as the radial offset
 *       and polar angles in the camera are not yet computed correctly (as
 *       the actual IRF implement does not need these values).
 ***************************************************************************/
double GCTAResponse::irf_diffuse(const GCTAInstDir& obsDir, const GEnergy& obsEng,
                                 const GTime& obsTime,
                                 const GModelSky& model, const GEnergy& srcEng,
                                 const GTime& srcTime, const GCTAObservation& obs) const
{
    // Get CTA pointing
    const GCTAPointing* pnt = obs.pointing(srcTime);

    // Get pointing direction zenith angle and azimuth
    double zenith  = 0.0;
    double azimuth = 0.0;

    // Get radial offset and polar angles in camera
    double theta = 0.0;
    double phi   = 0.0;

    // Get log10(E/TeV) of true and measured photon energies
    double srcLogEng = srcEng.log10TeV();
    double obsLogEng = obsEng.log10TeV();

    // Get PSF sigma
    double sigma = psf_dummy_sigma(srcLogEng);

    // Get maximum angular separation for which PSF is significant
    double delta_max = psf_delta_max(theta, phi, zenith, azimuth, srcLogEng);

    // Determine angular distance and position angle of measured photon
    // direction with respect to pointing direction (in radians)
    double dist = pnt->dir().dist(obsDir.skydir());
    double pa   = pnt->dir().posang(obsDir.skydir());

    // Setup integration kernel
    GCTAResponse::psf_kern_theta integrand(this, pnt, model.spatial(),
                                           dist, pa, delta_max,
                                           zenith, azimuth,
                                           srcLogEng, obsLogEng, sigma);

    // Integrate over theta
    GIntegral integral(&integrand);
    double    theta_min = (dist > delta_max) ? dist - delta_max : 0.0;
    double    theta_max = dist + delta_max;
    double    irf       = integral.romb(theta_min, theta_max);

    // Compile option: Show integration results
    #if G_DEBUG_IRF_DIFFUSE
    std::cout << "irf_diffuse";
    std::cout << " sigma=" << sigma;
    std::cout << " theta_min=" << theta_min;
    std::cout << " theta_max=" << theta_max;
    std::cout << " irf=" << irf << std::endl;
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return effective area (units: cm2)
 *
 * @param[in] theta Radial offset angle in camera (radians).
 * @param[in] phi Polar angle in camera (radians).
 * @param[in] zenith Zenith angle of telescope pointing (radians).
 * @param[in] azimuth Azimuth angle of telescope pointing (radians).
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 *
 * @todo So far the parameters theta, phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::aeff(const double& theta, const double& phi,
                          const double& zenith, const double& azimuth,
                          const double& srcLogEng) const
{
    // Interpolate effective area using node array and convert to cm^2
    double aeff = m_nodes.interpolate(srcLogEng, m_aeff) * 10000.0;

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
                         const double& theta, const double& phi,
                         const double& zenith, const double& azimuth,
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
double GCTAResponse::psf_delta_max(const double& theta, const double& phi,
                                   const double& zenith, const double& azimuth,
                                   const double& srcLogEng) const
{
    // Determine energy dependent width of PSF
    double sigma = psf_dummy_sigma(srcLogEng);

    // Set maximum angular separation
    double delta_max = sigma * 5.0;

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
                           const double& theta, const double& phi,
                           const double& zenith, const double& azimuth,
                           const double& srcLogEng) const
{
    // Dirac energy dispersion
    double edisp = (obsLogEng == srcLogEng) ? 1.0 : 0.0;

    // Return energy dispersion
    return edisp;
}


/***********************************************************************//**
 * @brief Return spatial integral of instrument response function
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GException::no_list
 *            Observation does not contain an event list.
 ***************************************************************************/
double GCTAResponse::npred(const GSkyDir& srcDir, const GEnergy& srcEng,
                           const GTime& srcTime, const GObservation& obs) const
{
    // Get pointer on CTA events list
    const GCTAEventList* ptr = dynamic_cast<const GCTAEventList*>(obs.events());
    if (ptr == NULL)
        throw GException::no_list(G_NPRED2);

    // Get pointing direction zenith angle and azimuth
    const GPointing *pnt = obs.pointing(srcTime);
    double zenith  = 0.0;
    double azimuth = 0.0;

    // Get radial offset and polar angles in camera
    double theta = 0.0;
    double phi   = 0.0;

    // Get log10(E/TeV) of true photon energy.
    double srcLogEng = srcEng.log10TeV();

    // Get IRF components
    double nirf = aeff(theta, phi, zenith, azimuth, srcLogEng);
    nirf       *= npsf(srcDir, srcLogEng, srcTime, *pnt, ptr->roi());

    // Multiply-in energy dispersion
    if (hasedisp()) {

        // Multiply-in energy dispersion
        nirf *= nedisp(srcDir, srcEng, srcTime, *pnt, ptr->ebounds());

    } // endif: had energy dispersion

    // Return integrated IRF value
    return nirf;
}


/***********************************************************************//**
 * @brief Return result of PSF integration over ROI.
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] pnt Pointer to instrument pointing information (nut used).
 * @param[in] roi Region of interest of data selection.
 *
 * This method integrates the PSF over the circular region of interest.
 * Integration is done in a polar coordinate system centred on the PSF since
 * the PSF is assumed to be azimuthally symmetric. The polar integration is
 * done using the method npsf_kern_rad_azsym() that computes analytically
 * the arclength that is comprised within the ROI. The radial integration
 * is done using a standard Romberg integration. Note that the integration
 * is only performed when the PSF is sufficiently close to the ROI border so
 * that part of the PSF may be outside the ROI. For all other cases, the
 * integral is simply 1.
 ***************************************************************************/
double GCTAResponse::npsf(const GSkyDir& srcDir, const double& srcLogEng,
                          const GTime& srcTime, const GPointing& pnt,
                          const GRoi& roi) const
{
    // Declare result
    double value = 0.0;

    // Get pointer to CTA ROI
    GCTARoi* ctaroi = (GCTARoi*)&roi;

    // Extract relevant parameters from arguments
    double radroi = ctaroi->radius() * deg2rad;
    double psf    = ctaroi->centre().dist(srcDir);
    double sigma  = psf_dummy_sigma(srcLogEng);

    // Get maximum PSF radius
    double rmax = 5.0*sigma;

    // If PSF is sufficiently enclosed by ROI, skip the numerical integration
    // and assume that the integral is 1.0
    if (psf+rmax < radroi)
        value = 1.0;

    // ... otherwise perform numerical integration
    else {

        // Setup integration function
        GCTAResponse::npsf_kern_rad_azsym integrand(this, radroi, psf, sigma);
        GIntegral                         integral(&integrand);

        // Integrate PSF
        value = integral.romb(0.0, rmax);

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
 * @param[in] pnt Pointer to instrument pointing information.
 * @param[in] ebds Energy boundaries of data selection.
 *
 * @todo Implement integration over energy range.
 ***************************************************************************/
double GCTAResponse::nedisp(const GSkyDir& srcDir, const GEnergy& srcEng,
                            const GTime& srcTime, const GPointing& pnt,
                            const GEbounds& ebds) const
{
    // Dummy
    double nedisp = 1.0;

    // Return integral
    return nedisp;
}


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
    double sigma = m_nodes.interpolate(srcLogEng, m_r68) * conv;

    // Return result
    return sigma;
}


/***********************************************************************//**
 * @brief Integration kernel for IRF radial offset angle integration
 *
 * @param[in] theta Radial offset angle (radians).
 *
 * This method provides the integration kernel for the IRF radial offset
 * angle integration. The radial offset angle is given in the camera system.
 ***************************************************************************/
double GCTAResponse::psf_kern_theta::eval(double theta)
{
    // Compute half interval length
    double delta_phi = 0.5 * cta_roi_arclength(theta, 
                                               m_dist, m_cos_dist, m_sin_dist,
                                               m_delta_max, m_cos_delta_max);

    // Set phi interval
    double phi_min = m_pa - delta_phi;
    double phi_max = m_pa + delta_phi;

    // Compute cosine and sine terms for PSF delta computation
    double cos_term = std::cos(theta) * m_cos_dist;
    double sin_term = std::sin(theta) * m_sin_dist;

    // Setup integration kernel
    GCTAResponse::psf_kern_phi integrand(m_rsp, m_pnt, m_spatial,
                                         theta, m_pa,
                                         m_zenith, m_azimuth,
                                         m_srcLogEng, m_obsLogEng, m_sigma,
                                         cos_term, sin_term);

    // Integrate over phi
    GIntegral integral(&integrand);
    double    value = integral.romb(phi_min, phi_max) * std::sin(theta);

    // Return result
    return value;
}


/***********************************************************************//**
 * @brief Integration kernel for IRF polar angle integration
 *
 * @param[in] phi Polar angle (radians).
 *
 * This method provides the integration kernel for the IRF polar angle
 * integration. The polar angle is given in the camera system.
 ***************************************************************************/
double GCTAResponse::psf_kern_phi::eval(double phi)
{
    // Initialise result
    double value = 0.0;

    // Transform theta and phi in sky coordinates
    double  cos_phi   = std::cos(-phi);
    double  sin_phi   = std::sin(-phi);
    double  cos_theta = std::cos(-m_theta);
    double  sin_theta = std::sin(-m_theta);
    GVector vector(cos_phi*sin_theta, sin_phi*sin_theta, cos_theta);
    GVector vector_rot = m_pnt->rot() * vector;
    GSkyDir skydir;
    skydir.celvector(vector_rot);

    // Evaluate sky model
    double model = m_spatial->eval(skydir);

    // Continue only if model is valid
    if (model > 0.0) {

        // Compute PSF offset angle
        double delta = arccos(m_cos_term + m_sin_term * std::cos(phi - m_pa));

        // Evaluate IRF
        double irf = m_rsp->aeff(m_theta, phi, m_zenith, m_azimuth, m_srcLogEng) *
                     m_rsp->psf_dummy(delta, m_sigma) *
                     m_rsp->edisp(m_obsLogEng, m_theta, phi, m_zenith, m_azimuth, m_srcLogEng);

        // Compute value
        value = model * irf;

    } // endif: model was valid

    // Return
    return value;
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
    // Get arclength for given radius in radians
    double phi = cta_roi_arclength(theta, m_psf, m_cospsf, m_sinpsf, m_roi, m_cosroi);

    // Get PSF value
    double value = m_parent->psf_dummy(theta, m_sigma) * phi * sin(theta);

    // Return
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
void GCTAResponse::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();
    m_logE.clear();
    m_aeff.clear();
    m_r68.clear();
    m_r80.clear();

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
    m_caldb   = rsp.m_caldb;
    m_rspname = rsp.m_rspname;
    m_nodes   = rsp.m_nodes;
    m_logE    = rsp.m_logE;
    m_aeff    = rsp.m_aeff;
    m_r68     = rsp.m_r68;
    m_r80     = rsp.m_r80;

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
 * distributed within the CTA collaboration.
 ***************************************************************************/
void GCTAResponse::read_performance_table(const std::string& filename)
{
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
        m_aeff.push_back(todouble(elements[1]));
        m_r68.push_back(todouble(elements[2]));
        m_r80.push_back(todouble(elements[3]));

    } // endwhile: looped over lines

    // If we have nodes then setup node array
    int num = m_logE.size();
    if (num > 0) {
        for (int i = 0; i < num; ++i)
            m_nodes.append(m_logE.at(i));
    }

    // Close file
    std::fclose(fptr);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
