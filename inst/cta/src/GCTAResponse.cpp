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
 * @brief GCTAResponse class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <vector>
#include <string>
#include <unistd.h>           // access() function
#include <stdio.h>            // fopen, fgets, fclose, etc...
#include "GCTAResponse.hpp"
#include "GCTAPointing.hpp"
#include "GCTAInstDir.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                           "GCTAResponse::caldb(std::string&)"
#define G_IRF_ATOM     "GCTAResponse::irf(GCTAEventAtom&, GModel&, GEnergy&,"\
                                                    " GTime&, GObservation&)"
#define G_IRF_BIN       "GCTAResponse::irf(GCTAEventBin&, GModel&, GEnergy&,"\
                                                    " GTime&, GObservation&)"
#define G_NPRED              "GCTAResponse::npred(GModel&, GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_READ           "GCTAResponse::read_performance_table(std::string&)"

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
 ***************************************************************************/
GCTAResponse::GCTAResponse(void) : GResponse()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief IRF constructor
 ***************************************************************************/
GCTAResponse::GCTAResponse(const std::string& rspname, const std::string& caldb)
                                                                   : GResponse()
{
    // Initialise class members for clean destruction
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
 * @param[in] rsp Response.
 ***************************************************************************/
GCTAResponse::GCTAResponse(const GCTAResponse& rsp) : GResponse(rsp)
{
    // Initialise class members for clean destruction
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
 * @param[in] rsp Response.
 ***************************************************************************/
GCTAResponse& GCTAResponse::operator= (const GCTAResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
 * @brief Return value of point source instrument response function.
 *
 * @param[in] event Event.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 ***************************************************************************/
double GCTAResponse::irf(const GEvent& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Get IRF value
    double rsp;
    if (event.isatom())
        rsp = irf(static_cast<const GCTAEventAtom&>(event), model,
                  srcEng, srcTime, obs);
    else
        rsp = irf(static_cast<const GCTAEventBin&>(event), model,
                  srcEng, srcTime, obs);

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return value of model IRF for event atom
 *
 * @param[in] event Event atom.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Diffuse IRF not yet implemented.
 ***************************************************************************/
double GCTAResponse::irf(const GCTAEventAtom& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Initialise response value
    double rsp = 0.0;

    // If model is a point source then return the point source IRF
    if (model.spatial()->isptsource()) {

        // Get point source location
        GSkyDir srcDir = static_cast<GModelSpatialPtsrc*>(model.spatial())->dir();

        // Compute IRF
        rsp = irf(event.dir(), event.energy(), event.time(),
                  srcDir, srcEng, srcTime, obs);

    } // endif: model was point source

    // ... otherwise return diffuse IRF
    else {

        // Feature not yet implemented
        throw GException::feature_not_implemented(G_IRF_ATOM,
              "Diffuse IRF not yet implemented.");

    } // endelse: model was not a point source

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return value of model IRF for event bin
 *        
 * @param[in] event Event bin.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Diffuse IRF not yet implemented.
 ***************************************************************************/
double GCTAResponse::irf(const GCTAEventBin& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Initialise response value
    double rsp = 0.0;

    // If model is a point source then return the point source IRF
    if (model.spatial()->isptsource()) {

        // Get point source location
        GSkyDir srcDir = static_cast<GModelSpatialPtsrc*>(model.spatial())->dir();

        // Compute IRF
        rsp = irf(event.dir(), event.energy(), event.time(),
                  srcDir, srcEng, srcTime, obs);

    } // endif: model was point source

    // ... otherwise return diffuse IRF
    else {

        // Feature not yet implemented
        throw GException::feature_not_implemented(G_IRF_BIN,
              "Diffuse IRF not yet implemented.");

    } // endelse: model was not a point source

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return value of point source instrument response function.
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observations.
 ***************************************************************************/
double GCTAResponse::irf(const GInstDir& obsDir, const GEnergy& obsEng,
                         const GTime& obsTime,
                         const GSkyDir& srcDir, const GEnergy& srcEng,
                         const GTime& srcTime, const GObservation& obs) const
{
    // Get pointing
    const GPointing *pnt = obs.pointing(srcTime);

    // Get point source IRF components
    double rsp  =  live(srcDir,  srcEng, srcTime, *pnt);
    rsp        *=  aeff(srcDir,  srcEng, srcTime, *pnt);
    rsp        *=   psf(obsDir,  srcDir, srcEng, srcTime, *pnt);
    rsp        *= edisp(obsEng,  srcDir, srcEng, srcTime, *pnt);
    rsp        *= tdisp(obsTime, srcDir, srcEng, srcTime, *pnt);

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return integral of instrument response function.
 *
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @param[in] roi Region of interest of data selection.
 * @param[in] ebds Energy boundaries of data selection.
 * @param[in] gti Good Time Intervals of data selection.
 ***************************************************************************/
double GCTAResponse::npred(const GModelSky& model, const GEnergy& srcEng,
                           const GTime& srcTime, const GObservation& obs) const
{
    // Initialise response value
    double rsp = 0.0;

    // If model is a point source then return the point source IRF
    if (model.spatial()->isptsource()) {

        // Get point source location
        GSkyDir srcDir = static_cast<GModelSpatialPtsrc*>(model.spatial())->dir();

        // Compute IRF
        rsp = npred(srcDir, srcEng, srcTime, obs);

    } // endif: model was point source

    // ... otherwise return diffuse IRF
    else {

        // Feature not yet implemented
        throw GException::feature_not_implemented(G_NPRED,
              "Diffuse IRF not yet implemented.");

    } // endelse: model was not a point source

    // Return response value
    return rsp;
}


/***********************************************************************//**
 * @brief Return integral of instrument response function.
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @param[in] roi Region of interest of data selection.
 * @param[in] ebds Energy boundaries of data selection.
 * @param[in] gti Good Time Intervals of data selection.
 ***************************************************************************/
double GCTAResponse::npred(const GSkyDir&  srcDir, const GEnergy& srcEng,
                           const GTime& srcTime, const GObservation& obs) const
{
    // Get pointers
    const GPointing *pnt  = obs.pointing(srcTime);
    GRoi            *roi  = ((GObservation*)&obs)->roi();
    GEbounds        *ebds = ((GObservation*)&obs)->ebounds();
    GGti            *gti  = ((GObservation*)&obs)->gti();

    // Get IRF components
    double nirf  =   live(srcDir, srcEng, srcTime, *pnt);
    nirf        *=   aeff(srcDir, srcEng, srcTime, *pnt);
    nirf        *=   npsf(srcDir, srcEng, srcTime, *pnt, *roi);
    nirf        *= nedisp(srcDir, srcEng, srcTime, *pnt, *ebds);
    nirf        *= ntdisp(srcDir, srcEng, srcTime, *pnt, *gti);

    // Return integrated IRF value
    return nirf;
}


/***********************************************************************//**
 * @brief Return livetime fraction (units: s s-1).
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * Dummy livetime fraction of 0.8.
 ***************************************************************************/
double GCTAResponse::live(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GPointing& pnt) const
{
    // Dummy
    double live = 0.8;

    // Return effective area
    return live;
}


/***********************************************************************//**
 * @brief Return effective area (units: cm2).
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * The actual implementation of this method assumes an effective area that
 * depends only on the true photon energy. No dependence on the photon's
 * arrival direction, observatory pointing and arrival time is assumed.
 ***************************************************************************/
double GCTAResponse::aeff(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GPointing& pnt) const
{
    // Get log(E)
    double logE = log10(srcEng.TeV());

    // Interpolate effective area using node array and convert to cm^2
    GNodeArray* nodes = (GNodeArray*)&m_nodes; // circumvent const correctness
    double aeff = nodes->interpolate(logE, m_aeff) * 10000.0;

    // Return effective area
    return aeff;
}


/***********************************************************************//**
 * @brief Return point spread function (units: sr^-1)
 *
 * @param[in] obsDir Pointer to observed photon direction.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] pnt Pointer to instrument pointing information (not used).
 *
 * The Point Spread Function defines the probability density 
 * \f$d^2P/d\theta d\phi\f$
 * that a photon coming from direction 'srcDir' is measured towards direction
 * 'obsDir'. The actual method implements a simple 2D Gaussian for the PSF.
 * The performance table quotes the size of the PSF as the 68%
 * containment radius \f$r_{68}\f$ in degrees. 
 * The containment radius \f$r\f$ is related to the 2D Gaussian 
 * \f$\sigma\f$ by the relation \f$r=\sigma \sqrt{-2 \ln (1-P)}\f$, where
 * \f$P\f$ is the containment fraction. For 68% one obtains
 * \f$\sigma=0.6624 \times r_{68}\f$.
 ***************************************************************************/
double GCTAResponse::psf(const GInstDir& obsDir,
                         const GSkyDir& srcDir, const GEnergy& srcEng,
                         const GTime& srcTime, const GPointing& pnt) const
{
    // Determine energy dependent width of PSF
    double sigma = psf_sigma(srcEng);

    // Determine angular separation between true and measured photon
    // direction in radians
    double theta = ((GCTAInstDir*)&obsDir)->dist(srcDir);

    // Get PSF value
    double value = psf(theta, sigma);

    // Return PSF value
    return value;
}


/***********************************************************************//**
 * @brief Return energy dispersion (units: MeV^-1)
 *
 * @param[in] obsEng Observed energy of photon.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * The actual implementation of this method assumes no energy dispersion,
 * which is equivalent of having a Dirac type energy dispersion.
 ***************************************************************************/
double GCTAResponse::edisp(const GEnergy& obsEng,
                           const GSkyDir& srcDir, const GEnergy& srcEng,
                           const GTime& srcTime, const GPointing& pnt) const
{
    // Dirac energy dispersion
    double edisp = (obsEng == srcEng) ? 1.0 : 0.0;

    // Return energy dispersion
    return edisp;
}


/***********************************************************************//**
 * @brief Return time dispersion (units: s^-1)
 *
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * The actual implementation of this method assumes no time dispersion,
 * which is equivalent of having a Dirac type time dispersion.
 ***************************************************************************/
double GCTAResponse::tdisp(const GTime& obsTime,
                           const GSkyDir& srcDir, const GEnergy& srcEng,
                           const GTime& srcTime, const GPointing& pnt) const
{
    // Dirac time dispersion
    double tdisp = (obsTime == srcTime) ? 1.0 : 0.0;

    // Return time dispersion
    return tdisp;
}


/***********************************************************************//**
 * @brief Return result of PSF integration over ROI.
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] pnt Pointer to instrument pointing information (nut used).
 * @param[in] roi Region of interest of data selection.
 ***************************************************************************/
double GCTAResponse::npsf(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GPointing& pnt,
                          const GRoi& roi) const
{
    // Get pointer to CTA ROI (bypass const correctness)
    GCTARoi* ctaroi = (GCTARoi*)&roi;

    // Extract relevant parameters from arguments
    double radroi = ctaroi->radius() * deg2rad;
    double psf    = ctaroi->centre().dist(srcDir);
    double sigma  = psf_sigma(srcEng);

    // Compute integrated PSF
    double value = npsf(psf, radroi, sigma);

    // Return integral
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
 * @brief Return result of time dispersion integral over GTIs
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 * @param[in] gti Good Time Intervals of data selection.
 *
 * @todo Implement integration over GTIs.
 ***************************************************************************/
double GCTAResponse::ntdisp(const GSkyDir& srcDir, const GEnergy& srcEng,
                            const GTime& srcTime, const GPointing& pnt,
                            const GGti& gti) const
{
    // Dummy
    double ntdisp = 1.0;

    // Return integral
    return ntdisp;
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
 * @todo Implement energy dispersion.
 ***************************************************************************/
GCTAEventAtom* GCTAResponse::mc(const double& area, const GPhoton& photon,
                                const GPointing& pnt, GRan& ran) const
{
    // Initialise event
    GCTAEventAtom* event = NULL;

    // Compute effective area for photon
    double effective_area = aeff(photon.dir(), photon.energy(), photon.time(), pnt);

    // Compute livetime fraction for photon
    double livetime_fraction = live(photon.dir(), photon.energy(), photon.time(), pnt);

    // Compute limiting value
    double ulimite = (effective_area*livetime_fraction) / area;

    // Continue only if event is detected
    if (ran.uniform() <= ulimite) {

        // Simulate offset from photon arrival direction
        double theta = psf_sigma(photon.energy()) * ran.chisq2() * rad2deg;
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
 * @brief Print CTA response information
 ***************************************************************************/
std::string GCTAResponse::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAResponse ===\n");
    result.append(parformat("Calibration database")+m_caldb+"\n");
    result.append(parformat("Response name")+m_rspname+"\n");
    result.append(parformat("Response definiton"));
    for (int i = 0; i < m_logE.size(); ++i) {
        result.append("\n"+parformat("logE="+str(m_logE.at(i))));
        result.append("Aeff="+str(m_aeff.at(i))+" m2");
        result.append(", r68="+str(m_r68.at(i))+" deg");
        result.append(", r80="+str(m_r68.at(i))+" deg");
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                         CTA specific IRF methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return point spread function (in units of sr^-1)
 *
 * @param[in] theta Angular separation between true and measured photon
 *            directions (radians).
 * @param[in] sigma Width of point spread function (radians).
 *
 * A simple Gaussian function is assumed to describe the CTA point spread
 * function.
 ***************************************************************************/
double GCTAResponse::psf(const double& theta, const double& sigma) const
{
    // Compute Psf value
    double sigma2 = sigma * sigma;
    double value  = exp(-0.5 * theta * theta / sigma2) / (twopi * sigma2);

    // Return PSF value
    return value;
}


/***********************************************************************//**
 * @brief Return width parameter of point spread function (in radians)
 *
 * @param[in] srcEng True energy of photon.
 *
 * This method returns the Gaussian sigma of the CTA PSF as function of
 * incident photon energy.
 ***************************************************************************/
double GCTAResponse::psf_sigma(const GEnergy& srcEng) const
{
    // Get log(E)
    double logE = log10(srcEng.TeV());

    // Determine Gaussian sigma in radians
    GNodeArray* nodes = (GNodeArray*)&m_nodes; // bypass const correctness
    double sigma  = nodes->interpolate(logE, m_r68) * 0.6624 * deg2rad;

    // Return result
    return sigma;
}


/***********************************************************************//**
 * @brief Integrate the PSF over the ROI.
 *
 * @param[in] psf Angular separation between PSF-ROI centres (radians).
 * @param[in] radroi Radius of ROI (radians).
 * @param[in] sigma PSF width parameter (radians).
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
double GCTAResponse::npsf(const double& psf, const double& radroi,
                          const double& sigma) const
{
    // Declare result
    double value;

    // Get maximum PSF radius
    double rmax = 5.0*sigma;

    // If PSF is sufficiently enclosed by ROI, skip the numerical integration
    // and assume that the integral is 1.0
    if (psf+rmax < radroi)
        value = 1.0;

    // ... otherwise perform numerical integration
    else {

        // Pre-computations
        double cosroi = cos(radroi);
        double cospsf = cos(psf);
        double sinpsf = sin(psf);

        // Setup integration function
        GCTAResponse::npsf_kern_rad_azsym integrand(this, radroi, cosroi,
                                                    psf, cospsf, sinpsf, sigma);
        GIntegral integral(&integrand);

        // Integrate PSF
        value = integral.romb(0.0, rmax);

    } // endelse: numerical integration required

    // Return integral
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
    double phi = m_parent->npsf_kern_azsym(theta, m_roi, m_cosroi, m_psf, 
                                           m_cospsf, m_sinpsf);

    // Get PSF value
    double value = m_parent->psf(theta, m_sigma) * phi * theta;

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Integration kernel for npsf integration for azimuthally symmetric PSF
 *
 * @param[in] rad Radial distance of PSF centre in radians (<pi).
 * @param[in] roi ROI radius in radians.
 * @param[in] cosroi Cosine of ROI radius.
 * @param[in] psf PSF distance to ROI centre in radians (<pi).
 * @param[in] cospsf Cosine of PSF distance to ROI centre.
 * @param[in] sinpsf Sinus of PSF distance to ROI centre.
 *
 * This method returns the arclength in radians of a circle of radius 'rad'
 * with a centre that is offset by 'psf' from the ROI centre, where the ROI
 * radius is given by 'roi'. To speed-up computations, the cosines and sinus
 * of 'roi' and 'psf' should be calculated by the client and be passed to
 * the method.
 ***************************************************************************/
double GCTAResponse::npsf_kern_azsym(const double& rad,
                                     const double& roi, const double& cosroi,
                                     const double& psf, const double& cospsf,
                                     const double& sinpsf) const
{
    // Declare arclength
    double arclength;

    // Handle special case of identical PSF and ROI centres
    if (psf == 0.0) {
        if (rad > roi) arclength = 0.0;   // PSF radius outside ROI
        else           arclength = twopi; // PSF radius inside ROI
    }

    // ... PSF and ROI centres are not identical
    else {

        // Handle special case of zero radial distance to PSF centre
        if (rad == 0.0) {
            if (psf > roi) arclength = 0.0;   // PSF centre outside ROI
            else           arclength = twopi; // PSF centre inside ROI
        }

        // Handle general case
        else {
            double dist = roi - psf;
            if (-rad >= dist) 
                arclength = 0.0;
            else if (rad <= dist) 
                arclength = twopi;
            else {
                double cosrad = cos(rad);
                double sinrad = sin(rad);
                double cosang = (cosroi - cospsf*cosrad) / (sinpsf*sinrad);
                arclength = acos(cosang);
            }
        }

    } // endelse: PSF and ROI centres were not identical

    // Return arclength
    return arclength;
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
    FILE* fptr = fopen(filename.c_str(), "r");
    if (fptr == NULL)
        throw GCTAException::file_open_error(G_READ, filename);

    // Read lines
    while (fgets(line, n, fptr) != NULL) {

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
    fclose(fptr);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
