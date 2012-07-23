/***************************************************************************
 *                  GCTAResponse.cpp  -  CTA Response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
#include <cmath>
#include <vector>
#include <string>
#include "GModelSpatialPtsrc.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GCaldb.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GCTAResponse_helpers.hpp"
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
#define G_NPRED_EXTENDED                      "GCTAResponse::npred_extended("\
                       "GModelExtendedSource&,GEnergy&,GTime&,GObservation&)"
#define G_NPRED_DIFFUSE                        "GCTAResponse::npred_diffuse("\
                        "GModelDiffuseSource&,GEnergy&,GTime&,GObservation&)"
#define G_READ           "GCTAResponse::read_performance_table(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_READ_ARF                         //!< Debug read_arf method
//#define G_DEBUG_IRF_EXTENDED                 //!< Debug irf_extended method
//#define G_DEBUG_IRF_DIFFUSE                  //!< Debug irf_diffuse method
//#define G_DEBUG_NPRED_EXTENDED             //!< Debug npred_extended method
//#define G_DEBUG_NPRED_DIFFUSE               //!< Debug npred_diffuse method
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
 * Creates an empty instance of a CTA response object.
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
 *
 * Copies an instance of a CTA response object. Note that a deep copy is
 * performed, hence the original object can be destroyed without any loss
 * of information in the copy.
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
 * Create instance of a CTA response object by specifying the response file
 * name and the calibration database path.
 *
 * If an empty string is passed as calibration database path, the method will
 * use the CALDB environment variable to determine calibration database path.
 * This is done in the method GCTAResponse::caldb which makes use of a
 * GCaldb object to locate the calibration database.
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
 *
 * Destroys instance of CTA response object.
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
 *
 * Assigns CTA response object to another CTA response object. The assignment
 * performs a deep copy of all information, hence the original object from
 * which the assignment has been performed can be destroyed after this
 * operation without any loss of information.
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
 *
 * Clears CTA response object by resetting all members to an initial state.
 * Any information that was present in the object before will be lost.
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
 *
 * Creates a clone (deep copy) of a CTA response object.
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
    if (ctaobs == NULL) {
        throw GCTAException::bad_observation_type(G_IRF);
    }

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing();
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_IRF);
    }

    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&obsDir);
    if (dir == NULL) {
        throw GCTAException::bad_instdir_type(G_IRF);
    }

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
 * @todo Implement Phi dependence in CTA IRF
 * @todo Write method documentation
 ***************************************************************************/
double GCTAResponse::npred(const GSkyDir&      srcDir,
                           const GEnergy&      srcEng,
                           const GTime&        srcTime,
                           const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        throw GCTAException::bad_observation_type(G_NPRED);
    }

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing();
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_NPRED);
    }

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(ctaobs->events());
    if (events == NULL) {
        throw GException::no_list(G_NPRED);
    }

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
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 *
 * @exception GCTAException::no_pointing
 *            No CTA pointing found in observation.
 *
 * Simulates a CTA event using the response function from an incident photon.
 * If the event is not detected a NULL pointer is returned.
 *
 * The method also applies a deadtime correction using a Monte Carlo process,
 * taking into account temporal deadtime variations. For this purpose, the
 * method makes use of the time dependent GObservation::deadc method.
 *
 * @todo Implement Phi dependence in CTA IRF
 * @todo Implement energy dispersion
 * @todo Implement the method for a 3 Gaussian PSF
 ***************************************************************************/
GCTAEventAtom* GCTAResponse::mc(const double& area, const GPhoton& photon,
                                const GObservation& obs, GRan& ran) const
{
    // Initialise event
    GCTAEventAtom* event = NULL;

    // Get pointer on CTA pointing
    GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing());
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_MC);
    }

    // Get pointing direction zenith angle and azimuth [radians]
    double zenith  = pnt->zenith();
    double azimuth = pnt->azimuth();

    // Get radial offset and polar angles of true photon in camera [radians]
    double theta_p = pnt->dir().dist(photon.dir());
    double phi_p   = 0.0;  //TODO Implement Phi dependence

    // Compute effective area for photon
    double srcLogEng      = photon.energy().log10TeV();
    double effective_area = aeff(theta_p, phi_p, zenith, azimuth, srcLogEng);

    // Compute limiting value
    double ulimite = effective_area / area;

    // Continue only if event is detected
    if (ran.uniform() <= ulimite) {
        
        // Apply deadtime correction
        double deadc = obs.deadc(photon.time());
        if (deadc >= 1.0 || ran.uniform() <= deadc) {

            // Simulate offset from photon arrival direction
            //TODO: Make a proper implementation depending on the response
            // version. For now, the first Gaussian is used.
            GCTAPsfPars pars  = psf_dummy_sigma(srcLogEng, theta_p);
            double      theta = pars[1] * ran.chisq2() * rad2deg;
            double      phi   = 360.0 * ran.uniform();

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

        } // endif: detector was alive

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
 *
 * @todo Add support for other response versions.
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

    // Set PSF version
    m_psf_version = -10;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load CTA ARF vector
 *
 * @param[in] filename FITS file name.
 *
 * This method loads a CTA ARF vector from a FITS file. See the
 * GCTAResponse::read_arf method for more information on ARF vector reading.
 ***************************************************************************/
void GCTAResponse::load_arf(const std::string& filename)
{
    // Open ARF FITS file
    GFits file(filename);

    // Get ARF table
    GFitsTable* table = file.table("SPECRESP");

    // Read ARF
    read_arf(table);

    // Close ARF FITS file
    file.close();

    // Store filename
    m_arffile = filename;

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
 ***************************************************************************/
void GCTAResponse::load_psf(const std::string& filename)
{
    // Open PSF FITS file
    GFits file(filename);

    // Get PSF table. We assure here that the PSF information is stored
    // in extension number 1 (the second HDU of the FITS file)
    GFitsTable* table = file.table(1);

    // If the table has a single row we have a response table
    if (table->nrows() == 1) {

        // Read PSF table
        m_psf_table.read(table);

        // Set energy axis to logarithmic scale
        m_psf_table.axis_log10(0);

        // Set offset angle axis to radians
        m_psf_table.axis_radians(1);

        // Set PSF version
        m_psf_version = -8;
    }

    // ... otherwise we have a response vector
    else {
    
        // Read PSF vector
        read_psf(table);
        
        // Set PSF version
        m_psf_version = -9;
    }

    // Close PSF FITS file
    file.close();

    // Store filename
    m_psffile = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA ARF vector
 *
 * @param[in] hdu FITS table pointer.
 *
 * This method reads a CTA ARF vector from the FITS HDU. Note that the
 * energies are converted to TeV and the effective area is converted to cm2.
 * Conversion is done based on the units provided for the energy and
 * effective area columns. Units that are recognized are 'keV', 'MeV', 'GeV',
 * 'TeV', 'm^2', 'm2', 'cm^2' and 'cm^2' (case independent).
 *
 * @todo Assign appropriate theta angle for PSF. So far we use onaxis.
 *       For appropriate theta angle assignment, we would need this
 *       information in the response header.
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

    // Determine unit conversion factors (default: TeV and cm^2)
    std::string u_energy_lo = tolower(strip_whitespace(energy_lo->unit()));
    std::string u_energy_hi = tolower(strip_whitespace(energy_hi->unit()));
    std::string u_specresp  = tolower(strip_whitespace(specresp->unit()));
    double c_energy_lo = 1.0;
    double c_energy_hi = 1.0;
    double c_specresp  = 1.0;
    if (u_energy_lo == "kev") {
        c_energy_lo = 1.0e-9;
    }
    else if (u_energy_lo == "mev") {
        c_energy_lo = 1.0e-6;
    }
    else if (u_energy_lo == "gev") {
        c_energy_lo = 1.0e-3;
    }
    if (u_energy_hi == "kev") {
        c_energy_hi = 1.0e-9;
    }
    else if (u_energy_hi == "mev") {
        c_energy_hi = 1.0e-6;
    }
    else if (u_energy_hi == "gev") {
        c_energy_hi = 1.0e-3;
    }
    if (u_specresp == "m^2" || u_specresp == "m2") {
        c_specresp = 10000.0;
    }

    // Extract number of energy bins
    int num = energy_lo->length();

    // Set nodes
    for (int i = 0; i < num; ++i) {
    
        // Compute log10 mean energy in TeV
        double e_min = energy_lo->real(i) * c_energy_lo;
        double e_max = energy_hi->real(i) * c_energy_hi;
        double logE  = 0.5 * (log10(e_min) + log10(e_max));

        // Initialise scale factor
        double scale = m_arf_scale;

        // Optionally compute scaling factor from thetacut. This is done
        // by computing the containment fraction for the specified thetacut.
        if (m_arf_thetacut > 0.0) {

            // Get PSF parameters for node energy and theta angle
            //TODO: Implement theta angle computation
            GCTAPsfPars pars = psf_dummy_sigma(logE, 0.0);

            // Get maximum integration radius
            double rmax = m_arf_thetacut * deg2rad;

            // Setup integration kernel
            cta_npsf_kern_rad_azsym integrand(this,
                                              rmax,
                                              0.0,
                                              pars);

            // Setup integration
            GIntegral integral(&integrand);
            integral.eps(m_eps);

            // Perform integration
            double fraction = integral.romb(0.0, rmax);

            // Update scale factor
            if (fraction > 0.0) {
                scale /= fraction;
                #if defined(G_DEBUG_READ_ARF)
                std::cout << "GCTAResponse::read_arf:";
                std::cout << " e_min=" << e_min;
                std::cout << " e_max=" << e_max;
                std::cout << " logE=" << logE;
                std::cout << " scale=" << scale;
                std::cout << " fraction=" << fraction;
                std::cout << std::endl;
                #endif
            }
            else {
                std::cout << "WARNING: GCTAResponse::read_arf:";
                std::cout << " Non-positive integral occured in";
                std::cout << " PSF integration in GCTAResponse::read_arf.";
                std::cout << std::endl;
            }
        }

        // Compute effective area in cm2
        double aeff = specresp->real(i) * c_specresp * scale;
        
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
 * @brief Read CTA PSF vector
 *
 * @param[in] hdu FITS table pointer.
 *
 * This method reads a CTA PSF vector from the FITS HDU. Note that the
 * energies are converted to TeV. Conversion is done based on the units
 * provided for the energy columns. Units that are recognized are 'keV',
 * 'MeV', 'GeV', and 'TeV' (case independent).
 *
 * The Gaussian width parameter may be either given in 68% containment radius
 * (R68) or in sigma (ANGRES40; corresponding to a 38% containment radius).
 * The former format is used for CTA and H.E.S.S. data, the latter for
 * MAGIC data. All these things should be more uniform once we have a well
 * defined format.
 ***************************************************************************/
void GCTAResponse::read_psf(const GFitsTable* hdu)
{
    // Clear arrays
    m_psf_logE.clear();
    m_r68.clear();

    // Get pointers to table columns
    const GFitsTableCol* energy_lo = &(*hdu)["ENERG_LO"];
    const GFitsTableCol* energy_hi = &(*hdu)["ENERG_HI"];

    // Handle various data formats (H.E.S.S. and MAGIC)
    const GFitsTableCol* r68;
    double               r68_scale = 1.0;
    if (hdu->hascolumn("R68")) {
        r68 = &(*hdu)["R68"];
    }
    else {
        r68 = &(*hdu)["ANGRES40"];
        r68_scale = 1.0 / 0.6624305; // MAGIC PSF is already 1 sigma
    }

    // Determine unit conversion factors (default: TeV)
    std::string u_energy_lo = tolower(strip_whitespace(energy_lo->unit()));
    std::string u_energy_hi = tolower(strip_whitespace(energy_hi->unit()));
    double c_energy_lo = 1.0;
    double c_energy_hi = 1.0;
    if (u_energy_lo == "kev") {
        c_energy_lo = 1.0e-9;
    }
    else if (u_energy_lo == "mev") {
        c_energy_lo = 1.0e-6;
    }
    else if (u_energy_lo == "gev") {
        c_energy_lo = 1.0e-3;
    }
    if (u_energy_hi == "kev") {
        c_energy_hi = 1.0e-9;
    }
    else if (u_energy_hi == "mev") {
        c_energy_hi = 1.0e-6;
    }
    else if (u_energy_hi == "gev") {
        c_energy_hi = 1.0e-3;
    }

    // Extract number of energy bins
    int num = energy_lo->length();

    // Set nodes
    for (int i = 0; i < num; ++i) {
    
        // Compute log10 mean energy in TeV
        double e_min = energy_lo->real(i) * c_energy_lo;
        double e_max = energy_hi->real(i) * c_energy_hi;
        double logE  = 0.5 * (log10(e_min) + log10(e_max));
        
        // Extract r68 value and scale as required
        double r68_value = r68->real(i) * r68_scale;
        
        // Store log10 mean energy and r68 value
        m_psf_logE.append(logE);
        m_r68.push_back(r68_value);

    } // endfor: looped over nodes
    
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
    result.append("\n"+parformat("ARF file name")+m_arffile);
    result.append("\n"+parformat("RMF file name")+m_rmffile);
    result.append("\n"+parformat("PSF file name")+m_psffile);
    result.append("\n"+parformat("ARF theta angle cut"));
    if (m_arf_thetacut > 0) {
        result.append(str(m_arf_thetacut));
    }
    else {
        result.append("none");
    }
    result.append("\n"+parformat("ARF scaling")+str(m_arf_scale));
    if (m_offset_sigma == 0) {
        result.append("\n"+parformat("Offset angle dependence")+"none");
    }
    else {
        std::string txt = "Fixed sigma="+str(m_offset_sigma);
        result.append("\n"+parformat("Offset angle dependence")+txt);
    }
    result.append("\n"+parformat("Effective area nodes")+str(m_aeff_logE.size()));
    
    // Append PSF information
    result.append("\n"+parformat("PSF version"));
    if (m_psf_version == -8) {
        result.append("Response table");
    }
    else if (m_psf_version == -9) {
        result.append("Response vector");
    }
    else if (m_psf_version == -10) {
        result.append("ASCII table");
    }
    else {
        result.append("Unknown version");
    }
    result.append(" ("+str(m_psf_version)+")");
    if (m_psf_version == -8) {
        result += "\n" + m_psf_table.print();
    }
    else {
        result.append("\n"+parformat("PSF nodes")+str(m_psf_logE.size()));
    }

    // Debug option: Plot Aeff
    #if defined(G_DEBUG_PRINT_AEFF)
    result.append("\n"+parformat("Effective area"));
    for (int i = 0; i < m_aeff_logE.size(); ++i) {
        result.append("\n"+parformat("logE="+str(m_aeff_logE[i])));
        result.append("Aeff="+str(m_aeff.at(i))+" m2");
    }
    #endif

    // Debug option: Plot PSF
    #if defined(G_DEBUG_PRINT_PSF)
    result.append("\n"+parformat("Point spread function"));
    for (int i = 0; i < m_psf_logE.size(); ++i) {
        result.append("\n"+parformat("logE="+str(m_psf_logE[i])));
        result.append("r68="+str(m_r68.at(i))+" deg");
    }
    #endif

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
 *
 * Note that this method approximates the true theta angle (angle between
 * incident photon and pointing direction) by the measured theta angle
 * (angle between the measured photon arrival direction and the pointing
 * direction). Given the slow variation of the PSF shape over the field of
 * view, this approximation should be fine. It helps in fact a lot in
 * speeding up the computations.
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
    if (ctaobs == NULL) {
        throw GCTAException::bad_observation_type(G_IRF_EXTENDED);
    }

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing();
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_IRF_EXTENDED);
    }

    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&obsDir);
    if (dir == NULL) {
        throw GCTAException::bad_instdir_type(G_IRF_EXTENDED);
    }

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

    // Assign the observed theta angle (eta) as the true theta angle
    // between the source and the pointing directions. This is a (not
    // too bad) approximation which helps to speed up computations.
    // If we want to do this correctly, however, we would need to move
    // the psf_dummy_sigma down to the integration kernel, and we would
    // need to make sure that psf_delta_max really gives the absolute
    // maximum (this is certainly less critical)
    double srcTheta = eta;

    // Get PSF parameters.
    GCTAPsfPars psf_parameters = psf_dummy_sigma(srcLogEng, srcTheta);

    // Get maximum PSF and source radius in radians.
    double delta_max = psf_delta_max(srcTheta, 0.0, 0.0, 0.0, srcLogEng);
    double src_max   = model.radial()->theta_max();

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
                                          model.radial(),
                                          pnt->zenith(),
                                          pnt->azimuth(),
                                          srcLogEng,
                                          obsLogEng,
                                          psf_parameters,
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
    std::cout << "GCTAResponse::irf_extended:";
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
 * @param[in] obs Observation.
 *
 * Performs integration of the spatial component of the diffuse source model
 * multiplied by the IRF over the true photon arrival directions for a given
 * observed arrival direction.
 *
 * The integration is performed in the reference of the observed arrival
 * direction. Integration is done first over the azimuth angle \f$\phi\f$ and
 * then over the offset angle \f$\theta\f$. Specifically, the method computes
 *
 * \f[irf = \int_{0}^{\theta_{\rm max}}
 *          \int_{0}^{2\pi}
 *          M(\theta, \phi) IRF(\theta, \phi) d\phi \sin \theta d\theta\f]
 *
 * Here, \f$M(\theta, \phi)\f$ is the spatial component of the diffuse source
 * model and \f$IRF(\theta, \phi)\f$ is the point spread function.
 *
 * The integration kernels for this method are implemented by the response
 * helper classes cta_irf_diffuse_kern_theta and cta_irf_diffuse_kern_phi.
 ***************************************************************************/
double GCTAResponse::irf_diffuse(const GInstDir&            obsDir,
                                 const GEnergy&             obsEng,
                                 const GTime&               obsTime,
                                 const GModelDiffuseSource& model,
                                 const GEnergy&             srcEng,
                                 const GTime&               srcTime,
                                 const GObservation&        obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        throw GCTAException::bad_observation_type(G_IRF_DIFFUSE);
    }

    // Get pointer on CTA pointing
    const GCTAPointing *pnt = ctaobs->pointing();
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_IRF_DIFFUSE);
    }

    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&obsDir);
    if (dir == NULL) {
        throw GCTAException::bad_instdir_type(G_IRF_DIFFUSE);
    }

    // Determine angular distance between measured photon direction and
    // pointing direction [radians]
    double eta = pnt->dir().dist(dir->skydir());

    // Get log10(E/TeV) of true and measured photon energies
    double srcLogEng = srcEng.log10TeV();
    double obsLogEng = obsEng.log10TeV();

    // Assign the observed theta angle (eta) as the true theta angle
    // between the source and the pointing directions. This is a (not
    // too bad) approximation which helps to speed up computations.
    // If we want to do this correctly, however, we would need to move
    // the psf_dummy_sigma down to the integration kernel, and we would
    // need to make sure that psf_delta_max really gives the absolute
    // maximum (this is certainly less critical)
    double srcTheta = eta;

    // Get PSF parameters.
    GCTAPsfPars pars = psf_dummy_sigma(srcLogEng, srcTheta);

    // Get maximum PSF radius in radians
    double delta_max = psf_delta_max(srcTheta, 0.0, 0.0, 0.0, srcLogEng);

    // Initialise IRF value
    double irf = 0.0;

    // Perform zenith angle integration if interval is valid
    if (delta_max > 0.0) {

        // Compute rotation matrix to convert from coordinates (theta,phi)
        // in the reference frame of the observed arrival direction into
        // celestial coordinates
        GMatrix ry;
        GMatrix rz;
        GMatrix rot;
        ry.eulery(dir->dec_deg() - 90.0);
        rz.eulerz(-dir->ra_deg());
        rot = transpose(ry * rz);

        // Setup integration kernel
        cta_irf_diffuse_kern_theta integrand(this,
                                             model.spatial(),
                                             pnt->zenith(),
                                             pnt->azimuth(),
                                             srcLogEng,
                                             obsLogEng,
                                             &rot,
                                             pars,
                                             eta);

        // Integrate over zenith angle
        GIntegral integral(&integrand);
        integral.eps(1.0e-2);
        irf = integral.romb(0.0, delta_max);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(irf) || isinfinite(irf)) {
            std::cout << "*** ERROR: GCTAResponse::irf_diffuse:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", delta_max=" << delta_max << ")";
            std::cout << std::endl;
        }
        #endif
    }

    // Compile option: Show integration results
    #if defined(G_DEBUG_IRF_DIFFUSE)
    std::cout << "GCTAResponse::irf_diffuse:";
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
 * @brief Return spatial integral of extended source model
 *
 * @param[in] model Extended source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * Note that we estimate the integration radius based on the size of the
 * onaxis PSF in this method. This should be fine as long as the offaxis
 * PSF is not considerably larger than the onaxis PSF. We should verify
 * this, however.
 *
 * @todo Verify that offaxis PSF is not considerably larger than onaxis
 *       PSF. 
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
    if (ctaobs == NULL) {
        throw GCTAException::bad_observation_type(G_NPRED_EXTENDED);
    }

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(ctaobs->events());
    if (events == NULL) {
        throw GException::no_list(G_NPRED_EXTENDED);
    }

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius (radians). We do this for the onaxis PSF only,
    // as this allows us doing this computation in the outer loop. This
    // should be sufficient here, unless the offaxis PSF becomes much worse
    // than the onaxis PSF. In this case, we may add a safety factor here
    // to make sure we encompass the entire PSF.
    GCTAPsfPars psf_parameters = psf_dummy_sigma(srcLogEng, 0.0);
    double      psf_max_radius = psf_dummy_max(psf_parameters);

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
        cta_npred_radial_kern_theta integrand(this,
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
        std::cout << "GCTAResponse::npred_extended:";
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


/***********************************************************************//**
 * @brief Return spatial integral of diffuse source model
 *
 * @param[in] model Diffuse source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * This method provides the intergal of the diffuse source model over the
 * region of interest. The integration is performed in the ROI system.
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
double GCTAResponse::npred_diffuse(const GModelDiffuseSource& model,
                                   const GEnergy&             srcEng,
                                   const GTime&               srcTime,
                                   const GObservation&        obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        throw GCTAException::bad_observation_type(G_NPRED_DIFFUSE);
    }

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(ctaobs->events());
    if (events == NULL) {
        throw GException::no_list(G_NPRED_DIFFUSE);
    }

    // Get log10(E/TeV) of true photon energy
    double srcLogEng = srcEng.log10TeV();

    // Get maximum PSF radius (radians). We do this for the onaxis PSF only,
    // as this allows us doing this computation in the outer loop. This
    // should be sufficient here, unless the offaxis PSF becomes much worse
    // than the onaxis PSF. In this case, we may add a safety factor here
    // to make sure we encompass the entire PSF.
    GCTAPsfPars psf_parameters = psf_dummy_sigma(srcLogEng, 0.0);
    double      psf_max_radius = psf_dummy_max(psf_parameters);

    // Extract ROI radius (radians)
    double roi_radius = events->roi().radius() * deg2rad;

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
        GMatrix rot;
        ry.eulery(events->roi().centre().dec_deg() - 90.0);
        rz.eulerz(-events->roi().centre().ra_deg());
        rot = transpose(ry * rz);
        
        // Setup integration kernel
        cta_npred_diffuse_kern_theta integrand(this,
                                               model.spatial(),
                                               &srcEng,
                                               &srcTime,
                                               ctaobs,
                                               &rot);

        // Integrate over theta
        GIntegral integral(&integrand);
        integral.eps(1.0e-4);
        npred = integral.romb(0.0, roi_psf_radius);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NPRED_DIFFUSE)
        std::cout << "GCTAResponse::npred_diffuse:";
        std::cout << " srcEng=" << srcEng;
        std::cout << " srcTime=" << srcTime;
        std::cout << " roi_psf_radius=" << roi_psf_radius;
        std::cout << " npred=" << npred << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GCTAResponse::npred_diffuse:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", roi_psf_radius=" << roi_psf_radius;
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
    if (aeff < 0) {
        aeff = 0.0;
    }
    
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
 * @todo So far the parameters phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::psf(const double& delta,
                         const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const double& srcLogEng) const
{
    // Determine energy dependent PSF parameters
    GCTAPsfPars psf_parameters = psf_dummy_sigma(srcLogEng, theta);

    // Compute PSF
    double psf = psf_dummy(delta, psf_parameters);

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
 * @todo So far the parameters phi, zenith, and azimuth are not used.
 ***************************************************************************/
double GCTAResponse::psf_delta_max(const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth,
                                   const double& srcLogEng) const
{
    // Determine energy dependent width of PSF
    GCTAPsfPars psf_parameters = psf_dummy_sigma(srcLogEng, theta);

    // Set maximum angular separation
    double delta_max = psf_dummy_max(psf_parameters);

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
 *       (see comment about kluge below)
 ***************************************************************************/
double GCTAResponse::npsf(const GSkyDir&      srcDir,
                          const double&       srcLogEng,
                          const GTime&        srcTime,
                          const GCTAPointing& pnt,
                          const GCTARoi&      roi) const
{
    // Declare result
    double value = 0.0;
    
    // Compute offset angle of source direction in camera system
    double srcTheta = pnt.dir().dist(srcDir);

    // Extract relevant parameters from arguments
    double      roi_radius       = roi.radius() * deg2rad;
    double      roi_psf_distance = roi.centre().dist(srcDir);
    GCTAPsfPars psf_parameters   = psf_dummy_sigma(srcLogEng, srcTheta);
    double      rmax             = psf_dummy_max(psf_parameters);

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
                                              psf_parameters);

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
 * @param[in] pars PSF parameters.
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
 * @todo Correct method documentation.
 ***************************************************************************/
double GCTAResponse::psf_dummy(const double& delta, const GCTAPsfPars& pars) const
{
    // Initialise response value
    double value = 0.0;

    // Case A: Response is stored in a response table
    if (m_psf_version == -8) {

        // Compute distance squared
        double delta2 = delta * delta;
        
        // Compute Psf value
        value  = exp(pars[6] * delta2);
        value += exp(pars[7] * delta2) * pars[2];
        value += exp(pars[8] * delta2) * pars[4];
        value *= pars[0];
    }
    
    // Case B: Response is stored in a response vector
    else {
        
        // Compute Psf value
        value = pars[0] * exp(pars[2] * delta * delta);
    }

    // Return PSF value
    return value;
}


/***********************************************************************//**
 * @brief Return PSF parameter vector
 *
 * @param[in] srcLogEng Log10 of true photon energy (E/TeV).
 * @param[in] srcTheta  Offset angle of source in camera system (radians).
 *
 * This method returns the PSF parameter vector as function of incident
 * photon energy and offset angle. Two response types are supported.
 *
 * If the response is stored in a one dimensional vector, the PSF parameter
 * vector is composed of 3 parameters:
 * \f[{\tt SCALE} = \frac{1}{2 \pi \sigma^2}\f]
 * \f[{\tt SIGMA} = \sigma\f]
 * \f[{\tt WIDTH} = -\frac{1}{2 \sigma^2}\f]
 * where
 * \f$\sigma\f$ is the Gaussian sigma in radians.
 *
 * If the response is stored in a 2D response table, the PSF parameter
 * vector is composed of 9 parameters
 * \f[{\tt SCALE}\f]
 * \f[{\tt SIGMA\_1} = \sigma_1\f]
 * \f[{\tt AMPL\_2} = a_2\f]
 * \f[{\tt SIGMA\_2} = \sigma_2\f]
 * \f[{\tt AMPL\_3} = a_3\f]
 * \f[{\tt SIGMA\_3} = \sigma_3\f]
 * \f[{\tt WIDTH\_1} = -\frac{1}{2 \sigma_1^2}\f]
 * \f[{\tt WIDTH\_2} = -\frac{1}{2 \sigma_2^2}\f]
 * \f[{\tt WIDTH\_3} = -\frac{1}{2 \sigma_3^2}\f]
 * where
 * \f$\sigma_i\f$ are the Gaussian sigma values of the three components
 * in radians and
 * \f$a_i\f$ are the relative amplitudes of the 2nd and 3rd Gaussians (the
 * relative amplitude of the first Gaussian is by definition 1).
 *
 * @todo Convert sigma parameter to radians upon loading of the response
 *       table. This saves some operations here.
 * @todo Verify parameter validity. Interpolation may lead for example to
 *       negative sigmas or scales and amplitudes. This should be checked
 *       and avoided (if possible). We need to develop a strategy to deal
 *       with such cases.
 ***************************************************************************/
GCTAPsfPars GCTAResponse::psf_dummy_sigma(const double& srcLogEng,
                                          const double& srcTheta) const
{
    // Allocate PSF parameter vector
    GCTAPsfPars pars;
    
    // Case A: Response is stored in a response table
    if (m_psf_version == -8) {
    
        // Interpolate response parameters
        pars = m_psf_table(srcLogEng, srcTheta);

        // Convert sigma parameters to radians
        pars[1] *= deg2rad;
        pars[3] *= deg2rad;
        pars[5] *= deg2rad;

        // Compute normalization scale
        double sigma1   = pars[1] * pars[1];
        double sigma2   = pars[3] * pars[3];
        double sigma3   = pars[5] * pars[5];
        double integral = twopi * (sigma1 +
                                   sigma2 * pars[2] +
                                   sigma3 * pars[4]);
        double scale = (integral > 0.0) ? 1.0 / integral : 0.0;

        // Compile option: Show scaling results
        #if defined(G_DEBUG_PSF_DUMMY_SIGMA)
        std::cout << "GCTAResponse::psf_dummy_sigma:";
        std::cout << " srcLogEng=" << srcLogEng;
        std::cout << " srcTheta=" << srcTheta;
        std::cout << " old scale=" << pars[0];
        std::cout << " new scale=" << scale << std::endl;
        #endif

        // Update scale
        pars[0] = scale;

        // Compute widths
        double width1 = -0.5 / sigma1;
        double width2 = -0.5 / sigma2;
        double width3 = -0.5 / sigma3;

        // Append widths
        pars.push_back(width1);
        pars.push_back(width2);
        pars.push_back(width3);
    }
    
    // Case B: Response is stored in a response vector
    else {

        // Set conversion factor from 68% containment radius to 1 sigma
        const double conv = 0.6624305 * deg2rad;

        // Determine Gaussian sigma in radians
        double sigma = m_psf_logE.interpolate(srcLogEng, m_r68) * conv;
        
        // Derive width=-0.5/(sigma*sigma) and scale=1/(twopi*sigma*sigma)
        double sigma2 = sigma * sigma;
        double scale  =  1.0 / (twopi * sigma2);
        double width  = -0.5 / sigma2;
        
        // Store PSF parameters in vector (3 elements)
        pars.reserve(3);
        pars.push_back(scale);
        pars.push_back(sigma);
        pars.push_back(width);

    }

    // Return result
    return pars;
}


/***********************************************************************//**
 * @brief Returns radius beyond which PSF is negligible (in radians)
 *
 * @param[in] pars PSF parameters.
 *
 * Determine the radius beyond which the PSF becomes negligible. This radius
 * is set by this method to \f$5 \times \sigma\f$, where \f$\sigma\f$ is the
 * Gaussian width of the largest PSF component.
 *
 * The method supports both response vectors and response tables. For a
 * response vector, the Gaussian width \f$\sigma\f$ is used for the
 * computation. For the response table, which is composed of 3 Gaussians,
 * the maximum
 * \f[\sigma = \max{ \sigma1, \sigma_2, \sigma_3 }\f]
 * is used.
 *
 * For a definition of the PSF parameter vector refer to the description of
 * the GCTAResponse::psf_dummy_sigma method.
 ***************************************************************************/
double GCTAResponse::psf_dummy_max(const GCTAPsfPars& pars) const
{
    // Initialise size of PSF
    double sigma = 0.0;
    
    // Case A: Response is stored in a response table
    if (m_psf_version == -8) {
        sigma = pars[1];
        if (pars[3] > sigma) sigma = pars[3];
        if (pars[5] > sigma) sigma = pars[5];
    }

    // Case B: Response is stored in a response vector
    else {
        sigma = pars[1];
    }

    // Return radius
    return (5.0 * sigma);
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
    m_arf_thetacut = 0.0;    // Default: no thetacut
    m_arf_scale    = 1.0;    // Default: ARF is unscaled

    // Initialise PSF information
    m_psf_version = -99; // Unknown PSF version
    m_psf_table.clear();
    
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
    m_arf_thetacut = rsp.m_arf_thetacut;
    m_arf_scale    = rsp.m_arf_scale;

    // Copy PSF information
    m_psf_version = rsp.m_psf_version;
    m_psf_table   = rsp.m_psf_table;

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
    if (fptr == NULL) {
        throw GCTAException::file_open_error(G_READ, filename);
    }

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements. Strip empty elements from vector.
        std::vector<std::string> elements = split(line, " ");
        for (int i = elements.size()-1; i >= 0; i--) {
            if (strip_whitespace(elements[i]).length() == 0) {
                elements.erase(elements.begin()+i);
            }
        }

        // Skip header
        if (elements[0].find("log(E)") != std::string::npos) {
            continue;
        }

        // Break loop if end of data table has been reached
        if (elements[0].find("----------") != std::string::npos) {
            break;
        }

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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
