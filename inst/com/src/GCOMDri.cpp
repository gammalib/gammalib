/***************************************************************************
 *                  GCOMDri.cpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMDri.cpp
 * @brief COMPTEL Data Space class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GWcs.hpp"
#include "GFits.hpp"
#include "GFitsImage.hpp"
#include "GCOMDri.hpp"
#include "GCOMOads.hpp"
#include "GCOMTim.hpp"
#include "GCOMStatus.hpp"
#include "GCOMSupport.hpp"
#include "GCOMEventList.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_DRE    "GCOMDri::dre(GCOMEventList&, GCOMOads&, GCOMTim&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_CHECK_EHA_COMPUTATION

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_DRE
#define G_DEBUG_DRG


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCOMDri::GCOMDri(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File name constructor
 *
 * @param[in] filename COMPTEL Data Space FITS file name.
 ***************************************************************************/
GCOMDri::GCOMDri(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load DRI FITS file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dri COMPTEL Data Space.
 ***************************************************************************/
GCOMDri::GCOMDri(const GCOMDri& dri)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dri);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMDri::~GCOMDri(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] dri COMPTEL Data Space.
 * @return COMPTEL Data Space.
 ***************************************************************************/
GCOMDri& GCOMDri::operator=(const GCOMDri& dri)
{
    // Execute only if object is not identical
    if (this != &dri) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dri);

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
 * @brief Clear COMPTEL Data Space
 ***************************************************************************/
void GCOMDri::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Data Space
 *
 * @return Pointer to deep copy of COMPTEL Data Space.
 ***************************************************************************/
GCOMDri* GCOMDri::clone(void) const
{
    return new GCOMDri(*this);
}


/***********************************************************************//**
 * @brief Compute event cube
 *
 * @param[in] events Events list.
 * @param[in] oads Orbit Aspect Data.
 * @param[in] tim Good Time Intervals.
 * @param[in] zeta Minimum Earth horizon - Phibar cut (deg).
 *
 * @exception GException::invalid_argument
 *            DRE cube has a non-positive Phibar bin size.
 *
 * Compute DRE event cube from event list (EVP), Good Time Intervals (TIM)
 * and Orbit Aspect Data (OAD).
 ***************************************************************************/
void GCOMDri::dre(const GCOMEventList& events,
                  const GCOMOads&      oads,
                  const GCOMTim&       tim,
                  const double&        zeta)
{
    // Debug
    #if defined(G_DEBUG_DRE)
    std::cout << "GCOMDri::dre" << std::endl;
    std::cout << "============" << std::endl;
    #endif

    // Check for positive Phibar bin size
    if (m_phibin <= 0.0) {
        std::string msg = "DRE cube has a non-positive Phibar bin size. Please "
                          "specify a DRE cube with a positive Phibar bin size.";
        throw GException::invalid_argument(G_DRE, msg);
    }

    // Initialise variables
    int        i_evt = 0;
    GTime      tstart;
    GTime      tstop;
    GCOMStatus status;

    // Initialise statistics
    int num_used_superpackets     = 0;
    int num_skipped_superpackets  = 0;
    int num_used_events           = 0;
    int num_event_outside_sp      = 0;
    int num_energy_too_low        = 0;
    int num_energy_too_high       = 0;
    int num_no_scatter_angle      = 0;
    int num_bad_modcom            = 0;
    int num_eha_too_small         = 0;
    int num_phibar_too_low        = 0;
    int num_phibar_too_high       = 0;
    int num_outside_dre           = 0;
    int num_outside_e1            = 0;
    int num_outside_e2            = 0;
    int num_outside_tof           = 0;
    int num_outside_psd           = 0;
    int num_bad_reflag            = 0;
    int num_event_before_dre      = 0;
    int num_event_after_dre       = 0;
    int num_d1module_off          = 0;
    int num_d2module_off          = 0;
    int num_processed             = 0;
    int num_superpackets          = 0;
    int last_superpackets         = 0;
    int last_skipped_superpackets = 0;

    // Set all DRE bins to zero
    for (int i = 0; i < size(); ++i) {
        (*this)[i] = 0.0;
    }

    // Signal that loop should be terminated
    bool terminate = false;

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < oads.size(); ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = oads[i_oad];

        // Increment superpacket counter
        num_superpackets++;

        // Skip superpacket if it is not within Good Time Intervals.
        // According to the original routine dalsel17.p7time.f only the
        // superpacket start time is checked.
        if (!tim.contains(oad.tstart())) {
            num_skipped_superpackets++;
            continue;
        }

        // Update superpackets statistics
        num_used_superpackets++;

        // Prepare Earth horizon angle comparison
        GSkyDir sky_geocentre;
        double  theta_geocentre = double(oad.gcel());
        double  phi_geocentre   = double(oad.gcaz());
        sky_geocentre.radec_deg(phi_geocentre, 90.0-theta_geocentre);

        // Update validity interval so that DRE will have the correct
        // interval
        if (num_used_superpackets == 1) {
            tstart = oad.tstart();
            tstop  = oad.tstop();
        }
        else {
            tstop  = oad.tstop();
        }

        // Collect all events within superpacket. Break if the end
        // of the event list was reached.
        for (; i_evt < events.size(); ++i_evt) {

            // Get pointer to event
            const GCOMEventAtom* event = events[i_evt];

            // Break loop if the end of the superpacket was reached
            if (event->time() > oad.tstop()) {
                break;
            }

            // Increase number of processed events
            num_processed++;

            // Skip event if it lies before the DRE start
            if (event->time() < m_gti.tstart()) {
                num_event_before_dre++;
                continue;
            }

            // Break if event lies after the DRE stop
            if (event->time() > m_gti.tstop()) {
                num_event_after_dre++;
                terminate = true;
                break;
            }

            // Skip event if it lies before the superpacket start
            if (event->time() < oad.tstart()) {
                num_event_outside_sp++;
                continue;
            }

            // Skip event if it lies outside energy range
            if (event->energy() < m_ebounds.emin()) {
                num_energy_too_low++;
                continue;
            }
            else if (event->energy() > m_ebounds.emax()) {
                num_energy_too_high++;
                continue;
            }

            // Apply event selection, see PSSEVP.F
            // TODO: Hard coded for the moment, seems like EVP files
            //       contain already these thresholds !!!
            if (event->e1() < 0.070 || event->e1() > 20.0) {
                num_outside_e1++;
                continue;
            }
            if (event->e2() < 0.650 || event->e2() > 30.0) {
                num_outside_e2++;
                continue;
            }
            if (event->tof() < 115.0 || event->tof() > 130.0) {
                num_outside_tof++;
                continue;
            }
            if (event->psd() < 0.0 || event->psd() > 110.0) {
                num_outside_psd++;
                continue;
            }
            if (event->reflag() < 1) {
                num_bad_reflag++;
                continue;
            }

            // Check whether the event has a scatter angle determined.
            // This is signaled by a scatter angle -10e20 in radians.
            if (event->theta() < -1.0e3) {
                num_no_scatter_angle++;
                continue;
            }

            // Check for valid module IDs from MODCOM
            if (event->modcom() < 1 && event->modcom() > 98) {
                num_bad_modcom++;
                continue;
            }

            // Extract module IDs from MODCOM
            int id2 = (event->modcom()-1)/7 + 1;      // [1-14]
            int id1 =  event->modcom() - (id2-1) * 7; // [1-7]

            // Check module status. Skip all events where module is
            // signalled as not on.
            if (status.d1status(oad.tjd(), id1) != 1) {
                num_d1module_off++;
                continue;
            }
            if (status.d2status(oad.tjd(), id2) != 1) {
                num_d2module_off++;
                continue;
            }

            // Compute Compton scatter angle index. Skip if it's invalid.
            int iphibar = (event->phibar() - m_phimin) / m_phibin;
            if (iphibar < 0) {
                num_phibar_too_low++;
                continue;
            }
            else if (iphibar >= nphibar()) {
                num_phibar_too_high++;
                continue;
            }

            // Compute Earth horizon angle for comparison with EVP value
            GSkyDir sky_event;
            double  theta_event = double(event->theta());
            double  phi_event   = double(event->phi());
            sky_event.radec_deg(-phi_event, 90.0-theta_event);
            double eha = sky_geocentre.dist_deg(sky_event) - oad.georad();

            // Option: Earth horizon angle comparison
            #if defined(G_CHECK_EHA_COMPUTATION)
            if (std::abs(eha - event->eha()) > 1.5) {
                std::string msg = "Earth horizon angle from EVP dataset ("+
                                  gammalib::str(event->eha())+" deg) "
                                  "differs from Earth horizon angle "
                                  "computed from Orbit Aspect Data ("+
                                  gammalib::str(eha)+" deg). Use the EVP "
                                  "value.";
                gammalib::warning(G_DRE, msg);
            }
            #endif

            // Check for Earth horizon angle. The limit is fixed here to
            // the values of 5, 7, ..., 53 deg for 25 Phibar layers that
            // was usually used for COMPTEL.
            // zeta = 5.0
            double ehamin = double(iphibar) * m_phibin + zeta;
            if (event->eha() < ehamin) {
                num_eha_too_small++;
                continue;
            }

            // Now fill the DRE event array
            GSkyPixel pixel = m_dri.dir2pix(event->dir().dir());
            int       ichi  = int(pixel.x());
            int       ipsi  = int(pixel.y());
            if ((ichi >= 0) && (ichi < nchi()) &&
                (ipsi >= 0) && (ipsi < npsi())) {
                int inx       = ichi + (ipsi + iphibar * npsi()) * nchi();
                (*this)[inx] += 1.0;
                num_used_events++;
            }
            else {
                num_outside_dre++;
            }

        } // endfor: collected events

        // Break if termination was signalled or if there are no more events
        if (terminate || i_evt >= events.size()) {
            break;
        }

    } // endfor: looped over Orbit Aspect Data

    // Set Good Time interval for DRE
    m_gti = GGti(tstart, tstop);

    // Debug
    #if defined(G_DEBUG_DRE)
    std::cout << "Total number of superpackets .: " << num_superpackets << std::endl;
    std::cout << "Used superpackets ............: " << num_used_superpackets << std::endl;
    std::cout << "Skipped superpackets .........: " << num_skipped_superpackets << std::endl;
    std::cout << "Total number of events .......: " << events.size() << std::endl;
    std::cout << "Processed events .............: " << num_processed << std::endl;
    std::cout << "Used events ..................: " << num_used_events << std::endl;
    std::cout << "Events outside superpacket ...: " << num_event_outside_sp << std::endl;
    std::cout << "Events before DRE GTI ........: " << num_event_before_dre << std::endl;
    std::cout << "Events after DRE GTI .........: " << num_event_after_dre << std::endl;
    std::cout << "Events with D1 module off ....: " << num_d1module_off << std::endl;
    std::cout << "Events with D2 module off ....: " << num_d2module_off << std::endl;
    std::cout << "Energy too low ...............: " << num_energy_too_low << std::endl;
    std::cout << "Energy too high ..............: " << num_energy_too_high << std::endl;
    std::cout << "Phibar too low ...............: " << num_phibar_too_low << std::endl;
    std::cout << "Phibar too high ..............: " << num_phibar_too_high << std::endl;
    std::cout << "Outside D1 selection .........: " << num_outside_e1 << std::endl;
    std::cout << "Outside D2 selection .........: " << num_outside_e2 << std::endl;
    std::cout << "Outside TOF selection ........: " << num_outside_tof << std::endl;
    std::cout << "Outside PSD selection ........: " << num_outside_psd << std::endl;
    std::cout << "Outside rejection flag sel. ..: " << num_bad_reflag << std::endl;
    std::cout << "Outside DRE cube .............: " << num_outside_dre << std::endl;
    std::cout << "No scatter angle .............: " << num_no_scatter_angle << std::endl;
    std::cout << "Bad module combination .......: " << num_bad_modcom << std::endl;
    std::cout << "Earth horizon angle too small : " << num_eha_too_small << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute geometry cube
 *
 * @param[in] oads Orbit Aspect Data.
 * @param[in] tim Good Time Intervals.
 * @param[in] zeta Minimum Earth horizon - Phibar cut (deg).
 *
 * Compute DRG cube from Orbit Aspect Data (OAD) and Good Time Intervals
 * (TIM).
 ***************************************************************************/
void GCOMDri::drg(const GCOMOads& oads, const GCOMTim& tim, const double& zeta)
{
    // Debug
    #if defined(G_DEBUG_DRG)
    std::cout << "GCOMDri::drg" << std::endl;
    std::cout << "============" << std::endl;
    #endif

    // Initialise variables
    int        npix = nchi() * npsi(); // Number of pixels in (Chi, Psi)
    GTime      tstart;                 // Start time
    GTime      tstop;                  // Stop time
    GCOMStatus status;                 // D1 & D2 module status

    // Initialise superpacket statistics
    int num_sp         = 0;
    int num_sp_used    = 0;
    int num_sp_skipped = 0;

    // Set all DRG bins to zero
    for (int i = 0; i < size(); ++i) {
        (*this)[i] = 0.0;
    }

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < oads.size(); ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = oads[i_oad];

        // Increment superpacket counter
        num_sp++;

        // Skip superpacket if it is not within Good Time Intervals.
        // According to the original routine dalsel17.p7time.f only the
        // superpacket start time is checked.
        if (!tim.contains(oad.tstart())) {
            num_sp_skipped++;
            continue;
        }

        //TODO: Add GTI check to be compliant with DRE

        // Update superpackets statistics
        num_sp_used++;

        // Prepare Earth horizon angle computation. The celestial system
        // is reinterpreted as the COMPTEL coordinate system, where the
        // 90 degrees - zenith angle becomes the declination and the azimuth
        // angle becomes Right Ascension. This allows us later to use the
        // GSkyDir::dist method to compute the distance between the geocentre
        // and the telescope Z-axis.
        double  theta_geocentre = double(oad.gcel());
        double  phi_geocentre   = double(oad.gcaz());
        GSkyDir geocentre_comptel;
        geocentre_comptel.radec_deg(phi_geocentre, 90.0-theta_geocentre);

        // Update validity interval so that DRG will have the correct time
        // interval
        if (num_sp_used == 1) {
            tstart = oad.tstart();
            tstop  = oad.tstop();
        }
        else {
            tstop  = oad.tstop();
        }

        // Compute geometrical factors for all (Chi, Psi)
        for (int index = 0; index < npix; ++index) {

            // Get sky direction for (Chi, Psi)
            GSkyDir sky = m_dri.inx2dir(index);

            // Convert sky direction to COMPTEL coordinates
            double theta = oad.theta(sky);
            double phi   = oad.phi(sky);

            // Compute geometric factor summed over all D1, D2
            double geometry = compute_geometry(oad.tjd(), theta, phi, status);

            // Compute Earth horizon angle as the distance between the Earth
            // centre in COMPTEL coordinates and the (Chi, Psi) pixel in
            // COMPTEL coordinates minus the radius of the Earth.
            GSkyDir chipsi_comptel;
            chipsi_comptel.radec_deg(-phi, 90.0-theta);
            double eha = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

            // Loop over Phibar
            for (int iphibar = 0; iphibar < nphibar(); ++iphibar) {

                // Compute minimum Earth horizon angle for Phibar layer
                double ehamin = double(iphibar) * m_phibin + zeta;

                // Add up geometry if the Earth horizon angle is equal or
                // larger than the minimum
                if (eha >= ehamin) {
                    int inx       = index + iphibar * npix;
                    (*this)[inx] += geometry;
                }

            } // endfor: looped over Phibar

        } // endfor: looped over (Chi, Psi)

    } // endfor: looped over Orbit Aspect Data

    // Divide DRG by total number of superpackets
    if (num_sp > 0) {
        double norm = 1.0 / double(num_sp);
        for (int i = 0; i < size(); ++i) {
            (*this)[i] *= norm;
        }
    }

    // Clear energy boundaries for DRG since DRG does not depend on energy)
    m_ebounds.clear();

    // Set the Good Time interval for DRG
    m_gti = GGti(tstart, tstop);

    // Debug
    #if defined(G_DEBUG_DRG)
    std::cout << "Total number of superpackets .: " << num_sp << std::endl;
    std::cout << "Used superpackets ............: " << num_sp_used << std::endl;
    std::cout << "Skipped superpackets .........: " << num_sp_skipped << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load COMPTEL Data Space from DRI FITS file
 *
 * @param[in] filename DRI FITS file name.
 ***************************************************************************/
void GCOMDri::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get HDU (pointer is always valid)
    const GFitsImage& hdu = *fits.image(0);

    // Read DRI file
    read(hdu);

    // Close FITS file
    fits.close();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Save COMPTEL Data Space into DRI FITS file
 *
 * @param[in] filename DRI FITS file name.
 ***************************************************************************/
void GCOMDri::save(const GFilename& filename, const bool& clobber) const
{
    // Create FITS file
    GFits fits;

    // Write data space into FITS file
    write(fits, filename.extname(gammalib::extname_dri));

    // Save FITS file
    fits.saveto(filename, clobber);

    // Close FITS file
    fits.close();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL Data Space from DRI FITS image
 *
 * @param[in] image DRI FITS image.
 ***************************************************************************/
void GCOMDri::read(const GFitsImage& image)
{
    // Clear
    clear();

    // Read sky map
    m_dri.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    com_wcs_mer2car(m_dri);

    // Read attributes
    read_attributes(&image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL Data Space into FITS image
 *
 * @param[in] fits FITS file.
 * @param[in] extname Extension name.
 ***************************************************************************/
void GCOMDri::write(GFits& fits, const std::string& extname) const
{
    // Write sky map into FITS file
    GFitsHDU *image = m_dri.write(fits, extname);

    // Write DRI attributes
    write_attributes(image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL Data Space
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Data Space information.
 ***************************************************************************/
std::string GCOMDri::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute scatter angle dimensions
        double      chimin = 0.0;
        double      chimax = 0.0;
        double      chibin = 0.0;
        double      psimin = 0.0;
        double      psimax = 0.0;
        double      psibin = 0.0;
        const GWcs* wcs    = dynamic_cast<const GWcs*>(m_dri.projection());
        if (wcs != NULL) {
            chibin = wcs->cdelt(0);
            chimin = wcs->crval(0) - (wcs->crpix(0)-0.5) * chibin;
            chimax = chimin + m_dri.nx() * chibin;
            psibin = wcs->cdelt(1);
            psimin = wcs->crval(1) - (wcs->crpix(1)-0.5) * psibin;
            psimax = psimin + m_dri.ny() * psibin;
        }

        // Compute Phibar maximum
        double phimax = m_phimin + m_dri.nmaps() * m_phibin;

        // Append header
        result.append("=== GCOMDri ===");

        // Append Phibar information
        result.append("\n"+gammalib::parformat("Chi range"));
        result.append(gammalib::str(chimin)+" - ");
        result.append(gammalib::str(chimax)+" deg");
        result.append("\n"+gammalib::parformat("Chi bin size"));
        result.append(gammalib::str(chibin)+" deg");
        result.append("\n"+gammalib::parformat("Psi range"));
        result.append(gammalib::str(psimin)+" - ");
        result.append(gammalib::str(psimax)+" deg");
        result.append("\n"+gammalib::parformat("Psi bin size"));
        result.append(gammalib::str(psibin)+" deg");
        result.append("\n"+gammalib::parformat("Phibar range"));
        result.append(gammalib::str(m_phimin)+" - ");
        result.append(gammalib::str(phimax)+" deg");
        result.append("\n"+gammalib::parformat("Phibar bin size"));
        result.append(gammalib::str(m_phibin)+" deg");

        // Append energy boundaries
        result.append("\n"+m_ebounds.print(gammalib::reduce(chatter)));

        // Append GTI
        result.append("\n"+m_gti.print(gammalib::reduce(chatter)));

        // EXPLICIT: Append sky map
        if (chatter >= EXPLICIT) {
            result.append("\n"+m_dri.print(gammalib::reduce(chatter)));
        }

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
void GCOMDri::init_members(void)
{
    // Initialise members
    m_dri.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_phimin = 0.0;
    m_phibin = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dri COMPTEL Data Space.
 ***************************************************************************/
void GCOMDri::copy_members(const GCOMDri& dri)
{
    // Copy members
    m_dri     = dri.m_dri;
    m_ebounds = dri.m_ebounds;
    m_gti     = dri.m_gti;
    m_phimin  = dri.m_phimin;
    m_phibin  = dri.m_phibin;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMDri::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read DRI attributes from FITS HDU
 *
 * @param[in] hdu FITS HDU pointer.
 ***************************************************************************/
void GCOMDri::read_attributes(const GFitsHDU* hdu)
{
    // Get phibar attributes
    m_phibin = hdu->real("CDELT3");
    m_phimin = hdu->real("CRVAL3") - (hdu->real("CRPIX3")-0.5) * m_phibin;
    
    // Get time attributes
    GTime tstart = com_time(hdu->integer("VISDAY"), hdu->integer("VISTIM"));
    GTime tstop  = com_time(hdu->integer("VIEDAY"), hdu->integer("VIETIM"));

    // Set Good Time Intervals
    m_gti = GGti(tstart, tstop);

    // Optionally read energy attributes
    if (hdu->has_card("E_MIN") && hdu->has_card("E_MAX")) {

        // Get energy attributes
        GEnergy emin(hdu->real("E_MIN"), "MeV");
        GEnergy emax(hdu->real("E_MAX"), "MeV");

        // Set energy boundaries
        m_ebounds = GEbounds(emin, emax);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write DRI attributes into FITS HDU
 *
 * @param[in] hdu FITS HDU pointer.
 ***************************************************************************/
void GCOMDri::write_attributes(GFitsHDU* hdu) const
{
    // Set Phibar keywords
    double crval3 = m_phimin + 0.5 * m_phibin;

    // Write Phibar keywords
    hdu->card("CTYPE3", "Phibar", "Compton scatter angle");
    hdu->card("CRPIX3", 1.0, "Pixel coordinate of reference point (starting from 1)");
    hdu->card("CRVAL3", crval3, "[deg] Coordinate value at reference point");
    hdu->card("CDELT3", m_phibin, "[deg] Coordinate increment at reference point");

    // Write OGIP time keywords
    m_gti.reference().write(*hdu);
    hdu->card("TSTART", m_gti.tstart().secs(), "[s] Start time");
    hdu->card("TSTOP",  m_gti.tstop().secs(), "[s] Stop time");
    hdu->card("DATE-OBS", m_gti.tstart().utc(), "Start of observation in UTC");
    hdu->card("DATE-END", m_gti.tstop().utc(), "Stop of observation in UTC");

    // Set time keywords
    int visday = com_tjd(m_gti.tstart());
    int vistim = com_tics(m_gti.tstart());
    int vieday = com_tjd(m_gti.tstop());
    int vietim = com_tics(m_gti.tstop());

    // Write COMPTEL time keywords
    hdu->card("VISDAY", visday, "[TJD] Data validity interval start day");
    hdu->card("VISTIM", vistim, "[tics] Data validity interval start time");
    hdu->card("VIEDAY", vieday, "[TJD] Data validity interval end day");
    hdu->card("VIETIM", vietim, "[tics] Data validity interval start time");

    // If there are energy boundaries then write them
    if (m_ebounds.size() > 0) {
        hdu->card("E_MIN", m_ebounds.emin().MeV(), "[MeV] Lower bound of energy range");
        hdu->card("E_MAX", m_ebounds.emax().MeV(), "[MeV] Upper bound of energy range");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute DRG geometry factor
 *
 * @param[in] tjd TJD for module status
 * @param[in] theta Zenith angle in COMPTEL coordinates (deg).
 * @param[in] phi Azimuth angle in COMPTEL coordinates (deg).
 * @param[in] status D1 and D2 module status
 * @return Geometry factor.
 *
 * Computes the DRG geometry factor as function of zenith and azimuth angles
 * given in the COMPTEL coordinate system.
 ***************************************************************************/
double GCOMDri::compute_geometry(const int&        tjd,
                                 const double&     theta,
                                 const double&     phi,
                                 const GCOMStatus& status) const
{
    // Set D1 module positions (from COM-RP-MPE-M10-123, Issue 1, Page 3-3)
    const double xd1[] = {0.0,-42.3,-26.0,26.0, 42.3,26.0,-26.0};
    const double yd1[] = {0.0,  0.0, 39.1,39.1, 0.0,-39.1,-39.1};

    // Set D2 module positions (from COM-RP-MPE-M10-123, Issue 1, Page 4-4)
    const double xd2[] = { 30.2,    0.0,  -30.2,   45.3,  15.1, -15.1, -45.3,
                           45.3,   15.1,  -15.1,  -45.3,  30.2,   0.0, -30.2};
    const double yd2[] = {-41.254,-41.254,-41.254,-15.1, -15.1, -15.1, -15.1,
                           15.1,   15.1,   15.1 ,  15.1,  41.254,41.254,41.254};

    // Set distance between D1 and D2 levels in cm (from COM-SP-MPE-M10-123,
    // Issue 1, page 8-5
    const double delz = 158.0;

    // Set D1 module radius in cm (from COM-TN-UNH-F70-051)
    const double r1 = 13.8;

    // Set D2 module radius in cm (from SIM-AL-005)
    const double r2 = 14.085;

    // Derive some constants
    const double r1sq     = r1 * r1;
    const double r2sq     = r2 * r2;
    const double drsq     = r1sq - r2sq;
    const double norm_geo = 1.0 / (gammalib::pi * r1sq);

    // Initialise geometry factor
    double geometry = 0.0;

    // Precompute results
    double theta_rad = theta * gammalib::deg2rad;
    double phi_rad   = phi   * gammalib::deg2rad;
    double cosphi    = std::cos(phi_rad);
    double sinphi    = std::sin(phi_rad);
    double tantheta  = std::tan(theta_rad);

    // Loop over all D1 modules
    for (int id1 = 0; id1 < 7; ++id1) {
    
        // Skip D1 module if it's off
        if (status.d1status(tjd, id1+1) != 1) {
            continue;
        }

        // Compute coordinates of D1 projected on D2 plane
        double xd1_prj = xd1[id1] - delz * tantheta * cosphi;
        double yd1_prj = yd1[id1] - delz * tantheta * sinphi;

        // Loop over all D2 modules
        for (int id2 = 0; id2 < 14; ++id2) {

            // Skip D2 module if it's off
            if (status.d2status(tjd, id2+1) != 1) {
                continue;
            }

            // Compute distance between center of D2 and projected centre
            // of D1
            double dx  = xd2[id2] - xd1_prj;
            double dy  = yd2[id2] - yd1_prj;
            double dsq = dx * dx + dy * dy;
            double d   = std::sqrt(dsq);

            // If there is no overlap then skip the module
            if (d >= r1 + r2) {
                continue;
            }

            // If there is total overlap (within 0.1 cm) then add 1 to the
            // geometry factor
            else if (d <= r2 - r1 + 0.1) {
                geometry += 1.0;
            }

            // ... otherwise if there is a partial overlap then compute the
            // semiangle subtended by overlap sector of D2 and projected
            // D1 module
            else {

                // Cosine beta
                double d2     = 2.0 * d;
                double cbeta1 = (dsq + drsq) / (d2 * r1);
                double cbeta2 = (dsq - drsq) / (d2 * r2);

                // Sin beta
                double sbeta1 = std::sqrt(1.0 - cbeta1 * cbeta1);
                double sbeta2 = std::sqrt(1.0 - cbeta2 * cbeta2);

                // Beta
                double beta1 = std::acos(cbeta1);
                double beta2 = std::acos(cbeta2);

                // Projection
                geometry += (r1sq * (beta1 - sbeta1 * cbeta1) +
                             r2sq * (beta2 - sbeta2 * cbeta2)) * norm_geo;

            } // endelse: there was partial overlap

        } // endfor: looped over D2 modules

    } // endfor: looped over D1 modules

    // Now divide by 7 since we want the probability of hitting a given
    // D1 module
    geometry /= 7.0;

    // Return geometry factor
    return geometry;
}