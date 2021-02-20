/***************************************************************************
 *                  GCOMDri.cpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2019 by Juergen Knoedlseder                         *
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
#include "GModelSky.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GCOMDri.hpp"
#include "GCOMOad.hpp"
#include "GCOMOads.hpp"
#include "GCOMTim.hpp"
#include "GCOMStatus.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventList.hpp"
#include "GCOMSelection.hpp"
#include "GCOMSupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COMPUTE_DRE               "GCOMDri::compute_dre(GCOMObservation&, "\
                                                   "GCOMSelection&, double&)"
#define G_COMPUTE_DRM       "GCOMDri::compute_drm(GCOMObservation&, GModel&)"
#define G_COMPUTE_DRE_PTSRC   "GCOMDri::compute_drm_ptsrc(GCOMObservation&, "\
                                                                "GModelSky&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_CHECK_EHA_COMPUTATION

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_COMPUTE_DRE
//#define G_DEBUG_COMPUTE_DRG
//#define G_DEBUG_COMPUTE_DRX

/* __ Constants __________________________________________________________ */
const double superpacket_duration = 16.384; // Duration of superpacket (s)
const double r1                   = 13.8;   // D1 module radius (cm)

/* __ Derived constants __________________________________________________ */
const double r1sq    = r1 * r1;
const double d1_area = 7.0 * gammalib::pi * r1sq;


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
 * @brief Sky map constructor
 *
 * @param[in] map Sky map defining the DRI cube.
 * @param[in] phimin Minimum Phibar angle (deg).
 * @param[in] phibin Bin size of Phibar angle (deg).
 * @param[in] nphibin Number of Phibar bins.
 *
 * Constructs a DRI cube from a sky map and a Phibar binning definition.
 ***************************************************************************/
GCOMDri::GCOMDri(const GSkyMap& map,
                 const double&  phimin,
                 const double&  phibin,
                 const int&     nphibin)
{
    // Initialise class members
    init_members();

    // Set sky map
    m_dri = map;

    // Make sure that sky map can hold the required number of maps
    if (nphibin > 0) {
        m_dri.nmaps(nphibin);
    }

    // Set all sky map pixels to zero
    m_dri = 0;

    // Store minimum Phibar and bin size
    m_phimin = phimin;
    m_phibin = phibin;

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
 * @param[in] obs COMPTEL observation.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 *
 * @exception GException::invalid_argument
 *            DRE cube has a non-positive Phibar bin size.
 *
 * Compute DRE event cube for a COMPTEL observation.
 ***************************************************************************/
void GCOMDri::compute_dre(const GCOMObservation& obs,
                          const GCOMSelection&   select,
                          const double&          zetamin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRE)
    std::cout << "GCOMDri::compute_dre" << std::endl;
    std::cout << "====================" << std::endl;
    #endif

    // Check if observation contains a COMPTEL event list
    const GCOMEventList* events = dynamic_cast<const GCOMEventList*>(obs.events());
    if (events == NULL) {
        std::string msg = "No event list found in COMPTEL observation. Please "
                          "specify an observation that contains an event list.";
        throw GException::invalid_argument(G_COMPUTE_DRE, msg);
    }

    // Check for positive Phibar bin size
    if (m_phibin <= 0.0) {
        std::string msg = "DRE cube has a non-positive Phibar bin size. Please "
                          "specify a DRE cube with a positive Phibar bin size.";
        throw GException::invalid_argument(G_COMPUTE_DRE, msg);
    }

    // Initialise event counter
    int i_evt = 0;

    // Initialise D1 & D2 module status
    GCOMStatus status;

    // Initialise superpacket statistics
    init_statistics();

    // Initialise selection statistics
    select.init_statistics();

    // Signal that event selection was used (this will enable writing into
    // the FITS HDU) and store Earth horizon cut
    m_has_selection = true;
    m_zetamin       = zetamin;

    // Initialise statistics
    int num_used_events           = 0;
    int num_event_outside_sp      = 0;
    int num_energy_too_low        = 0;
    int num_energy_too_high       = 0;
    int num_eha_too_small         = 0;
    int num_phibar_too_low        = 0;
    int num_phibar_too_high       = 0;
    int num_outside_dre           = 0;
    int num_event_before_dre      = 0;
    int num_event_after_dre       = 0;
    int num_d1module_off          = 0;
    int num_d2module_off          = 0;
    int num_processed             = 0;

    // Set all DRE bins to zero
    init_cube();

    // Signal that loop should be terminated
    bool terminate = false;

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < obs.oads().size(); ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = obs.oads()[i_oad];

        // Check superpacket usage
        if (!use_superpacket(oad, obs.tim())) {
            continue;
        }

        // Prepare Earth horizon angle comparison
        GSkyDir sky_geocentre;
        double  theta_geocentre = double(oad.gcel());
        double  phi_geocentre   = double(oad.gcaz());
        sky_geocentre.radec_deg(phi_geocentre, 90.0-theta_geocentre);

        // Collect all events within superpacket. Break if the end
        // of the event list was reached.
        for (; i_evt < events->size(); ++i_evt) {

            // Get pointer to event
            const GCOMEventAtom* event = (*events)[i_evt];

            // Break loop if the end of the superpacket was reached
            if (event->time() > oad.tstop()) {
                break;
            }

            // Increase number of processed events
            num_processed++;

            // Check GTIs if the DRE has GTIs
            if (m_gti.size() > 0) {

                // Skip event if it lies before the DRE start.
                if (event->time() < m_gti.tstart()) {
                    num_event_before_dre++;
                    continue;
                }

                // Break if event lies after the DRE stop
                else if (event->time() > m_gti.tstop()) {
                    num_event_after_dre++;
                    terminate = true;
                    break;
                }
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

            // Apply event selection
            if (!select.use_event(*event)) {
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

            // Option: Earth horizon angle comparison
            #if defined(G_CHECK_EHA_COMPUTATION)
            GSkyDir sky_event;
            double  theta_event = double(event->theta());
            double  phi_event   = double(event->phi());
            sky_event.radec_deg(-phi_event, 90.0-theta_event);
            double eha = sky_geocentre.dist_deg(sky_event) - oad.georad();
            if (std::abs(eha - event->eha()) > 1.5) {
                std::string msg = "Earth horizon angle from EVP dataset ("+
                                  gammalib::str(event->eha())+" deg) "
                                  "differs from Earth horizon angle "
                                  "computed from Orbit Aspect Data ("+
                                  gammalib::str(eha)+" deg). Use the EVP "
                                  "value.";
                gammalib::warning(G_COMPUTE_DRE, msg);
            }
            #endif

            // Check for Earth horizon angle. There is a constant EHA limit
            // over a Phibar layer to be compliant with the DRG.
            double ehamin = double(iphibar) * m_phibin + zetamin;
            if (event->eha() < ehamin) {
                num_eha_too_small++;
                continue;
            }

            // Now fill the event into the DRE. Put this in a try-catch
            // block so that any invalid transformations are catched.
            try {
                GSkyPixel pixel = m_dri.dir2pix(event->dir().dir());
                if (m_dri.contains(pixel)) {
                    int inx        = m_dri.pix2inx(pixel) + iphibar * m_dri.npix();
                    (*this)[inx] += 1.0;
                    num_used_events++;
                }
                else {
                    num_outside_dre++;
                }
            }
            catch (GException::wcs_invalid_phi_theta) {
                num_outside_dre++;
            }

        } // endfor: collected events

        // Break if termination was signalled or if there are no more events
        if (terminate || i_evt >= events->size()) {
            break;
        }

    } // endfor: looped over Orbit Aspect Data

    // Set Good Time interval for DRE
    if (m_num_used_superpackets > 0) {
        m_gti = GGti(m_tstart, m_tstop);
    }

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRE)
    std::cout << "Total number of superpackets .: " << m_num_superpackets << std::endl;
    std::cout << "Used superpackets ............: " << m_num_used_superpackets << std::endl;
    std::cout << "Skipped superpackets .........: " << m_num_skipped_superpackets << std::endl;
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
    std::cout << "Outside DRE cube .............: " << num_outside_dre << std::endl;
    std::cout << "Earth horizon angle too small : " << num_eha_too_small << std::endl;
    std::cout << select << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute geometry cube
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] select Selection set (not used so far).
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 *
 * Compute DRG cube for a COMPTEL observation.
 ***************************************************************************/
void GCOMDri::compute_drg(const GCOMObservation& obs,
                          const GCOMSelection&   select,
                          const double&          zetamin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRG)
    std::cout << "GCOMDri::compute_drg" << std::endl;
    std::cout << "====================" << std::endl;
    #endif

    // Initialise D1 & D2 module status
    GCOMStatus status;

    // Initialise superpacket statistics
    init_statistics();

    // Initialise selection statistics
    select.init_statistics();

    // Set all DRG bins to zero
    init_cube();

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < obs.oads().size(); ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = obs.oads()[i_oad];

        // Check superpacket usage
        if (!use_superpacket(oad, obs.tim())) {
            continue;
        }

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

        // Compute geometrical factors for all (Chi, Psi)
        for (int index = 0; index < m_dri.npix(); ++index) {

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
            chipsi_comptel.radec_deg(phi, 90.0-theta);
            double eha = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

            // Loop over Phibar
            for (int iphibar = 0; iphibar < nphibar(); ++iphibar) {

                // Compute minimum Earth horizon angle for Phibar layer
                double ehamin = double(iphibar) * m_phibin + zetamin;

                // Add up geometry if the Earth horizon angle is equal or
                // larger than the minimum
                if (eha >= ehamin) {
                    int inx       = index + iphibar * m_dri.npix();
                    (*this)[inx] += geometry;
                }

            } // endfor: looped over Phibar

        } // endfor: looped over (Chi, Psi)

    } // endfor: looped over Orbit Aspect Data

    // Divide DRG by number of used superpackets
    if (m_num_used_superpackets > 0) {
        double norm = 1.0 / double(m_num_used_superpackets);
        for (int i = 0; i < size(); ++i) {
            (*this)[i] *= norm;
        }
    }

    // Clear energy boundaries for DRG since DRG does not depend on energy)
    m_ebounds.clear();

    // Set the Good Time interval for DRG
    if (m_num_used_superpackets > 0) {
        m_gti = GGti(m_tstart, m_tstop);
    }

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRG)
    std::cout << "Total number of superpackets .: " << m_num_superpackets << std::endl;
    std::cout << "Used superpackets ............: " << m_num_used_superpackets << std::endl;
    std::cout << "Skipped superpackets .........: " << m_num_skipped_superpackets << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute DRX exposure map
 *
 * @param[in] obs COMPTEL observation.
 *
 * Compute DRX exposure map for a COMPTEL observation.
 *
 * For a given superpacket, the exposure is computed using
 *
 * \f[
 *    X_i(\theta_c) = 7 \pi r_1^2 \cos \theta_c
 *    \frac{1 - \exp \left( -\tau \ \cos \theta_c \right)}
 *         {1 - \exp \left( -\tau \right)}
 * \f]
 *
 * where
 * \f$\tau=0.2\f$ is the typical thickness of a D1 module in radiation
 * lengths,
 * \f$r_1=13.8\f$ cm is the radius of a D1 module, and
 * \f$\theta_c\f$ is the zenith angle in COMPTEL coordinates.
 ***************************************************************************/
void GCOMDri::compute_drx(const GCOMObservation& obs)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRX)
    std::cout << "GCOMDri::compute_drx" << std::endl;
    std::cout << "====================" << std::endl;
    #endif

    // Initialise constants
    const double tau = 0.2; // Thickness in radiation lenghts

    // Initialise superpacket statistics
    init_statistics();

    // Set all DRX bins to zero
    init_cube();

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < obs.oads().size(); ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = obs.oads()[i_oad];

        // Check superpacket usage
        if (!use_superpacket(oad, obs.tim())) {
            continue;
        }

        // Loop over all DRX pixels
        for (int index = 0; index < m_dri.npix(); ++index) {

            // Get sky direction for DRX pixel
            GSkyDir sky = m_dri.inx2dir(index);

            // Compute zenith angle of pixel in COMPTEL coordinates
            double theta = oad.theta(sky);

            // Initialise exposure
            double exposure = 0.0;

            // Skip pixel if zenith angle is beyond 90 degrees
            if (theta >= 90.0) {
                continue;
            }

            // ... otherwise compute the exposure
            else {
                double costheta = std::cos(theta * gammalib::deg2rad);
                if (theta < 89.0) {
                    exposure = d1_area * costheta *
                               (1.0 - std::exp(-tau / costheta)) /
                               (1.0 - std::exp(-tau));
                }
                else {
                    exposure = d1_area * costheta / (1.0 - std::exp(-tau));
                }
            }

            // Accumulate exposure
            (*this)[index] += exposure;

        } // endfor: looped over all DRX pixels

    } // endfor: looped over Orbit Aspect Data

    // Multiply by time per superpacket to give the result in cm^2 s
    for (int i = 0; i < size(); ++i) {
        (*this)[i] *= superpacket_duration;
    }

    // Clear energy boundaries for DRX since DRX does not depend on energy)
    m_ebounds.clear();

    // Set the Good Time interval for DRX
    if (m_num_used_superpackets > 0) {
        m_gti = GGti(m_tstart, m_tstop);
    }

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRX)
    std::cout << "Total number of superpackets .: " << m_num_superpackets << std::endl;
    std::cout << "Used superpackets ............: " << m_num_used_superpackets << std::endl;
    std::cout << "Skipped superpackets .........: " << m_num_skipped_superpackets << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute DRM model
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] model Model.
 *
 * @exception GException::feature_not_implemented
 *            Method is only implemented for a point source sky model.
 *
 * Compute DRM model cube for a COMPTEL observation.
 ***************************************************************************/
void GCOMDri::compute_drm(const GCOMObservation& obs,
                          const GModel&          model)
{
    // Check if model is a sky model
    const GModelSky* skymodel = dynamic_cast<const GModelSky*>(&model);
    if (skymodel == NULL) {
        std::string msg = "Method is only implement for sky models.";
        throw GException::feature_not_implemented(G_COMPUTE_DRM, msg);
    }

    // Extract spatial model component
    const GModelSpatial* spatial = skymodel->spatial();

    // Select computation depending on the spatial model type
    switch (spatial->code()) {
        case GMODEL_SPATIAL_POINT_SOURCE:
            {
            compute_drm_ptsrc(obs, *skymodel);
            }
            break;
        case GMODEL_SPATIAL_RADIAL:
        case GMODEL_SPATIAL_ELLIPTICAL:
        case GMODEL_SPATIAL_DIFFUSE:
            {
            std::string msg = "Method is not yet implemented for spatial model "
                              "type \""+spatial->type()+"\".";
            throw GException::feature_not_implemented(G_COMPUTE_DRM, msg);
            }
            break;
        default:
            break;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute content in cone
 *
 * @param[in] dir Sky direction of cone apex.
 * @param[in] armmin Minimum Angular Resolution Measure (deg).
 * @param[in] armmax Maximum Angular Resolution Measure (deg).
 * @return Content in cone.
 *
 * Compute the sum of the DRI bins within an event cone with apex at a given
 * sky direction. All bins with an Angular Resolution Measure comprised
 * between @p armmin (inclusive) and @p armmax (exclusive) will be
 * considered. The bin centres will be used for the computation of the
 * Angular Resolution Measure. The Angular Resolution Measure is defined as
 * phibar - phigeo.
 ***************************************************************************/
double GCOMDri::cone_content(const GSkyDir& dir,
                             const double&  armmin,
                             const double&  armmax) const
{
    // Initialise content
    double content = 0.0;

    // Create Phigeo map in degrees
    GSkyMap phigeo = m_dri.extract(0);
    for (int i = 0; i < phigeo.npix(); ++i) {
        phigeo(i) = dir.dist_deg(phigeo.inx2dir(i));
    }

    // Loop over phibar layers
    for (int iphibar = 0; iphibar < this->nphibar(); ++iphibar) {

        // Compute phibar
        double phibar = this->phimin() + (iphibar+0.5) * this->phibin();

        // Compute index offset
        int offset = iphibar * phigeo.npix();

        // Loop over bins in phibar layer
        for (int i = 0; i < phigeo.npix(); ++i) {

            // Compute ARM
            double arm = phibar - phigeo(i);

            // If ARM is within specified interval then add bin content
            if (arm >= armmin && arm < armmax) {
                content += (*this)[i + offset];
            }

        } // endfor: looped over bins in phibar layers

    } // endfor: looped over phibar layers

    // Return content
    return content;
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
 * @param[in] clobber Overwrite existing file?
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
    gammalib::com_wcs_mer2car(m_dri);

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
    if (image != NULL) {
        write_attributes(image);
    }

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
        if (m_tofcor > 1.0) {
            result.append("\n"+gammalib::parformat("ToF correction"));
            result.append(gammalib::str(m_tofcor));
        }

        // Append energy boundaries
        result.append("\n"+m_ebounds.print(gammalib::reduce(chatter)));

        // Append GTI
        result.append("\n"+m_gti.print(gammalib::reduce(chatter)));

        // EXPLICIT: Append sky map
        if (chatter >= EXPLICIT) {
            result.append("\n"+m_dri.print(gammalib::reduce(chatter)));
        }

        // Append computation statistics
        result.append("\n"+gammalib::parformat("Input superpackets"));
        result.append(gammalib::str(m_num_superpackets));
        result.append("\n"+gammalib::parformat("Used superpackets"));
        result.append(gammalib::str(m_num_used_superpackets));
        result.append("\n"+gammalib::parformat("Skipped superpackets"));
        result.append(gammalib::str(m_num_skipped_superpackets));

        // Append selection
        if (m_has_selection) {
            result.append("\n"+m_selection.print(gammalib::reduce(chatter)));
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
    m_name.clear();
    m_dri.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_phimin = 0.0;
    m_phibin = 0.0;
    m_tofcor = 1.0;

    // Initialise statistics
    init_statistics();

    // Initialise selection parameters
    m_has_selection = false;
    m_selection.clear();
    m_zetamin = -90.0;   // Signals no selection

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
    m_name    = dri.m_name;
    m_dri     = dri.m_dri;
    m_ebounds = dri.m_ebounds;
    m_gti     = dri.m_gti;
    m_phimin  = dri.m_phimin;
    m_phibin  = dri.m_phibin;
    m_tofcor  = dri.m_tofcor;

    // Copy statistics
    m_tstart                   = dri.m_tstart;
    m_tstop                    = dri.m_tstop;
    m_num_superpackets         = dri.m_num_superpackets;
    m_num_used_superpackets    = dri.m_num_used_superpackets;
    m_num_skipped_superpackets = dri.m_num_skipped_superpackets;

    // Copy selection parameters
    m_has_selection = dri.m_has_selection;
    m_selection     = dri.m_selection;
    m_zetamin       = dri.m_zetamin;

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
 * @brief Initialise DRI cube
 *
 * Sets all DRI cube bins to zero.
 ***************************************************************************/
void GCOMDri::init_cube(void)
{
    // Set all cube bins to zero
    for (int i = 0; i < size(); ++i) {
        (*this)[i] = 0.0;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise computation statistics
 ***************************************************************************/
void GCOMDri::init_statistics(void)
{
    // Initialise statistics
    m_tstart.clear();
    m_tstop.clear();
    m_num_superpackets         = 0;
    m_num_used_superpackets    = 0;
    m_num_skipped_superpackets = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if superpacket should be used
 *
 * @param[in] oad Orbit Aspect Data record (i.e. superpacket).
 * @param[in] tim Good Time Intervals.
 * @return True if superpacket should be used, false otherwise.
 *
 * Checks if a superpacket should be used and updated the superpacket
 * statistics and selected time interval. A superpacket will be used if
 * it is fully enclosed within the COMPTEL Good Time Intervals and the
 * Good Time Intervals of the DRI dataset.
 ***************************************************************************/
bool GCOMDri::use_superpacket(const GCOMOad &oad, const GCOMTim& tim)
{
    // Initialise usage flag
    bool use = true;

    // Increment superpacket counter
    m_num_superpackets++;

    // Skip superpacket if it is not fully enclosed within the COMPTEL
    // Good Time Intervals
    if (!(tim.contains(oad.tstart()) && tim.contains(oad.tstop()))) {
        m_num_skipped_superpackets++;
        use = false;
    }

    // Skip superpacket if it is not fully enclosed within the DRI Good
    // Time Intervals. Only check if there are Good Time Intervals in the
    // DRI.
    else if ((m_gti.size() > 0) &&
             !(m_gti.contains(oad.tstart()) && m_gti.contains(oad.tstop()))) {
        m_num_skipped_superpackets++;
        use = false;
    }

    // ... otherwise use superpacket
    else {
        m_num_used_superpackets++;
    }

    // Update selection validity interval
    if (m_num_used_superpackets == 1) {
        m_tstart = oad.tstart();
        m_tstop  = oad.tstop();
    }
    else {
        m_tstop  = oad.tstop();
    }

    // Return usage
    return use;
}


/***********************************************************************//**
 * @brief Read DRI attributes from FITS HDU
 *
 * @param[in] hdu FITS HDU pointer.
 *
 * Reads the time interval from the FITS header and sets the Phibar definiton
 * and energy boundaries from the header keywords if they are provided.
 ***************************************************************************/
void GCOMDri::read_attributes(const GFitsHDU* hdu)
{
    // Get time attributes
    GTime tstart = gammalib::com_time(hdu->integer("VISDAY"), hdu->integer("VISTIM"));
    GTime tstop  = gammalib::com_time(hdu->integer("VIEDAY"), hdu->integer("VIETIM"));

    // Set Good Time Intervals
    if (tstop > tstart) {
        m_gti = GGti(tstart, tstop);
    }

    // Optionally read Phibar attributes
    if (hdu->has_card("CDELT3") &&
        hdu->has_card("CRVAL3") &&
        hdu->has_card("CRPIX3")) {

        // Get phibar attributes
        m_phibin = hdu->real("CDELT3");
        m_phimin = hdu->real("CRVAL3") - (hdu->real("CRPIX3")-0.5) * m_phibin;

    }

    // ... otherwise set Phibar attributes to zero
    else {
        m_phimin = 0.0;
        m_phibin = 0.0;
    }

    // Optionally read energy attributes
    if (hdu->has_card("E_MIN") && hdu->has_card("E_MAX")) {

        // Get energy attributes
        GEnergy emin(hdu->real("E_MIN"), "MeV");
        GEnergy emax(hdu->real("E_MAX"), "MeV");

        // Set energy boundaries
        m_ebounds = GEbounds(emin, emax);

    }

    // ... otherwise clear energy boundaries
    else {
        m_ebounds.clear();
    }

    // Read selection set
    m_selection.read(*hdu);

    // Optionally read Earth horizon cut
    if (hdu->has_card("EHAMIN")) {
        m_zetamin = hdu->real("EHAMIN");
    }

    // Optionally read ToF correction
    if (hdu->has_card("TOFCOR")) {
        m_tofcor = hdu->real("TOFCOR");
    }

    // Optionally read superpacket statistics
    if (hdu->has_card("NSPINP")) {
        m_num_superpackets = hdu->integer("NSPINP");
    }
    if (hdu->has_card("NSPUSE")) {
        m_num_used_superpackets = hdu->integer("NSPUSE");
    }
    if (hdu->has_card("NSPSKP")) {
        m_num_skipped_superpackets = hdu->integer("NSPSKP");
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
    // Write 3rd dimension keywords in case that a third dimension exists.
    // This avoids writing a Phibar layer for DRX.
    if (m_phibin > 0.0) {

        // Set Phibar keywords
        double crval3 = m_phimin + 0.5 * m_phibin;

        // Write Phibar keywords
        hdu->card("CTYPE3", "Phibar", "Compton scatter angle");
        hdu->card("CRPIX3", 1.0, "Pixel coordinate of reference point (starting from 1)");
        hdu->card("CRVAL3", crval3, "[deg] Coordinate value at reference point");
        hdu->card("CDELT3", m_phibin, "[deg] Coordinate increment at reference point");

    } // endif: DRI was 3D

    // Write OGIP time keywords
    m_gti.reference().write(*hdu);
    hdu->card("TSTART", m_gti.tstart().secs(), "[s] Start time");
    hdu->card("TSTOP",  m_gti.tstop().secs(), "[s] Stop time");
    hdu->card("DATE-OBS", m_gti.tstart().utc(), "Start of observation in UTC");
    hdu->card("DATE-END", m_gti.tstop().utc(), "Stop of observation in UTC");

    // Set time keywords
    int visday = gammalib::com_tjd(m_gti.tstart());
    int vistim = gammalib::com_tics(m_gti.tstart());
    int vieday = gammalib::com_tjd(m_gti.tstop());
    int vietim = gammalib::com_tics(m_gti.tstop());

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

    // If there was an event seletion then write selection set
    if (m_has_selection) {
        m_selection.write(*hdu);
    }

    // If there was a valid Earth horizon cut then write it
    if (m_zetamin != -90.0) {
        hdu->card("EHAMIN", m_zetamin, "[deg] Minimum Earth horizon - Phibar cut");
    }

    // If there is a ToF correction then write it
    if (m_tofcor > 1.0) {
        hdu->card("TOFCOR", m_tofcor, "ToF correction");
    }

    // Write superpacket statistics
    hdu->card("NSPINP", m_num_superpackets, "Number of input superpackets");
    hdu->card("NSPUSE", m_num_used_superpackets, "Number of used superpackets");
    hdu->card("NSPSKP", m_num_skipped_superpackets, "Number of skipped superpackets");

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
    //const double r1 = 13.8; // Defined globally

    // Set D2 module radius in cm (from SIM-AL-005)
    const double r2 = 14.085;

    // Derive some constants
    //const double r1sq     = r1 * r1;  // Defined globally
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


/***********************************************************************//**
 * @brief Compute DRM model for a point source
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] model Sky model.
 *
 * @exception GException::invalid_argument
 *            Model is not a point source model.
 *
 * Compute point source DRM model cube for a COMPTEL observation.
 ***************************************************************************/
void GCOMDri::compute_drm_ptsrc(const GCOMObservation& obs,
                                const GModelSky&       model)
{
    // Extract source direction from spatial model component
    const GModelSpatialPointSource* source =
        dynamic_cast<const GModelSpatialPointSource*>(model.spatial());
    if (source == NULL) {
        std::string msg = "Spatial component of model \""+model.name()+"\" "
                          "is not a point source.";
        throw GException::invalid_argument(G_COMPUTE_DRE_PTSRC, msg);
    }
    const GSkyDir& srcDir = source->dir();

    // Extract response
    const GCOMResponse* rsp = obs.response();

    // Set all DRM cube bins to zero
    init_cube();

    // Get DRX value (units: cm^2 sec)
    double drx = obs.drx()(srcDir);

    // Get ontime
    double ontime = obs.ontime(); // sec

    // Get deadtime correction
    double deadc = obs.deadc(GTime());

    // Compute multiplicate IAQ normalisation
    double norm = drx * deadc / ontime;

    // Loop over all scatter angle pixels
    for (int index = 0; index < m_dri.npix(); ++index) {

        // Get sky direction of pixel
        GSkyDir obsDir = m_dri.inx2dir(index);

        // Compute angle between true photon arrival direction and scatter
        // direction (Chi,Psi)
        double phigeo = srcDir.dist_deg(obsDir);

        // Get Phigeo interpolation factor
        double phirat  = phigeo / rsp->m_phigeo_bin_size; // 0.5 at bin centre
        int    iphigeo = int(phirat);                     // index into which Phigeo falls
        double eps     = phirat - iphigeo - 0.5;          // 0.0 at bin centre

        // Continue only if Phigeo is inside IAQ
        if (iphigeo < rsp->m_phigeo_bins) {

            // Initialise IAQ pixel index
            int i = iphigeo;

            // Interpolate towards left
            if (eps < 0.0) {

                // Not the first bin
                if (iphigeo > 0) {
                    for (int iphibar = 0; iphibar < nphibar(); ++iphibar,
                         i += rsp->m_phigeo_bins) {
                        double iaq = obs.drg()(index, iphibar);
                        iaq       *= (1.0 + eps) * rsp->m_iaq[i] - eps * rsp->m_iaq[i-1];
                        m_dri(index, iphibar) = iaq * norm;
                    }
                }
                else {
                    for (int iphibar = 0; iphibar < nphibar(); ++iphibar,
                         i += rsp->m_phigeo_bins) {
                        double iaq = obs.drg()(index, iphibar);
                        iaq       *= (1.0 - eps) * rsp->m_iaq[i] + eps * rsp->m_iaq[i+1];
                        m_dri(index, iphibar) = iaq * norm;
                    }
                }

            }

            // Interpolate towards right
            else {

                // Not the last IAQ bin
                if (iphigeo < rsp->m_phigeo_bins-1) {
                    for (int iphibar = 0; iphibar < nphibar(); ++iphibar,
                         i += rsp->m_phigeo_bins) {
                        double iaq = obs.drg()(index, iphibar);
                        iaq       *= (1.0 - eps) * rsp->m_iaq[i] + eps * rsp->m_iaq[i+1];
                        m_dri(index, iphibar) = iaq * norm;
                    }
                }
                else {
                    for (int iphibar = 0; iphibar < nphibar(); ++iphibar,
                         i += rsp->m_phigeo_bins) {
                        double iaq = obs.drg()(index, iphibar);
                        iaq       *= (1.0 + eps) * rsp->m_iaq[i] - eps * rsp->m_iaq[i-1];
                        m_dri(index, iphibar) = iaq * norm;
                    }
                }

            } // endfor: looped over all Compton scatter angles

        } // endif: Phigeo was inside IAQ

    } // endfor: looped over all scatter angle pixels

    // Return
    return;
}
