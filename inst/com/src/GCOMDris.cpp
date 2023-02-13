/***************************************************************************
 *             GCOMDris.cpp - COMPTEL Data Space container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMDris.hpp
 * @brief COMPTEL Data Space container class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMDri.hpp"
#include "GCOMDris.hpp"
#include "GCOMOad.hpp"
#include "GCOMOads.hpp"
#include "GCOMTim.hpp"
#include "GCOMStatus.hpp"
#include "GCOMObservation.hpp"
#include "GCOMSelection.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMDris::at(int&)"
#define G_INSERT                           "GCOMDris::insert(int&, GCOMDri&)"
#define G_REMOVE                                     "GCOMDris::remove(int&)"
#define G_COMPUTE_DRWS            "GCOMDris::compute_drws(GCOMObservation&, "\
                                                   "GCOMSelection&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_COMPUTE_DRWS


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty Data Space container.
 ***************************************************************************/
GCOMDris::GCOMDris(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dris Data Space container.
 ***************************************************************************/
GCOMDris::GCOMDris(const GCOMDris& dris)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dris);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMDris::~GCOMDris(void)
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
 * @param[in] dris Data Space container.
 * @return Data Space container.
 ***************************************************************************/
GCOMDris& GCOMDris::operator=(const GCOMDris& dris)
{
    // Execute only if object is not identical
    if (this != &dris) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dris);

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
 * @brief Clear Data Space container
 ***************************************************************************/
void GCOMDris::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Data Space container
 *
 * @return Pointer to deep copy of Data Space container.
 ***************************************************************************/
GCOMDris* GCOMDris::clone(void) const
{
    return new GCOMDris(*this);
}


/***********************************************************************//**
 * @brief Return reference to Data Space
 *
 * @param[in] index Data Space index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Returns a reference to the Data Space with the specified @p index.
 ***************************************************************************/
GCOMDri& GCOMDris::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Data Space index",
                                       index, size());
    }

    // Return reference
    return m_dris[index];
}


/***********************************************************************//**
 * @brief Return reference to Data Space (const version)
 *
 * @param[in] index Data Space index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Returns a reference to the Data Space with the specified @p index.
 ***************************************************************************/
const GCOMDri& GCOMDris::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Data Space index",
                                       index, size());
    }

    // Return reference
    return m_dris[index];
}


/***********************************************************************//**
 * @brief Append Data Space to container
 *
 * @param[in] dri Data Space.
 * @return Reference to appended Data Space.
 *
 * Appends Data Space to the container by making a deep copy of the Data
 * Space.
 ***************************************************************************/
GCOMDri& GCOMDris::append(const GCOMDri& dri)
{
    // Append dri to list
    m_dris.push_back(dri);

    // Return reference
    return m_dris[size()-1];
}


/***********************************************************************//**
 * @brief Insert Data Space into container
 *
 * @param[in] index Data Space index (0,...,size()-1).
 * @param[in] dri Data Space.
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Inserts a Data Space into the container before the Data Space with the
 * specified @p index.
 ***************************************************************************/
GCOMDri& GCOMDris::insert(const int& index, const GCOMDri& dri)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Data Space index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Data Space index",
                                           index, size());
        }
    }
    #endif

    // Inserts Data Space
    m_dris.insert(m_dris.begin()+index, dri);

    // Return reference
    return m_dris[index];
}


/***********************************************************************//**
 * @brief Remove Data Space from container
 *
 * @param[in] index Data Space index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Remove Data Space of specified @p index from container.
 ***************************************************************************/
void GCOMDris::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Data Space index",
                                       index, size());
    }
    #endif

    // Erase Data Space from container
    m_dris.erase(m_dris.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Data Space container
 *
 * @param[in] oads Data Space container.
 *
 * Append Data Space container to the container.
 ***************************************************************************/
void GCOMDris::extend(const GCOMDris& dris)
{
    // Do nothing if Data Space container is empty
    if (!dris.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = dris.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_dris.push_back(dris[i]);
        }

    } // endif: Data Space container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute backgroud weighting cubes
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 *
 * @exception GException::invalid_value
 *            DRW with mismatching definition encountered in container.
 * @exception GException::invalid_argument
 *            No event list found in COMPTEL observation.
 *
 * Compute DRW cubes for a COMPTEL observation. DRW cubes are event-rate
 * weighted geometry cubes, where event rate are time and energy-dependent.
 *
 * Three hard-wired energy ranges are defined for which the event rate is
 * determined on a 5 minutes basis:
 *
 *       E < 1.8 MeV
 *       1.8 < E < 4.3 MeV
 *       E > 4.3 MeV
 *
 * For each superpacket, these event rates are linearly interpolated to the
 * relevant superpacket time.
 *
 * The DRW cubes are normalised so that the sum of each Phibar layer equals
 * to unity.
 ***************************************************************************/
void GCOMDris::compute_drws(const GCOMObservation& obs,
                            const GCOMSelection&   select,
                            const double&          zetamin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::compute_drws" << std::endl;
    std::cout << "======================" << std::endl;
    #endif

    // Set constants
    const double dtime = 5.0 * 60.0; // Event rate integration time in seconds

    // Continue only if there are DRWs in container
    if (!is_empty()) {

        // Debug
        #if defined(G_DEBUG_COMPUTE_DRWS)
        std::cout << "Initialisation" << std::endl;
        std::cout << "--------------" << std::endl;
        #endif

        // Initialise D1 & D2 module status
        GCOMStatus status;

        // Get pointer to first DRW in container (for convenient access to
        // members of GCOMDri; we won't use it for modifications, hence we
        // can declare it as constant)
        const GCOMDri* dri = &(m_dris[0]);

        // Initialise all DRWs
        for (int i = 0; i < size(); ++i) {

            // Throw an exception if the DRW definition differs
            if (dri->m_dri != m_dris[i].m_dri) {
                std::string msg = "Mismatch between definition of DRW "+
                                  gammalib::str(i)+" and the first DRW. "
                                  "Method only works on DRWs with identical "
                                  "definitions.";
                throw GException::invalid_value(G_COMPUTE_DRWS, msg);
            }

            // Initialise superpacket statistics
            m_dris[i].init_statistics();

            // Store selection set so that handling of failed PMT flag is
            // correctly ritten into the FITS header
            m_dris[i].m_selection = select;

            // Set all DRG bins to zero
            m_dris[i].init_cube();

            // Debug
            #if defined(G_DEBUG_COMPUTE_DRWS)
            std::cout << m_dris[i].print() << std::endl;
            #endif

        } // endfor: looped over all DRWs

        // Get pointer to event list. Throw an exception if the observation
        // did not contain an event list
        const GCOMEventList* events = dynamic_cast<const GCOMEventList*>(obs.events());
        if (events == NULL) {
            std::string msg = "No event list found in COMPTEL observation. Please "
                              "specify an observation that contains an event list.";
            throw GException::invalid_argument(G_COMPUTE_DRWS, msg);
        }

        // Debug
        #if defined(G_DEBUG_COMPUTE_DRWS)
        std::cout << "Determine energy-dependent event rates" << std::endl;
        std::cout << "--------------------------------------" << std::endl;
        #endif

        // Initialise current time and vectors for computation of node function
        // and vectors for energy-dependent rate interpolation
        double              time;
        double              rate0 = 0.0;
        double              rate1 = 0.0;
        double              rate2 = 0.0;
        GNodeArray          times;
        std::vector<double> rates0;
        std::vector<double> rates1;
        std::vector<double> rates2;

        // Set time of first event
        if (events->size() > 0) {
            time = (*events)[0]->time().secs();
        }

        // Determine energy-dependent event rates. The standard event selection
        // criteria and earth horizon cut is applied so that the event rate
        // corresponds to the current selection, yet we do not consider any
        // Good Time Intervals or OAD time information as we later only access
        // relevant times (YET THERE COULD BE AN ISSUE WITH A GLITCH). Rates
        // are computed for three hard-coded energy intervals.
        for (int i_evt = 0; i_evt < events->size(); ++i_evt) {

            // Get pointer to event
            const GCOMEventAtom* event = (*events)[i_evt];

            // Apply event selection
            if (!select.use_event(*event)) {
                continue;
            }

            // Apply Earth horizon cut
            if ((event->eha() - event->phibar()) < zetamin) {
                continue;
            }

            // If time exceeds time bin then add rates
            while (event->time().secs() > (time + dtime)) {

                // If rates are not all zero then append rates to vectors
                if ((rate0 != 0.0) || (rate1 != 0.0) || (rate2 != 0.0)) {

                    // Append time at mid-point of interval
                    times.append(time + 0.5 * dtime);

                    // Append rates
                    rates0.push_back(rate0);
                    rates1.push_back(rate1);
                    rates2.push_back(rate2);

                    // Debug
                    #if defined(G_DEBUG_COMPUTE_DRWS)
                    std::cout << "t=" << time + 0.5 * dtime;
                    std::cout << " r0=" << rate0;
                    std::cout << " r1=" << rate1;
                    std::cout << " r2=" << rate2;
                    std::cout << std::endl;
                    #endif

                    // Initialise rates for new accumulation
                    rate0 = 0.0;
                    rate1 = 0.0;
                    rate2 = 0.0;

                } // endif: there were rates to append

                // Increment time
                time += dtime;

            } // endwhile: event time was after the current integration time stop

            // Accumulate rates
            double energy = event->energy().MeV();
            if (energy < 1.8) {
                rate0 += 1.0;
            }
            else if (energy < 4.3) {
                rate1 += 1.0;
            }
            else {
                rate2 += 1.0;
            }

        } // endfor: looped over events

        // Debug
        #if defined(G_DEBUG_COMPUTE_DRWS)
        std::cout << "Compute DRWs" << std::endl;
        std::cout << "------------" << std::endl;
        #endif

        // Get Good Time Intervals. If a pulsar selection is specified then
        // reduce the Good Time Intervals to the validity intervals of the
        // pulsar emphemerides.
        GCOMTim tim = obs.tim();
        if (select.has_pulsar()) {
            tim.reduce(select.pulsar().validity());
        }

        // Loop over Orbit Aspect Data
        for (int i_oad = 0; i_oad < obs.oads().size(); ++i_oad) {

            // Get reference to Orbit Aspect Data of superpacket
            const GCOMOad &oad = obs.oads()[i_oad];

            // Check superpacket usage for all DRWs so that their statistics
            // get correctly updated
            bool skip = false;
            for (int i = 0; i < size(); ++i) {
                if (!m_dris[i].use_superpacket(oad, tim, select)) {
                    skip = true;
                }
            }
            if (skip) {
                continue;
            }

            // Get Orbit Aspect Data mid-time in seconds
            double time = oad.tstart().secs() + 0.5 * (oad.tstop() - oad.tstart());

            // Compute energy-dependent rates
            std::vector<double> rates;
            for (int i = 0; i < size(); ++i) {
                if (m_dris[i].m_ebounds.emin().MeV() > 4.3) {
                    rates.push_back(times.interpolate(time, rates2));
                }
                else if (m_dris[i].m_ebounds.emin().MeV() > 1.8) {
                    rates.push_back(times.interpolate(time, rates1));
                }
                else {
                    rates.push_back(times.interpolate(time, rates0));
                }
            }

            // Debug
            #if defined(G_DEBUG_COMPUTE_DRWS)
            std::cout << "time=" << time;
            for (int i = 0; i < size(); ++i) {
                std::cout << " r[" << i << "]=" << rates[i];
            }
            std::cout << std::endl;
            #endif

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

            // Compute geometrical factors for all (Chi, Psi). We assume here
            // that their size and definition is identical, which was checked
            // before
            for (int index = 0; index < dri->m_dri.npix(); ++index) {

                // Get sky direction and solid angle for (Chi, Psi)
                GSkyDir sky   = dri->m_dri.inx2dir(index);
                double  omega = dri->m_dri.solidangle(index);

                // Convert sky direction to COMPTEL coordinates
                double theta = oad.theta(sky);
                double phi   = oad.phi(sky);

                // Compute geometric factor summed over all D1, D2 times
                // the solid angle of the weighting cube pixel
                double geometry = dri->compute_geometry(oad.tjd(), theta, phi, select, status) *
                                  omega;

                // Compute Earth horizon angle as the distance between the Earth
                // centre in COMPTEL coordinates and the (Chi, Psi) pixel in
                // COMPTEL coordinates minus the radius of the Earth.
                GSkyDir chipsi_comptel;
                chipsi_comptel.radec_deg(phi, 90.0-theta);
                double eha = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

                // Loop over Phibar
                for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {

                    // Compute minimum Earth horizon angle for Phibar layer
                    double ehamin = double(iphibar) * dri->m_phibin + zetamin;

                    // Add up weighted geometry if the Earth horizon angle is
                    // equal or larger than the minimum
                    if (eha >= ehamin) {

                        // Get data space index
                        int inx = index + iphibar * dri->m_dri.npix();

                        // Loop over DRWs and multiply geometry with appopriate
                        // rates
                        for (int i = 0; i < size(); ++i) {
                            m_dris[i][inx] += geometry * rates[i];
                        }

                    } // endif: Earth horizon cut passed

                } // endfor: looped over Phibar

            } // endfor: looped over (Chi, Psi)

        } // endfor: looped over Orbit Aspect Data

        // Normalise Phibar layers of all DRWs to unity
        for (int i = 0; i < size(); ++i) {

            // Loop over Phibar layers
            for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {

                // Compute sum in Phibar layer
                double sum    = 0.0;
                int    istart = iphibar * dri->m_dri.npix();
                int    istop  = istart  + dri->m_dri.npix();
                for (int inx = istart; inx < istop; ++inx) {
                    sum += m_dris[i][inx];
                }

                // Normalise by sum
                if (sum != 0.0) {
                    for (int inx = istart; inx < istop; ++inx) {
                        m_dris[i][inx] /= sum;
                    }
                }

            } // endfor: looped over Phibar layers

            // Set the Good Time interval for DRWs
            if (m_dris[i].m_num_used_superpackets > 0) {
                m_dris[i].m_gti = GGti(dri->m_tstart, dri->m_tstop);
            }

        } // endfor: looped over DRWs

    } // endif: container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Data Space container
 *
 * @param[in] chatter Chattiness.
 * @return String containing Data Space container information.
 ***************************************************************************/
std::string GCOMDris::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMDris ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Number of DRIs"));
        result.append(gammalib::str(size()));

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
void GCOMDris::init_members(void)
{
    // Initialise members
    m_dris.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dris Data Space container.
 ***************************************************************************/
void GCOMDris::copy_members(const GCOMDris& dris)
{
    // Copy members
    m_dris = dris.m_dris;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMDris::free_members(void)
{
    // Return
    return;
}
