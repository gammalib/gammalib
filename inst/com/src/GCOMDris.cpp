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
#include <fstream>
#include <iostream>
#include "GException.hpp"
#include "GMath.hpp"
#include "GFits.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GOptimizerPars.hpp"
#include "GOptimizerPar.hpp"
#include "GOptimizerLM.hpp"
#include "GFft.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMDri.hpp"
#include "GCOMDris.hpp"
#include "GCOMOad.hpp"
#include "GCOMOads.hpp"
#include "GCOMTim.hpp"
#include "GCOMStatus.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventList.hpp"
#include "GCOMEventAtom.hpp"
#include "GCOMSelection.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMDris::at(int&)"
#define G_INSERT                           "GCOMDris::insert(int&, GCOMDri&)"
#define G_REMOVE                                     "GCOMDris::remove(int&)"
#define G_COMPUTE_DRWS             "GCOMDris::compute_drws(GCOMObservation&,"\
                           " GCOMSelection&, double&, double&, std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_COMPUTE_DRWS_EHACORR  //!< Correct Earth horizon cut

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_COMPUTE_DRWS
//#define G_DEBUG_SAVE_WORKING_ARRAYS

/* __ Constants __________________________________________________________ */


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
 * @brief Compute background weighting cubes
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 * @param[in] timebin Time binning for rate determination (sec).
 * @param[in] method Method to compute DRWs ("energy", "phibar" or "vetorate").
 *
 * @exception GException::invalid_value
 *            DRW with mismatching definition encountered in container.
 * @exception GException::invalid_argument
 *            No event list found in COMPTEL observation.
 *
 * Compute DRW cubes for a COMPTEL observation. DRW cubes are event-rate
 * weighted geometry cubes which were multiplied by the solid angle of
 * the (Chi,Psi) pixels.
 *
 * The DRW cubes are normalised so that the sum of each Phibar layer equals
 * to unity.
 ***************************************************************************/
void GCOMDris::compute_drws(const GCOMObservation& obs,
                            const GCOMSelection&   select,
                            const double&          zetamin,
                            const double&          timebin,
                            const std::string&     method)
{
    // Continue only if there are DRWs in container
    if (!is_empty()) {

        // Debug
        #if defined(G_DEBUG_COMPUTE_DRWS)
        std::cout << "Initialisation" << std::endl;
        std::cout << "--------------" << std::endl;
        #endif

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

            // Set DRW header method parameter
            m_dris[i].m_drw_method = method;

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

        // Method A: energy-dependent weighting
        if (method == "energy") {
            compute_drws_energy(obs, events, select, zetamin, timebin);
        }

        // Method B: Phibar-dependent weighting
        else if (method == "phibar") {
            compute_drws_phibar(obs, events, select, zetamin, timebin);
        }

        // Method C: Veto rate weighting
        else if (method == "vetorate") {
            compute_drws_vetorate(obs, events, select, zetamin);
        }

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

    } // endif: there were DRWs in container

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

    // Initialise working arrays
    m_wrk_counts.clear();
    m_wrk_ehacutcorr.clear();
    m_wrk_vetorate.clear();
    m_wrk_activrate.clear();
    m_wrk_rate.clear();
    m_wrk_use_sp.clear();

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

    // Copy working arrays
    m_wrk_counts     = dris.m_wrk_counts;
    m_wrk_ehacutcorr = dris.m_wrk_ehacutcorr;
    m_wrk_vetorate   = dris.m_wrk_vetorate;
    m_wrk_activrate  = dris.m_wrk_activrate;
    m_wrk_rate       = dris.m_wrk_rate;
    m_wrk_use_sp     = dris.m_wrk_use_sp;

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


/***********************************************************************//**
 * @brief Compute background weighting cubes using energy dependent rates
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] events COMPTEL event list.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 * @param[in] timebin Time binning for rate determination (sec).
 *
 * Compute DRW cubes for a COMPTEL observation. DRW cubes are event-rate
 * weighted geometry cubes which were multiplied by the solid angle of
 * the (Chi,Psi) pixels.
 *
 * For this method the event rates depend on time and energy. Three
 * hard-wired energy ranges are defined for which the event rate is
 * determined:
 *
 *       E < 1.8 MeV
 *       1.8 < E < 4.3 MeV
 *       E > 4.3 MeV
 *
 * For each superpacket, the event rates are linearly interpolated to the
 * relevant superpacket time.
 *
 * The DRW cubes are normalised so that the sum of each Phibar layer equals
 * to unity.
 ***************************************************************************/
void GCOMDris::compute_drws_energy(const GCOMObservation& obs,
                                   const GCOMEventList*   events,
                                   const GCOMSelection&   select,
                                   const double&          zetamin,
                                   const double&          timebin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::compute_drws_energy" << std::endl;
    std::cout << "=============================" << std::endl;
    #endif

    // Initialise D1 & D2 module status
    GCOMStatus status;

    // Get pointer to first DRW in container (for convenient access to members
    // of GCOMDri; we won't use it for modifications, hence we can declare it
    // as constant)
    const GCOMDri* dri = &(m_dris[0]);

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
    bool                set = false;

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

        // Apply Earth horizon cut. Note that we need to correct later
        // for this cut to remove the additional time variation coming
        // from the cut.
        if ((event->eha() - event->phibar()) < zetamin) {
            continue;
        }

        // If time exceeds time bin then add rates
        while (event->time().secs() > (time + timebin)) {

            // If rates are not all zero then append rates to vectors
            if (set) {

                // Append time at mid-point of interval
                times.append(time + 0.5 * timebin);

                // Append rates
                rates0.push_back(rate0);
                rates1.push_back(rate1);
                rates2.push_back(rate2);

                // Debug
                #if defined(G_DEBUG_COMPUTE_DRWS)
                std::cout << "t=" << time + 0.5 * timebin;
                std::cout << " r0=" << rate0;
                std::cout << " r1=" << rate1;
                std::cout << " r2=" << rate2;
                std::cout << std::endl;
                #endif

                // Initialise rates for new accumulation
                rate0 = 0.0;
                rate1 = 0.0;
                rate2 = 0.0;
                set   = false;

            } // endif: there were rates to append

            // Increment time
            time += timebin;

        } // endwhile: event time was after the current integration time stop

        // Accumulate rates
        double energy = event->energy().MeV();
        if (energy < 1.8) {
            rate0 += 1.0;
            set    = true;
        }
        else if (energy < 4.3) {
            rate1 += 1.0;
            set    = true;
        }
        else {
            rate2 += 1.0;
            set    = true;
        }

    } // endfor: looped over events

    // Append remaining rates
    if (set) {

        // Append time at mid-point of interval
        times.append(time + 0.5 * timebin);

        // Append rates
        rates0.push_back(rate0);
        rates1.push_back(rate1);
        rates2.push_back(rate2);

    } // endif: appended remaining rates

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Compute DRWs" << std::endl;
    std::cout << "------------" << std::endl;
    #endif

    // Allocate geometry factor and Earth horizon angle working arrays
    std::vector<double> geometry(dri->m_dri.npix());
    std::vector<double> eha(dri->m_dri.npix());

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

        // Compute solid-angle-corrected geometry factors and Earth horizon
        // angles for all (Chi,Psi) pixels
        for (int index = 0; index < dri->m_dri.npix(); ++index) {

            // Get sky direction and solid angle for (Chi, Psi)
            GSkyDir sky   = dri->m_dri.inx2dir(index);
            double  omega = dri->m_dri.solidangle(index);

            // Convert sky direction to COMPTEL coordinates
            double theta = oad.theta(sky);
            double phi   = oad.phi(sky);

            // Compute geometric factor summed over all D1, D2 times
            // the solid angle of the weighting cube pixel
            geometry[index] = dri->compute_geometry(oad.tjd(), theta, phi,
                                                    select, status) * omega;

            // Compute Earth horizon angle as the distance between the Earth
            // centre in COMPTEL coordinates and the (Chi, Psi) pixel in
            // COMPTEL coordinates minus the radius of the Earth.
            GSkyDir chipsi_comptel;
            chipsi_comptel.radec_deg(phi, 90.0-theta);
            eha[index] = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

        } // endfor: computed solid-angle-corrected geometry factors

        // Loop over all DRWs
        for (int i = 0; i < size(); ++i) {

            // Compute Earth horizon cut correction factor
            double rates_sum_all = 0.0;
            double rates_sum_cut = 0.0;
            double ehacorr       = 1.0;
            for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {

                // Sum product of geometry and rates for all pixels
                for (int index = 0; index < dri->m_dri.npix(); ++index) {
                    rates_sum_all += geometry[index] * rates[i];
                }

                // Sum product of geometry and rates for pixels passing
                // the Earth horizon cut
                double ehamin = double(iphibar) * dri->m_phibin + zetamin;
                for (int index = 0; index < dri->m_dri.npix(); ++index) {
                    if (eha[index] >= ehamin) {
                        rates_sum_cut += geometry[index] * rates[i];
                    }
                }

            } // endfor: looped over Phibar

            // If something passes the Earth horizon cut then compute the
            // correction factor now
            if (rates_sum_cut != 0.0) {
                ehacorr = rates_sum_all / rates_sum_cut;
            }

            // Loop over all Phibar layers
            for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {

                // Compute minimum Earth horizon angle for Phibar layer
                double ehamin = double(iphibar) * dri->m_phibin + zetamin;

                // Loop over all (Chi,Psi) pixels
                for (int index = 0; index < dri->m_dri.npix(); ++index) {

                    // If Earth horizon angle is equal or larger than the minimum
                    // requirement then add up the rate weighted geometry factor
                    if (eha[index] >= ehamin) {

                        // Get data space index
                        int inx = index + iphibar * dri->m_dri.npix();

                        // Store product as DRW value
                        m_dris[i][inx] += geometry[index] * rates[i] * ehacorr;

                    } // endif: Earth horizon cut passed

                } // endfor: looped over (Chi, Psi)

            } // endfor: looped over Phibar

        } // endfor: looped over all DRWs

    } // endfor: looped over Orbit Aspect Data

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute background weighting cubes using Phibar dependent rates
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] events COMPTEL event list.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 * @param[in] timebin Time binning for rate determination (sec).
 *
 * Compute DRW cubes for a COMPTEL observation. DRW cubes are event-rate
 * weighted geometry cubes which were multiplied by the solid angle of
 * the (Chi,Psi) pixels.
 *
 * For this method the event rates depend on time and Phibar angle. For
 * each superpacket, the event rates are linearly interpolated to the
 * relevant superpacket time.
 *
 * The DRW cubes are normalised so that the sum of each Phibar layer equals
 * to unity.
 ***************************************************************************/
void GCOMDris::compute_drws_phibar(const GCOMObservation& obs,
                                   const GCOMEventList*   events,
                                   const GCOMSelection&   select,
                                   const double&          zetamin,
                                   const double&          timebin)
{
    // Set constants
    const double ehacorrmax = 10.0; //!< Maximum EHA correction factor

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::compute_drws_phibar" << std::endl;
    std::cout << "=============================" << std::endl;
    #endif

    // Initialise D1 & D2 module status
    GCOMStatus status;

    // Get pointer to first DRW in container (for convenient access to members
    // of GCOMDri; we won't use it for modifications, hence we can declare it
    // as constant)
    const GCOMDri* dri = &(m_dris[0]);

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Determine Phibar-dependent event rates" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    #endif

    // Initialise current time and vectors for computation of node function
    // and vectors for energy-dependent rate interpolation
    double                            time;
    std::vector<double>               rate(dri->nphibar());
    GNodeArray                        times;
    std::vector<std::vector<double> > rates(dri->nphibar());
    bool                              set = false;

    // Set time of first event
    if (events->size() > 0) {
        time = (*events)[0]->time().secs();
    }

    // Determine energy-dependent event rates. The standard event selection
    // criteria and earth horizon cut is applied so that the event rate
    // corresponds to the current selection, yet we do not consider any
    // Good Time Intervals or OAD time information as we later only access
    // relevant times.
    for (int i_evt = 0; i_evt < events->size(); ++i_evt) {

        // Get pointer to event
        const GCOMEventAtom* event = (*events)[i_evt];

        // Apply event selection
        if (!select.use_event(*event)) {
            continue;
        }

        // Apply Earth horizon cut. Note that we need to correct later
        // for this cut to remove the additional time variation coming
        // from the cut.
        #if defined(G_COMPUTE_DRWS_EHACORR)
        if ((event->eha() - event->phibar()) < zetamin) {
            continue;
        }
        #endif

        // If time exceeds time bin then add rates
        while (event->time().secs() > (time + timebin)) {

            // If rates are set then append rates to vectors
            if (set) {

                // Append time at mid-point of interval
                times.append(time + 0.5 * timebin);

                // Append rates and initialise for new accumulation
                for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {
                    rates[iphibar].push_back(rate[iphibar]);
                    rate[iphibar] = 0.0;
                }

                // Signal that rates were not yet set
                set = false;

            } // endif: there were rates to append

            // Increment time
            time += timebin;

        } // endwhile: event time was after the current integration time stop

        // Accumulate rates
        int iphibar = int((event->phibar() - dri->phimin()) / dri->phibin());
        if (iphibar >= 0 and iphibar < dri->nphibar()) {
            rate[iphibar] += 1.0;
            set            = true;
        }

    } // endfor: looped over events

    // Append remaining rates
    if (set) {

        // Append time at mid-point of interval
        times.append(time + 0.5 * timebin);

        // Append rates
        for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {
            rates[iphibar].push_back(rate[iphibar]);
        }

    } // endif: appended remaining rates

    // Debug: write debug file
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::ofstream debugfile;
    debugfile.open("compute_drws_phibar_tbin_rates.csv");
    debugfile.precision(12);
    for (int i = 0; i < times.size(); ++i) {
        debugfile << times[i];
        for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {
            debugfile << "," << rates[iphibar][i];
        }
        debugfile << std::endl;
    }
    debugfile.close();
    #endif

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Compute DRWs" << std::endl;
    std::cout << "------------" << std::endl;
    #endif

    // Allocate geometry factor and Earth horizon angle working arrays
    std::vector<double> geometry(dri->m_dri.npix());
    std::vector<double> eha(dri->m_dri.npix());

    // Get Good Time Intervals. If a pulsar selection is specified then
    // reduce the Good Time Intervals to the validity intervals of the
    // pulsar emphemerides.
    GCOMTim tim = obs.tim();
    if (select.has_pulsar()) {
        tim.reduce(select.pulsar().validity());
    }

    // Debug: open debug files file
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::ofstream debugfile1;
    std::ofstream debugfile2;
    debugfile1.open("compute_drws_phibar_sp_rates.csv");
    debugfile2.open("compute_drws_phibar_sp_ehacorr.csv");
    debugfile1.precision(12);
    debugfile2.precision(12);
    #endif

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

        // Compute Phibar-dependent rates by interpolation of precomputed
        // rate vectors. Make sure that rates never become negative.
        std::vector<double> rate_oad(dri->nphibar());
        for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {
            rate_oad[iphibar] = times.interpolate(time, rates[iphibar]);
            if (rate_oad[iphibar] < 0.0) {
                rate_oad[iphibar] = 0.0;
            }
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

        // Compute solid-angle-corrected geometry factors and Earth horizon
        // angles for all (Chi,Psi) pixels
        for (int index = 0; index < dri->m_dri.npix(); ++index) {

            // Get sky direction and solid angle for (Chi, Psi)
            GSkyDir sky   = dri->m_dri.inx2dir(index);
            double  omega = dri->m_dri.solidangle(index);

            // Convert sky direction to COMPTEL coordinates
            double theta = oad.theta(sky);
            double phi   = oad.phi(sky);

            // Compute geometric factor summed over all D1, D2 times
            // the solid angle of the weighting cube pixel
            geometry[index] = dri->compute_geometry(oad.tjd(), theta, phi,
                                                    select, status) * omega;

            // Compute Earth horizon angle as the distance between the Earth
            // centre in COMPTEL coordinates and the (Chi, Psi) pixel in
            // COMPTEL coordinates minus the radius of the Earth.
            GSkyDir chipsi_comptel;
            chipsi_comptel.radec_deg(phi, 90.0-theta);
            eha[index] = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

        } // endfor: computed solid-angle-corrected geometry factors

        // Debug: write time in debug files
        #if defined(G_DEBUG_COMPUTE_DRWS)
        debugfile1 << time;
        debugfile2 << time;
        #endif

        // Correct rates for Earth horizon cut
        for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {

            // Compute minimum Earth Horizon angle for Phibar layer
            double ehamin = double(iphibar) * dri->m_phibin + zetamin;

            // Initialise sums
            double rates_sum_all = 0.0;
            double rates_sum_cut = 0.0;

            // Sum geometry for all pixels and pixels passing the Earth horizon cut
            for (int index = 0; index < dri->m_dri.npix(); ++index) {
                rates_sum_all += geometry[index];
                if (eha[index] >= ehamin) {
                    rates_sum_cut += geometry[index];
                }
            }

            // Compute Earth horizon cut correction, limited to a maximum of ehacorrmax
            double ehacorr = (rates_sum_cut != 0.0) ? rates_sum_all / rates_sum_cut : ehacorrmax;
            if (ehacorr > ehacorrmax) {
                ehacorr = ehacorrmax;
            }

            // Correct rate
            #if defined(G_COMPUTE_DRWS_EHACORR)
            rate_oad[iphibar] *= ehacorr;
            #endif

            // Debug: write rate and correction in debug files
            #if defined(G_DEBUG_COMPUTE_DRWS)
            debugfile1 << "," << rate_oad[iphibar];
            debugfile2 << "," << ehacorr;
            #endif

        } // endfor: looped over Phibar

        // Debug: write linefeeds in debug files
        #if defined(G_DEBUG_COMPUTE_DRWS)
        debugfile1 << std::endl;
        debugfile2 << std::endl;
        #endif

        // Loop over all DRWs
        for (int i = 0; i < size(); ++i) {

            // Loop over all Phibar layers
            for (int iphibar = 0; iphibar < dri->nphibar(); ++iphibar) {

                // Compute minimum Earth horizon angle for Phibar layer
                double ehamin = double(iphibar) * dri->m_phibin + zetamin;

                // Loop over all (Chi,Psi) pixels
                for (int index = 0; index < dri->m_dri.npix(); ++index) {

                    // If Earth horizon angle is equal or larger than the minimum
                    // requirement then add up the rate weighted geometry factor
                    if (eha[index] >= ehamin) {

                        // Get data space index
                        int inx = index + iphibar * dri->m_dri.npix();

                        // Store product as DRW value
                        m_dris[i][inx] += geometry[index] * rate_oad[iphibar];

                    } // endif: Earth horizon cut passed

                } // endfor: looped over (Chi, Psi)

            } // endfor: looped over Phibar

        } // endfor: looped over all DRWs

    } // endfor: looped over Orbit Aspect Data

    // Debug: close debug files
    #if defined(G_DEBUG_COMPUTE_DRWS)
    debugfile1.close();
    debugfile2.close();
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute background weighting cubes using veto rates
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] events COMPTEL event list.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 *
 * Compute DRW cubes for a COMPTEL observation. DRW cubes are event-rate
 * weighted geometry cubes which were multiplied by the solid angle of
 * the (Chi,Psi) pixels.
 ***************************************************************************/
void GCOMDris::compute_drws_vetorate(const GCOMObservation& obs,
                                     const GCOMEventList*   events,
                                     const GCOMSelection&   select,
                                     const double&          zetamin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::compute_drws_vetorate" << std::endl;
    std::cout << "===============================" << std::endl;
    #endif

    // Setup working arrays
    vetorate_setup(obs, events, select, zetamin);

    // Debug: Save working arrays in FITS file
    #if defined(G_DEBUG_SAVE_WORKING_ARRAYS)
    vetorate_save("debug_gcomdris_compute_drws_vetorate.fits");
    vetorate_load("debug_gcomdris_compute_drws_vetorate.fits");
    #endif

    // Initialise fprompt for all energies
    for (int i = 0; i < size(); ++i) {
        m_dris[i].m_drw_fprompt = 0.5; //!< Initial fprompt value
    }

    // Fit working arrays to determine f_prompt parameter
    vetorate_fit();
    for (int iter = 0; iter < 4; ++iter) { //!< Perform 4 iterations
        vetorate_update_activ();
        vetorate_fit();
    }

    // Generate DRWs
    vetorate_generate(obs, select, zetamin);

    // Finish DRWs
    vetorate_finish(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup working arrays for vetorate computation
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] events COMPTEL event list.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 *
 * Setup working arrays for vetorate DRW computation. There are five working
 * arrays:
 *
 * @c m_wrk_counts is a 3D working array spanned by superpacket index,
 * phibar layer and energy bin, containing the number of events passing the
 * event selections for the specified bin.
 *
 * @c m_wrk_ehacutcorr is a 2D working array spanned by superpacket index
 * and phibar layer, containing the fraction of the solid-angle weighted
 * geometry that passes the Earth horizon cut for the specified bin. This
 * array represents the response to a given incident event rate.
 *
 * @c m_wrk_vetorate is a 1D working array, containing the veto rate as
 * function of superpacket index.
 *
 * @c m_wrk_activrate is a 2D working array, containing the activation rate
 * as function of superpacket index and energy range.
 *
 * @c m_wrk_use_sp is a 1D working array, signalling the usage of a given
 * superpacket.
 *
 * Note that the method also sets the superpacket usage and event selection
 * statistics of all DRWs.
 ***************************************************************************/
void GCOMDris::vetorate_setup(const GCOMObservation& obs,
                              const GCOMEventList*   events,
                              const GCOMSelection&   select,
                              const double&          zetamin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_setup" << std::endl;
    std::cout << "========================" << std::endl;
    #endif

    // Initialise D1 & D2 module status
    GCOMStatus status;

    // Get pointer to first DRW (assume their definition is identical)
    const GCOMDri* dri0 = &(m_dris[0]);

    // Extract dimensions
    int noads   = obs.oads().size();
    int npix    = dri0->m_dri.npix();
    int nphibar = dri0->nphibar();
    int neng    = size();

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Number of OAD superpackets: " << noads << std::endl;
    std::cout << "Number of pixels: " << npix << std::endl;
    std::cout << "Number of Phibar layers: " << nphibar << std::endl;
    std::cout << "Number of energies: " << neng << std::endl;
    #endif

    // Define working arrays
    m_wrk_counts     = GNdarray(noads, nphibar, neng);
    m_wrk_ehacutcorr = GNdarray(noads, nphibar);
    m_wrk_vetorate   = GNdarray(noads);
    m_wrk_activrate  = GNdarray(noads, neng);
    m_wrk_use_sp     = std::vector<bool>(noads);

    // Initialise superpacket usage array
    for (int i_oad = 0; i_oad < noads; ++i_oad) {
        m_wrk_use_sp[i_oad] = false;
    }

    // Allocate further working arrays
    GNdarray geo(npix);
    GNdarray eha(npix);

    // Setup node array and rate vector for veto rate interpolation
    const GCOMHkd&      hdk = obs.hkds()["SCV2M"];
    GNodeArray          veto_times;
    std::vector<double> veto_rates;
    veto_times.reserve(hdk.size());
    veto_rates.reserve(hdk.size());
    for (int i = 0; i < hdk.size(); ++i) {
        double rate = hdk.value(i);
        if (rate > 0.0) {
            veto_times.append(hdk.time(i).secs());
            veto_rates.push_back(rate);
        }
    }

    // Setup sky directions and solid angles for all (Chi, Psi)
    std::vector<GSkyDir> skys;
    std::vector<double>  omegas;
    for (int index = 0; index < dri0->m_dri.npix(); ++index) {
        GSkyDir sky   = dri0->m_dri.inx2dir(index);
        double  omega = dri0->m_dri.solidangle(index);
        skys.push_back(sky);
        omegas.push_back(omega);
    }

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Setup node array for veto rate interpolation: ";
    std::cout << veto_times.size() << " times, ";
    std::cout << veto_rates.size() << " rates" << std::endl;
    #endif

    // Get Good Time Intervals. If a pulsar selection is specified then
    // reduce the Good Time Intervals to the validity intervals of the
    // pulsar emphemerides.
    GCOMTim tim = obs.tim();
    if (select.has_pulsar()) {
        tim.reduce(select.pulsar().validity());
    }

    // Initialise event and superpacket counters
    int i_evt = 0;
    int i_sp  = 0;

    // Signal that loop should be terminated
    bool terminate = false;

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < noads; ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = obs.oads()[i_oad];

        // Check superpacket usage for all DRWs so that their statistics
        // get correctly updated
        for (int i = 0; i < size(); ++i) {
            if (m_dris[i].use_superpacket(oad, tim, select)) {
                m_wrk_use_sp[i_oad] = true;
            }
        }
        if (!m_wrk_use_sp[i_oad]) {
            continue;
        }

        // Get Orbit Aspect Data mid-time in seconds
        double time = oad.tstart().secs() + 0.5 * (oad.tstop() - oad.tstart());

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

        // Compute solid-angle-corrected geometry factors and Earth horizon
        // angles for all (Chi,Psi) pixels
        double total_geo = 0.0;
        for (int index = 0; index < npix; ++index) {

            // Convert sky direction to COMPTEL coordinates
            double theta = oad.theta(skys[index]);
            double phi   = oad.phi(skys[index]);

            // Compute geometric factor summed over all D1, D2 times
            // the solid angle of the weighting cube pixel
            double geometry = dri0->compute_geometry(oad.tjd(), theta, phi,
                                                     select, status) * omegas[index];

            // Store and sum up geometric factor
            geo(index)  = geometry;
            total_geo  += geometry;

            // Compute Earth horizon angle as the distance between the Earth
            // centre in COMPTEL coordinates and the (Chi, Psi) pixel in
            // COMPTEL coordinates minus the radius of the Earth.
            GSkyDir chipsi_comptel;
            chipsi_comptel.radec_deg(phi, 90.0-theta);
            eha(index) = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

        } // endfor: looped over (Chi, Psi)

        // Compute Earth horizon cut correction factor for all Phibar layers
        if (total_geo > 0.0) {
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
                double ehamin       = double(iphibar) * dri0->m_phibin + zetamin;
                double accepted_geo = 0.0;
                for (int index = 0; index < npix; ++index) {
                    if (eha(index) >= ehamin) {
                        accepted_geo += geo(index);
                    }
                }
                m_wrk_ehacutcorr(i_sp, iphibar) = accepted_geo / total_geo;
            }
        }

        // Compute veto rate and activation rates for superpacket time if rates
        // are available. The activation rates are initialised with a constant
        // rate.
        if (!veto_times.is_empty()) {
            double value = veto_times.interpolate(time, veto_rates);
            if (value > 0.0) {
                m_wrk_vetorate(i_sp) = value;
                for (int ieng = 0; ieng < neng; ++ieng) {
                    m_wrk_activrate(i_sp, ieng) = 1.0;
                }
            }
        }

        // Collect all events within superpacket and store result in 3D array.
        // Break if the end of the event list was reached.
        for (; i_evt < events->size(); ++i_evt) {

            // Get pointer to event
            const GCOMEventAtom* event = (*events)[i_evt];

            // Break loop if the end of the superpacket was reached
            if (event->time() > oad.tstop()) {
                break;
            }

            // Determine energy bin
            int ienergy = -1;
            for (int i = 0; i < size(); ++i) {
                if ((event->energy() >= m_dris[i].m_ebounds.emin()) &&
                    (event->energy() <= m_dris[i].m_ebounds.emax())) {
                    ienergy = i;
                    break;
                }
            }

            // Skip event if it is not contained in any energy bin
            if (ienergy == -1) {
                continue;
            }

            // Get pointer to appropriate DRW
            const GCOMDri* dri = &(m_dris[ienergy]);

            // Check GTIs if the DRW has GTIs
            if (dri->gti().size() > 0) {

                // Skip event if it lies before the GTI start
                if (event->time() < dri->gti().tstart()) {
                    continue;
                }

                // Break Orbit Aspect Data loop if event lies after the GTI stop
                else if (event->time() > dri->gti().tstop()) {
                    terminate = true;
                    break;
                }

            } // endif: DRW had GTIs

            // Skip event if it lies before the superpacket start
            if (event->time() < oad.tstart()) {
                continue;
            }

            // Apply event selection
            if (!dri->m_selection.use_event(*event)) {
                continue;
            }

            // Compute Compton scatter angle index. Skip if it's invalid.
            int iphibar = int((event->phibar() - dri->m_phimin) / dri->m_phibin);
            if (iphibar < 0) {
                continue;
            }
            else if (iphibar >= dri->nphibar()) {
                continue;
            }

            // Check for Earth horizon angle. There is a constant EHA limit
            // over a Phibar layer to be compliant with the DRG.
            double ehamin = double(iphibar) * dri->m_phibin + zetamin;
            if (event->eha() < ehamin) {
                continue;
            }

            // Now fill the event in the 3D counts array
            m_wrk_counts(i_sp, iphibar, ienergy) += 1.0;

        } // endfor: collected events

        // Debug
        #if defined(G_DEBUG_COMPUTE_DRWS)
        std::cout << "Superpacket " << i_oad;
        std::cout << " i_sp=" << i_sp;
        std::cout << " time=" << time;
        std::cout << " vetorate=" << m_wrk_vetorate(i_oad);
        std::cout << " total_geo=" << total_geo;
        std::cout << std::endl;
        #endif

        // Increment superpacket counter
        i_sp++;

        // Break Orbit Aspect Data loop if termination was signalled or if there
        // are no more events
        if (terminate || i_evt >= events->size()) {
            break;
        }

    } // endfor: looped over Orbit Aspect Data

    // Copy results in shortened Ndarrays
    GNdarray wrk_counts(i_sp, nphibar, neng);
    GNdarray wrk_ehacutcorr(i_sp, nphibar);
    GNdarray wrk_vetorate(i_sp);
    GNdarray wrk_activrate(i_sp, neng);
    for (int i = 0; i < i_sp; ++i) {
        for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
            for (int ieng = 0; ieng < neng; ++ieng) {
                wrk_counts(i, iphibar, ieng) = m_wrk_counts(i, iphibar, ieng);
            }
            wrk_ehacutcorr(i, iphibar) = m_wrk_ehacutcorr(i, iphibar);
        }
        for (int ieng = 0; ieng < neng; ++ieng) {
            wrk_activrate(i, ieng) = m_wrk_activrate(i, ieng);
        }
        wrk_vetorate(i) = m_wrk_vetorate(i);
    }
    m_wrk_counts     = wrk_counts;
    m_wrk_ehacutcorr = wrk_ehacutcorr;
    m_wrk_vetorate   = wrk_vetorate;
    m_wrk_activrate  = wrk_activrate;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fit working arrays for vetorate computation
 *
 * Determines the f_prompt parameter for each energy bin using a maximum
 * log-likelihood fit using the data in the working arrays for each energy
 * bin. Based on the fitted f_prompt parameter the method also computes the
 * expected background rate time variation for each energy bin.
 ***************************************************************************/
void GCOMDris::vetorate_fit(void)
{
    // Set constants
    const double norm = 1.0; // Model normalisation

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_fit" << std::endl;
    std::cout << "======================" << std::endl;
    #endif

    // Extract dimension from working arrays
    int nsp     = m_wrk_counts.shape()[0];
    int nphibar = m_wrk_counts.shape()[1];
    int neng    = m_wrk_counts.shape()[2];

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Number of used superpackets: " << nsp << std::endl;
    std::cout << "Number of Phibar layers: " << nphibar << std::endl;
    std::cout << "Number of energies: " << neng << std::endl;
    #endif

    // Setup rate result array
    m_wrk_rate = GNdarray(nsp, neng);

    // Normalise veto rates
    double nvetorate = 0.0;
    for (int isp = 0; isp < nsp; ++isp) {
        if (m_wrk_vetorate(isp) > 0.0) {
            nvetorate += m_wrk_vetorate(isp);
        }
    }
    if (nvetorate > 0.0) {
        for (int isp = 0; isp < nsp; ++isp) {
            m_wrk_vetorate(isp) /= nvetorate;
        }
    }

    // Normalise activation rates
    for (int ieng = 0; ieng < neng; ++ieng) {
        double nactivrate = 0.0;
        for (int isp = 0; isp < nsp; ++isp) {
            if (m_wrk_activrate(isp, ieng) > 0.0) {
                nactivrate += m_wrk_activrate(isp, ieng);
            }
        }
        if (nactivrate > 0.0) {
            for (int isp = 0; isp < nsp; ++isp) {
                m_wrk_activrate(isp, ieng) /= nactivrate;
            }
        }
    }

    // Loop over energy bins
    for (int ieng = 0; ieng < neng; ++ieng) {

        // Set fprompt parameter
        GOptimizerPars pars;
        GOptimizerPar  par("fprompt", m_dris[ieng].m_drw_fprompt);
        par.range(0.0, 1.0);
        par.factor_gradient(1.0); // To avoid zero gradient log message
        pars.append(par);

        // Set fit function
        GCOMDris::likelihood fct(this, ieng, norm);

        // Optimize function and compute errors
        #if defined(G_DEBUG_COMPUTE_DRWS)
        GLog log;
        log.cout(true);
        GOptimizerLM opt(&log);
        #else
        GOptimizerLM opt;
        #endif
        opt.eps(0.1);
        opt.max_iter(500);
        opt.optimize(fct, pars);
        opt.errors(fct, pars);

        // Retrieve logL
        double logL = opt.value();

        // Recover fprompt parameter and errors
        double fprompt   = pars[0]->value();
        double e_fprompt = pars[0]->error();
        double factiv    = 1.0 - fprompt;

        // Set resulting rate vector for this energy bin
        for (int isp = 0; isp < nsp; ++isp) {
            m_wrk_rate(isp, ieng) = fprompt * m_wrk_vetorate(isp) +
                                    factiv  * m_wrk_activrate(isp, ieng);
        }

        // Set DRW header members
        m_dris[ieng].m_drw_status    = opt.status_string();
        m_dris[ieng].m_drw_fprompt   = fprompt;
        m_dris[ieng].m_drw_e_fprompt = e_fprompt;
        m_dris[ieng].m_drw_iter      = opt.iter();

        // Debug: Print fit results
        #if defined(G_DEBUG_COMPUTE_DRWS)
        std::cout << ieng;
        std::cout << " logL=" << logL;
        std::cout << " fprompt=" << fprompt << " +/- " << e_fprompt;
        std::cout << std::endl;
        #endif

    } // endfor: looped over energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update activation rate
 ***************************************************************************/
void GCOMDris::vetorate_update_activ(void)
{
    // Set constants
    const double norm   = 1.0; // Model normalisation
    const double t_hour = 1.0; // Smoothing duration (hours)

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_update_activ" << std::endl;
    std::cout << "===============================" << std::endl;
    #endif

    // Extract dimension from working arrays
    int nsp     = m_wrk_counts.shape()[0];
    int nphibar = m_wrk_counts.shape()[1];
    int neng    = m_wrk_counts.shape()[2];

    // Setup Gaussian smoothing kernel
    double   nsmooth = 3600.0/16.384 * t_hour;
    double   weight  = -0.5 / (nsmooth * nsmooth);
    GNdarray kernel(nsp);
    double   sum     = 0.0;
    for (int isp = 0; isp < nsp/2; ++isp) {
        double kvalue = std::exp(weight*isp*isp);
        if (kvalue <= 0.0) {
            break;
        }
        kernel(isp) = kvalue;
        sum        += kvalue;
        if (isp > 0) {
            kernel(nsp-isp) = kvalue;
            sum            += kvalue;
        }
    }
    if (sum > 0.0) {
        for (int isp = 0; isp < nsp; ++isp) {
            kernel(isp) /= sum;
        }
    }
    GFft fft_kernel(kernel);

    // Loop over energy bins
    for (int ieng = 0; ieng < neng; ++ieng) {

        // Compute model normalised to number of events for each Phibar
        GNdarray model(nsp, nphibar);
        GNdarray n(nphibar);
        for (int isp = 0; isp < nsp; ++isp) {
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
                model(isp, iphibar) = m_wrk_rate(isp, ieng) * m_wrk_ehacutcorr(isp, iphibar);
            }
        }
        for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
            double nevents = 0.0;
            double nmodel  = 0.0;
            for (int isp = 0; isp < nsp; ++isp) {
                if (model(isp, iphibar) > 0.0) {
                    nevents += m_wrk_counts(isp, iphibar, ieng);
                    nmodel  += model(isp, iphibar);
                }
            }
            if (nmodel > 0.0) {
                n(iphibar) = norm * (nevents / nmodel);
                for (int isp = 0; isp < nsp; ++isp) {
                    model(isp, iphibar) *= n(iphibar);
                }
            }
        }

        // Compute residual and Earth horizon correction cut by summing over
        // all Phibar layers
        GNdarray residual(nsp);
        GNdarray ehacutcorr(nsp);
        for (int isp = 0; isp < nsp; ++isp) {
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
                residual(isp)   += m_wrk_counts(isp, iphibar, ieng) - model(isp, iphibar);
                ehacutcorr(isp) += n(iphibar) * m_wrk_ehacutcorr(isp, iphibar);
            }
        }

        // Smooth residual
        GFft fft_residual(residual);
        fft_residual = fft_residual * fft_kernel;
        residual     = fft_residual.backward();

        // Smooth ehacutcorr
        GFft fft_ehacutcorr(ehacutcorr);
        fft_ehacutcorr = fft_ehacutcorr * fft_kernel;
        ehacutcorr     = fft_ehacutcorr.backward();

        // Derive indicent residual background rate
        for (int isp = 0; isp < nsp; ++isp) {
            if (ehacutcorr(isp) > 0.0) {
                residual(isp) /= ehacutcorr(isp);
            }
        }

        // Update activation rate
        for (int isp = 0; isp < nsp; ++isp) {
            double activ = m_wrk_activrate(isp, ieng);
            if (activ > 0.0) {
                double update = activ + residual(isp);
                if (update < 0.05 * activ) { //!< Cap reduction at 5%
                    update = 0.05 * activ;
                }
                m_wrk_activrate(isp, ieng) = update;
            }
        }

        // Normalise activation rate
        double nactivrate = 0.0;
        for (int isp = 0; isp < nsp; ++isp) {
            nactivrate += m_wrk_activrate(isp, ieng);
        }
        if (nactivrate > 0.0) {
            for (int isp = 0; isp < nsp; ++isp) {
                m_wrk_activrate(isp, ieng) /= nactivrate;
            }
        }

    } // endfor: looped over all energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate DRWs
 *
 * @param[in] obs COMPTEL observation.
 * @param[in] select Selection set.
 * @param[in] zetamin Minimum Earth horizon - Phibar cut (deg).
 ***************************************************************************/
void GCOMDris::vetorate_generate(const GCOMObservation& obs,
                                 const GCOMSelection&   select,
                                 const double&          zetamin)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_generate" << std::endl;
    std::cout << "===========================" << std::endl;
    #endif

    // Initialise D1 & D2 module status
    GCOMStatus status;

    // Get pointer to first DRW (assume their definition is identical)
    const GCOMDri* dri0 = &(m_dris[0]);

    // Extract dimensions
    int noads   = obs.oads().size();
    int npix    = dri0->m_dri.npix();
    int nphibar = dri0->nphibar();
    int neng    = size();

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Number of OAD superpackets: " << noads << std::endl;
    std::cout << "Number of pixels: " << npix << std::endl;
    std::cout << "Number of Phibar layers: " << nphibar << std::endl;
    std::cout << "Number of energies: " << neng << std::endl;
    #endif

    // Allocate geometry factor and Earth horizon angle working arrays
    GNdarray geometry(npix);
    GNdarray eha(npix);

    // Setup sky directions and solid angles for all (Chi, Psi)
    std::vector<GSkyDir> skys;
    std::vector<double>  omegas;
    for (int index = 0; index < npix; ++index) {
        GSkyDir sky   = dri0->m_dri.inx2dir(index);
        double  omega = dri0->m_dri.solidangle(index);
        skys.push_back(sky);
        omegas.push_back(omega);
    }

    // Get Good Time Intervals. If a pulsar selection is specified then
    // reduce the Good Time Intervals to the validity intervals of the
    // pulsar emphemerides.
    GCOMTim tim = obs.tim();
    if (select.has_pulsar()) {
        tim.reduce(select.pulsar().validity());
    }

    // Initialise superpacket counter
    int i_sp = 0;

    // Loop over Orbit Aspect Data
    for (int i_oad = 0; i_oad < noads; ++i_oad) {

        // Get reference to Orbit Aspect Data of superpacket
        const GCOMOad &oad = obs.oads()[i_oad];

        // Skip not-to-be-used superpackets
        if (!m_wrk_use_sp[i_oad]) {
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

        // Compute solid-angle-corrected geometry factors and Earth horizon
        // angles for all (Chi,Psi) pixels
        for (int index = 0; index < npix; ++index) {

            // Convert sky direction to COMPTEL coordinates
            double theta = oad.theta(skys[index]);
            double phi   = oad.phi(skys[index]);

            // Compute geometric factor summed over all D1, D2 times
            // the solid angle of the weighting cube pixel
            geometry(index) = dri0->compute_geometry(oad.tjd(), theta, phi,
                                                     select, status) * omegas[index];

            // Compute Earth horizon angle as the distance between the Earth
            // centre in COMPTEL coordinates and the (Chi, Psi) pixel in
            // COMPTEL coordinates minus the radius of the Earth.
            GSkyDir chipsi_comptel;
            chipsi_comptel.radec_deg(phi, 90.0-theta);
            eha(index) = geocentre_comptel.dist_deg(chipsi_comptel) - oad.georad();

        } // endfor: computed solid-angle-corrected geometry factors

        // Loop over all DRWs
        for (int i = 0; i < size(); ++i) {

            // Get pointer to appropriate DRW
            const GCOMDri* dri = &(m_dris[i]);

            // Loop over all Phibar layers
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

                // Compute minimum Earth horizon angle for Phibar layer
                double ehamin = double(iphibar) * dri->m_phibin + zetamin;

                // Loop over all (Chi,Psi) pixels
                for (int index = 0; index < npix; ++index) {

                    // If Earth horizon angle is equal or larger than the minimum
                    // requirement then add up the rate weighted geometry factor
                    if (eha(index) >= ehamin) {

                        // Get data space index
                        int inx = index + iphibar * npix;

                        // Store geometry weighted with rate as DRW
                        m_dris[i][inx] += geometry(index) * m_wrk_rate(i_sp, i);

                    } // endif: Earth horizon cut passed

                } // endfor: looped over (Chi, Psi)

            } // endfor: looped over Phibar

            // Set binary table information

        } // endfor: looped over all DRWs

        // Increment superpacket counter
        i_sp++;

    } // endfor: looped over Orbit Aspect Data

    // Return
    return;
}


/***********************************************************************//**
 * @brief Finish DRWs
 *
 * @param[in] obs COMPTEL observation.
 *
 * Set DRW binary tables
 ***************************************************************************/
void GCOMDris::vetorate_finish(const GCOMObservation& obs)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_finish" << std::endl;
    std::cout << "=========================" << std::endl;
    #endif

    // Extract dimensions
    int noads = obs.oads().size();

    // Compute number of to-be-used superpackets
    int nsp = 0;
    for (int i_oad = 0; i_oad < noads; ++i_oad) {
        if (m_wrk_use_sp[i_oad]) {
            nsp++;
        }
    }

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Number of OAD superpackets: " << noads << std::endl;
    std::cout << "Number of used superpackets: " << nsp << std::endl;
    #endif

    // Setup DRW binary tables
    for (int i = 0; i < size(); ++i) {

        // Get pointer to DRW
        GCOMDri* dri = &(m_dris[i]);

        // Allocate columns
        GFitsTableLongCol   tjd("TJD",  nsp);
        GFitsTableLongCol   tics("TICS", nsp);
        GFitsTableDoubleCol rate("RATE", nsp);
        GFitsTableDoubleCol vetorate("VETORATE", nsp);
        GFitsTableDoubleCol activrate("ACTIVRATE", nsp);

        // Initialise superpacket counter
        int isp = 0;

        // Loop over Orbit Aspect Data
        for (int i_oad = 0; i_oad < noads; ++i_oad) {

            // Get reference to Orbit Aspect Data of superpacket
            const GCOMOad &oad = obs.oads()[i_oad];

            // Skip not-to-be-used superpackets
            if (!m_wrk_use_sp[i_oad]) {
                continue;
            }

            // Fill table columns
            tjd(isp)       = oad.tjd();
            tics(isp)      = oad.tics();
            rate(isp)      = m_wrk_rate(isp, i);
            vetorate(isp)  = m_wrk_vetorate(isp);
            activrate(isp) = m_wrk_activrate(isp, i);

            // Increment superpacket counter
            isp++;

        } // endfor: looped over Orbit Aspect Data

        // Append columns to binary table
        dri->m_drw_table.clear();
        dri->m_drw_table.extname("RATES");
        dri->m_drw_table.append(tjd);
        dri->m_drw_table.append(tics);
        dri->m_drw_table.append(rate);
        dri->m_drw_table.append(vetorate);
        dri->m_drw_table.append(activrate);

    } // endfor: looped over DRWs

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save working arrays for vetorate computation
 *
 * @param[in] filename FITS filename.
 *
 * Save working arrays for vetorate DRW computation in FITS file.
 ***************************************************************************/
void GCOMDris::vetorate_save(const GFilename& filename) const
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_save" << std::endl;
    std::cout << "=======================" << std::endl;
    #endif

    // Extract dimension from working arrays
    int noads   = m_wrk_use_sp.size();
    int nsp     = m_wrk_counts.shape()[0];
    int nphibar = m_wrk_counts.shape()[1];
    int neng    = m_wrk_counts.shape()[2];

    // Allocate FITS images
    GFitsImageDouble image_counts(nsp, nphibar, neng);
    GFitsImageDouble image_ehacutcorr(nsp, nphibar);
    GFitsImageDouble image_vetorate(nsp);
    GFitsImageDouble image_activrate(nsp, neng);
    GFitsImageShort  image_use_sp(noads);

    // Fill FITS images from working arrays
    for (int isp = 0; isp < nsp; ++isp) {
        for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
            for (int ieng = 0; ieng < neng; ++ieng) {
                image_counts(isp, iphibar, ieng) = m_wrk_counts(isp, iphibar, ieng);
            }
            image_ehacutcorr(isp, iphibar) = m_wrk_ehacutcorr(isp, iphibar);
        }
        for (int ieng = 0; ieng < neng; ++ieng) {
            image_activrate(isp, ieng) = m_wrk_activrate(isp, ieng);
        }
        image_vetorate(isp) = m_wrk_vetorate(isp);
    }
    for (int ioad = 0; ioad < noads; ++ioad) {
        image_use_sp(ioad) = (m_wrk_use_sp[ioad]) ? 1 : 0;
    }

    // Set image attributes
    image_counts.extname("COUNTS");
    image_ehacutcorr.extname("EHACUTCORR");
    image_vetorate.extname("VETORATE");
    image_activrate.extname("ACTIVRATE");
    image_use_sp.extname("SPUSAGE");

    // Create FITS file and append images
    GFits fits;
    fits.append(image_counts);
    fits.append(image_ehacutcorr);
    fits.append(image_vetorate);
    fits.append(image_activrate);
    fits.append(image_use_sp);

    // Save FITS file
    fits.saveto(filename, true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load working arrays for vetorate computation
 *
 * @param[in] filename FITS filename.
 *
 * Load working arrays for vetorate DRW computation from FITS file.
 ***************************************************************************/
void GCOMDris::vetorate_load(const GFilename& filename)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "GCOMDris::vetorate_load" << std::endl;
    std::cout << "=======================" << std::endl;
    #endif

    // Open FITS file
    GFits fits(filename);

    // Get images
    const GFitsImage* image_counts     = fits.image("COUNTS");
    const GFitsImage* image_ehacutcorr = fits.image("EHACUTCORR");
    const GFitsImage* image_vetorate   = fits.image("VETORATE");
    const GFitsImage* image_activrate  = fits.image("ACTIVRATE");
    const GFitsImage* image_use_sp     = fits.image("SPUSAGE");

    // Extract dimensions
    int noads   = image_use_sp->naxes(0);
    int nsp     = image_counts->naxes(0);
    int nphibar = image_counts->naxes(1);
    int neng    = image_counts->naxes(2);

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRWS)
    std::cout << "Number of OAD superpackets: " << noads << std::endl;
    std::cout << "Number of used superpackets: " << nsp << std::endl;
    std::cout << "Number of Phibar layers: " << nphibar << std::endl;
    std::cout << "Number of energies: " << neng << std::endl;
    #endif

    // Allocate working array
    m_wrk_counts     = GNdarray(nsp, nphibar, neng);
    m_wrk_ehacutcorr = GNdarray(nsp, nphibar);
    m_wrk_vetorate   = GNdarray(nsp);
    m_wrk_activrate  = GNdarray(nsp, neng);
    m_wrk_use_sp     = std::vector<bool>(noads);

    // Extract data from images
    for (int isp = 0; isp < nsp; ++isp) {
        for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
            for (int ieng = 0; ieng < neng; ++ieng) {
                m_wrk_counts(isp, iphibar, ieng) = image_counts->pixel(isp, iphibar, ieng);
            }
            m_wrk_ehacutcorr(isp, iphibar) = image_ehacutcorr->pixel(isp, iphibar);
        }
        for (int ieng = 0; ieng < neng; ++ieng) {
            m_wrk_activrate(isp, ieng) = image_activrate->pixel(isp, ieng);
        }
        m_wrk_vetorate(isp) = image_vetorate->pixel(isp);
    }
    for (int ioad = 0; ioad < noads; ++ioad) {
        m_wrk_use_sp[ioad] = (image_use_sp->pixel(ioad) == 1);
    }

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Log-likelihood function constructor
 *
 * @param[in] dris Parent calling likelihood function.
 * @param[in] ieng DRW energy bin.
 * @param[in] norm Normalisation constant.
 *
 * Constructs the log-likelihood function that is used to determine the
 * value of f_prompt.
 ***************************************************************************/
GCOMDris::likelihood::likelihood(GCOMDris     *dris,
                                 const int&    ieng,
                                 const double& norm) : GOptimizerFunction()
{
    // Store arguments
    m_this = dris;
    m_ieng = ieng;
    m_norm = norm;

    // Set gradient and curvature members
    m_gradient  = GVector(1);
    m_curvature = GMatrixSparse(1,1);

    // Extract dimension from working array
    m_nsp     = m_this->m_wrk_counts.shape()[0];
    m_nphibar = m_this->m_wrk_counts.shape()[1];

    // Multiply veto and constant rate by EHA cut correction. Rates are zero
    // in case that the vetorate is zero.
    m_vetorate  = GNdarray(m_nsp, m_nphibar);
    m_activrate = GNdarray(m_nsp, m_nphibar);
    m_diffrate  = GNdarray(m_nsp, m_nphibar);
    for (int isp = 0; isp < m_nsp; ++isp) {
        double vetorate = m_this->m_wrk_vetorate(isp);
        if (vetorate > 0.0) {
            double activrate = m_this->m_wrk_activrate(isp, ieng);
            for (int iphibar = 0; iphibar < m_nphibar; ++iphibar) {
                m_vetorate(isp, iphibar)  = vetorate  * m_this->m_wrk_ehacutcorr(isp, iphibar);
                m_activrate(isp, iphibar) = activrate * m_this->m_wrk_ehacutcorr(isp, iphibar);
                m_diffrate(isp, iphibar)  = m_vetorate(isp, iphibar) - m_activrate(isp, iphibar);
            }
        }
    }

    // Compute rate sums
    m_vetorate_sum  = GNdarray(m_nphibar);
    m_activrate_sum = GNdarray(m_nphibar);
    m_diffrate_sum  = GNdarray(m_nphibar);
    for (int iphibar = 0; iphibar < m_nphibar; ++iphibar) {
        for (int isp = 0; isp < m_nsp; ++isp) {
            m_vetorate_sum(iphibar)  += m_vetorate(isp, iphibar);
            m_activrate_sum(iphibar) += m_activrate(isp, iphibar);
        }
        m_diffrate_sum(iphibar) = m_vetorate_sum(iphibar) - m_activrate_sum(iphibar);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Log-likelihood function evaluation
 *
 * @param[in] pars Function parameters.
 *
 * Computes the log-likelihood function, its gradient and its curvature at
 * the specified function parameters.
 ***************************************************************************/
void GCOMDris::likelihood::eval(const GOptimizerPars& pars)
{
    // Recover parameters
    double fprompt = pars[0]->value();
    double factiv  = 1.0 - fprompt;

    // Allocate phibar normalisation arrays
    GNdarray nevents(m_nphibar);
    GNdarray nmodel(m_nphibar);
    GNdarray norm(m_nphibar);

    // Compute model
    GNdarray model = fprompt * m_vetorate + factiv * m_activrate;

    // Compute Phibar normalised model
    GNdarray model_norm = model;
    for (int iphibar = 0; iphibar < m_nphibar; ++iphibar) {
        for (int isp = 0; isp < m_nsp; ++isp) {
            if (model(isp, iphibar) > 0.0) {
                nevents(iphibar) += m_this->m_wrk_counts(isp, iphibar, m_ieng);
                nmodel(iphibar)  += model(isp, iphibar);
            }
        }
        if (nmodel(iphibar) > 0.0) {
            norm(iphibar) = m_norm * nevents(iphibar) / nmodel(iphibar);
            for (int isp = 0; isp < m_nsp; ++isp) {
                model_norm(isp, iphibar) *= norm(iphibar);
            }
        }
    }

    // Compute LogL
    m_value = 0.0;
    for (int isp = 0; isp < m_nsp; ++isp) {
        for (int iphibar = 0; iphibar < m_nphibar; ++iphibar) {
            if (model_norm(isp, iphibar) > 0.0) {
                m_value -= m_this->m_wrk_counts(isp, iphibar, m_ieng) *
                           std::log(model_norm(isp, iphibar)) - model_norm(isp, iphibar);
            }
        }
    }

    // Evaluate gradient and curvature
    for (int iphibar = 0; iphibar < m_nphibar; ++iphibar) {

        // Precompute normalisation gradient
        double nmodel2 = nmodel(iphibar) * nmodel(iphibar);
        double g_norm  = -m_norm * nevents(iphibar) / nmodel2 * m_diffrate_sum(iphibar);

        // Loop over superpackets
        for (int isp = 0; isp < m_nsp; ++isp) {

            // Continue only if model is positive
            if (model_norm(isp, iphibar) > 0.0) {

                // Pre computation
                double fb = m_this->m_wrk_counts(isp, iphibar, m_ieng) / model_norm(isp, iphibar);
                double fc = (1.0 - fb);
                double fa = fb / model_norm(isp, iphibar);

                // Compute parameter gradient
                double g = norm(iphibar) * m_diffrate(isp, iphibar) + g_norm * model(isp, iphibar);

                // Update gradient
                m_gradient[0] += fc * g;

                // Update curvature
                m_curvature(0,0) += fa * g * g;

            } // endif: model was positive

        } // endfor: looped over superpackets

    } // endfor: looped over Phibar

    // Return
    return;
}
