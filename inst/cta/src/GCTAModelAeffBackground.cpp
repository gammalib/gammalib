/***************************************************************************
 *       GCTAModelAeffBackground.cpp - CTA Aeff background model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Michael Mayer                               *
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
 * @file GCTAModelAeffBackground.cpp
 * @brief CTA Aeff background model class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelSpectralNodes.hpp"
#include "GCTAModelAeffBackground.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTABackground.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAModelAeffBackground g_cta_aeff_background_seed;
const GModelRegistry          g_cta_aeff_background_registry(&g_cta_aeff_background_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL "GCTAModelAeffBackground::eval(GEvent&, GObservation&, bool&)"
#define G_NPRED            "GCTAModelAeffBackground::npred(GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC              "GCTAModelAeffBackground::mc(GObservation&, GRan&)"
#define G_XML_SPECTRAL  "GCTAModelAeffBackground::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL  "GCTAModelAeffBackground::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_NPRED_CACHE

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC
//#define G_DEBUG_NPRED

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelAeffBackground::GCTAModelAeffBackground(void) : GModelData()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA effective area background model from the information
 * provided by an XML element. The XML element is expected to have the
 * following structure
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * Optionally, a temporal model may be provided using the following
 * syntax
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporalModel type="...">
 *         ...
 *       </temporalModel>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
GCTAModelAeffBackground::GCTAModelAeffBackground(const GXmlElement& xml) :
                         GModelData(xml)
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from spectral component
 *
 * @param[in] spectral Spectral model component.
 *
 * Constructs a CTA effective area background model from a @p spectral
 * model component. The temporal component is assumed to be constant.
 ***************************************************************************/
GCTAModelAeffBackground::GCTAModelAeffBackground(const GModelSpectral& spectral) :
                         GModelData()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from model components
 *
 * @param[in] spectral Spectral model component.
 * @param[in] temporal Temporal model component.
 *
 * Constructs a CTA effective area background model from a @p spectral and a
 * @p temporal component.
 ***************************************************************************/
GCTAModelAeffBackground::GCTAModelAeffBackground(const GModelSpectral& spectral,
                                                 const GModelTemporal& temporal) :
                         GModelData()
{
    // Initialise members
    init_members();

    // Clone model components
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] aeff CTA instrument background model.
 ***************************************************************************/
GCTAModelAeffBackground::GCTAModelAeffBackground(const GCTAModelAeffBackground& aeff) :
                         GModelData(aeff)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(aeff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelAeffBackground::~GCTAModelAeffBackground(void)
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
 * @param[in] aeff CTA effective area background model.
 * @return CTA effective area background model.
 ***************************************************************************/
GCTAModelAeffBackground& GCTAModelAeffBackground::operator=(const GCTAModelAeffBackground& aeff)
{
    // Execute only if object is not identical
    if (this != &aeff) {

        // Copy base class members
        this->GModelData::operator=(aeff);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(aeff);

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
 * @brief Clear CTA effective area background model
 *
 * This method properly resets the CTA effective area background model to an
 * initial state.
 ***************************************************************************/
void GCTAModelAeffBackground::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelData::free_members();

    // Initialise members
    this->GModelData::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone CTA effective area background model
 *
 * @return Pointer to deep copy of CTA effective area background model.
 ***************************************************************************/
GCTAModelAeffBackground* GCTAModelAeffBackground::clone(void) const
{
    return new GCTAModelAeffBackground(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            Specified observation is not of the expected type.
 *
 * If the @p gradients flag is true the method will also set the parameter
 * gradients of the model parameters.
 *
 * @todo Make sure that DETX and DETY are always set in GCTAInstDir.
 ***************************************************************************/
double GCTAModelAeffBackground::eval(const GEvent&       event,
                                     const GObservation& obs,
                                     const bool&         gradients) const
{
    // Get pointer on CTA observation
    const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);
    if (cta == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Get pointer on CTA IRF response
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(cta->response());
    if (rsp == NULL) {
        std::string msg = "Specified observation does not contain an IRF response.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Retrieve pointer to CTA Effective Area
    const GCTAAeff* aeff = rsp->aeff();
    if (aeff == NULL) {
        std::string msg = "Specified observation contains no effective area"
                          " information.\n" + obs.print();
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Extract CTA instrument direction from event
    const GCTAInstDir* dir  = dynamic_cast<const GCTAInstDir*>(&(event.dir()));
    if (dir == NULL) {
        std::string msg = "No CTA instrument direction found in event.";
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Set instrument direction
    GCTAInstDir inst_dir = cta->pointing().instdir(dir->dir());

    // Evaluate function
    double logE = event.energy().log10TeV();
    double spat = (*aeff)(logE,
                          inst_dir.theta(),
                          inst_dir.phi(),
                          cta->pointing().zenith(),
                          cta->pointing().azimuth(), false);
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time(), gradients)
                  : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time(), gradients) : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // Apply deadtime correction
    double deadc = obs.deadc(event.time());
    value       *= deadc;

    // Optionally compute partial derivatives
    if (gradients) {

        // Multiply factors to spectral gradients
        if (spectral() != NULL) {
            double fact = spat * temp * deadc;
            if (fact != 1.0) {
                for (int i = 0; i < spectral()->size(); ++i)
                    (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact );
            }
        }

        // Multiply factors to temporal gradients
        if (temporal() != NULL) {
            double fact = spat * spec * deadc;
            if (fact != 1.0) {
                for (int i = 0; i < temporal()->size(); ++i)
                    (*temporal())[i].factor_gradient((*temporal())[i].factor_gradient() * fact );
            }
        }

    } // endif: computed partial derivatives

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated background model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated model.
 *
 * @exception GException::invalid_argument
 *            The specified observation is not a CTA observation.
 *
 * Spatially integrates the effective area background model for a given
 * measured event energy and event time. This method also applies a deadtime
 * correction factor, so that the normalization of the model is a real rate
 * (counts/MeV/s).
 ***************************************************************************/
double GCTAModelAeffBackground::npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    //static const int iter_theta = 6;
    //static const int iter_phi   = 6;

    // Initialise result
    double npred     = 0.0;
    bool   has_npred = false;

    // Build unique identifier
    std::string id = obs.instrument() + "::" + obs.id();

    // Check if Npred value is already in cache
    #if defined(G_USE_NPRED_CACHE)
    if (!m_npred_names.empty()) {

        // Search for unique identifier, and if found, recover Npred value
        // and break
        for (int i = 0; i < m_npred_names.size(); ++i) {
            if (m_npred_names[i] == id && m_npred_energies[i] == obsEng) {
                npred     = m_npred_values[i];
                has_npred = true;
                #if defined(G_DEBUG_NPRED)
                std::cout << "GCTAModelAeffBackground::npred:";
                std::cout << " cache=" << i;
                std::cout << " npred=" << npred << std::endl;
                #endif
                break;
            }
        }

    } // endif: there were values in the Npred cache
    #endif

    // Continue only if no Npred cache value has been found
    if (!has_npred) {

        // Evaluate only if model is valid
        if (valid_model()) {

            // Get log10 of energy in TeV
            double logE = obsEng.log10TeV();

            // Spatially integrate effective area component
            npred = this->aeff_integral(obs, logE);

            // Store result in Npred cache
            #if defined(G_USE_NPRED_CACHE)
            m_npred_names.push_back(id);
            m_npred_energies.push_back(obsEng);
            m_npred_times.push_back(obsTime);
            m_npred_values.push_back(npred);
            #endif

            // Debug: Check for NaN
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
                std::string origin  = "GCTAModelAeffBackground::npred";
                std::string message = " NaN/Inf encountered (npred=" +
                                      gammalib::str(npred) + ")";
                gammalib::warning(origin, message);
            }
            #endif

        } // endif: model was valid

    } // endif: Npred computation required

    // Multiply in spectral and temporal components
    npred *= spectral()->eval(obsEng, obsTime);
    npred *= temporal()->eval(obsTime);

    // Apply deadtime correction
    npred *= obs.deadc(obsTime);

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return Pointer to list of simulated events (needs to be de-allocated by
 *         client)
 *
 * @exception GException::invalid_argument
 *            Specified observation is not a CTA observation.
 *
 * Draws a sample of events from the background model using a Monte
 * Carlo simulation. The region of interest, the energy boundaries and the
 * good time interval for the sampling will be extracted from the observation
 * argument that is passed to the method. The method also requires a random
 * number generator of type GRan which is passed by reference, hence the
 * state of the random number generator will be changed by the method.
 *
 * The method also applies a deadtime correction using a Monte Carlo process,
 * taking into account temporal deadtime variations. For this purpose, the
 * method makes use of the time dependent GObservation::deadc method.
 *
 * For each event in the returned event list, the sky direction, the nominal
 * coordinates (DETX and DETY), the energy and the time will be set.
 ***************************************************************************/
GCTAEventList* GCTAModelAeffBackground::mc(const GObservation& obs,
                                           GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Retrieve CTA observation
        const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);
        if (cta == NULL) {
            std::string msg = "Specified observation is not a CTA "
                              "observation.\n" + obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Get pointer on CTA IRF response
        const GCTAResponseIrf* rsp =
              dynamic_cast<const GCTAResponseIrf*>(cta->response());
        if (rsp == NULL) {
            std::string msg = "Specified observation does not contain"
                              " an IRF response.\n" + obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Retrieve CTA response and pointing
        const GCTAPointing& pnt = cta->pointing();

        // Get pointer to CTA effective area
        const GCTAAeff* aeff = rsp->aeff();
        if (aeff == NULL) {
            std::string msg = "Specified observation contains no effective area"
                              " information.\n" + obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Retrieve event list to access the ROI, energy boundaries and GTIs
        const GCTAEventList* events =
              dynamic_cast<const GCTAEventList*>(obs.events());
        if (events == NULL) {
            std::string msg = "No CTA event list found in observation.\n" +
                              obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Get simulation region
        const GCTARoi&  roi     = events->roi();
        const GEbounds& ebounds = events->ebounds();
        const GGti&     gti     = events->gti();

        // Get maximum offset value for simulations
        double max_theta     = pnt.dir().dist(roi.centre().dir()) +
                               roi.radius() * gammalib::deg2rad;
        double cos_max_theta = std::cos(max_theta);

        // Set simulation region for result event list
        list->roi(roi);
        list->ebounds(ebounds);
        list->gti(gti);

        // Set up spectral model to draw random energies from. Here we use
        // a fixed energy sampling for an instance of GModelSpectralNodes.
        // This is analogous to to the GCTAModelIrfBackground::mc method.
        // We make sure that only non-negative nodes get appended.
        GEbounds spectral_ebounds =
            GEbounds(m_n_mc_energies, ebounds.emin(), ebounds.emax(), true);
        GModelSpectralNodes spectral;
        for (int i = 0; i < spectral_ebounds.size(); ++i) {
            GEnergy energy    = spectral_ebounds.elogmean(i);
            double  intensity = aeff_integral(obs, energy.log10TeV());
            double  norm      = m_spectral->eval(energy, events->tstart());
            double  arg       = norm * intensity;
            if (arg > 0.0) {
                spectral.append(energy, arg);
            }
        }

        // Loop over all energy bins
        for (int ieng = 0; ieng < ebounds.size(); ++ieng) {

            // Compute the background rate in model within the energy
            // boundaries from spectral component (units: cts/s).
            // Note that the time here is ontime. Deadtime correction will
            // be done later.
            double rate = spectral.flux(ebounds.emin(ieng), ebounds.emax(ieng));

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GCTAModelAeffBackground::mc(\"" << name() << "\": ";
            std::cout << "rate=" << rate << " cts/s)" << std::endl;
            #endif

            // If the rate is not positive then skip this energy bins
            if (rate <= 0.0) {
                continue;
            }

            // Loop over all good time intervals
            for (int itime = 0; itime < gti.size(); ++itime) {

                // Get Monte Carlo event arrival times from temporal model
                GTimes times = m_temporal->mc(rate,
                                              gti.tstart(itime),
                                              gti.tstop(itime),
                                              ran);

                // Get number of events
                int n_events = times.size();

                // Reserve space for events
                if (n_events > 0) {
                    list->reserve(n_events);
                }

                // Debug option: provide number of times and initialize
                // statisics
                #if defined(G_DUMP_MC)
                std::cout << " Interval " << itime;
                std::cout << " times=" << n_events << std::endl;
                int n_killed_by_deadtime = 0;
                int n_killed_by_roi      = 0;
                #endif

                // Loop over events
                for (int i = 0; i < n_events; ++i) {

                    // Apply deadtime correction
                    double deadc = obs.deadc(times[i]);
                    if (deadc < 1.0) {
                        if (ran.uniform() > deadc) {
                            #if defined(G_DUMP_MC)
                            n_killed_by_deadtime++;
                            #endif
                            continue;
                        }
                    }

                    // Get Monte Carlo event energy from spectral model
                    GEnergy energy = spectral.mc(ebounds.emin(ieng),
                                                 ebounds.emax(ieng),
                                                 times[i],
                                                 ran);

                    // Get maximum effective area for rejection method
                    double max_aeff = aeff->max(energy.log10TeV(), pnt.zenith(),
                                                pnt.azimuth(), false);

                    // Skip event if the maximum effective area is not positive
                    if (max_aeff <= 0.0) {
                        continue;
                    }

                    // Initialise randomised coordinates
                    double offset = 0.0;
                    double phi    = 0.0;

                    // Initialise acceptance fraction and counter of zeros for
                    // rejection method
                    double acceptance_fraction = 0.0;
                    int    zeros               = 0;

                    // Start rejection method loop
                    do {

                        // Throw random offset and azimuth angle in
                        // considered range
                        offset = std::acos(1.0 - ran.uniform() *
                                           (1.0 - cos_max_theta));
                        phi    = ran.uniform() * gammalib::twopi;

                        // Compute function value at this offset angle
                        double value = (*aeff)(energy.log10TeV(), offset, phi,
                                               pnt.zenith(), pnt.azimuth(),
                                               false);

                        // If the value is not positive then increment the
                        // zeros counter and fall through. The counter assures
                        // that this loop does not lock up.
                        if (value <= 0.0) {
                            zeros++;
                            continue;
                        }

                        // Value is non-zero so reset the zeros counter
                        zeros = 0;

                        // Compute acceptance fraction
                        acceptance_fraction = value / max_aeff;
			
                    } while ((ran.uniform() > acceptance_fraction) &&
                             (zeros < 1000));

                    // If the zeros counter is non-zero then the loop was
                    // exited due to exhaustion and the event is skipped
                    if (zeros > 0) {
                        continue;
                    }

                    // Convert CTA pointing direction in instrument system
                    GCTAInstDir mc_dir(pnt.dir());

                    // Rotate pointing direction by offset and azimuth angle
                    mc_dir.dir().rotate_deg(phi    * gammalib::rad2deg,
                                            offset * gammalib::rad2deg);

                    // Compute DETX and DETY coordinates
                    double detx(0.0);
                    double dety(0.0);
                    if (offset > 0.0 ) {
                        detx = offset * std::cos(phi);
                        dety = offset * std::sin(phi);
                    }

                    // Set DETX and DETY coordinates
                    mc_dir.detx(detx);
                    mc_dir.dety(dety);

                    // Allocate event
                    GCTAEventAtom event;

                    // Set event attributes
                    event.dir(mc_dir);
                    event.energy(energy);
                    event.time(times[i]);

                    // Append event to list if it falls in ROI
                    if (events->roi().contains(event)) {
                        list->append(event);
                    }
                    #if defined(G_DUMP_MC)
                    else {
                        n_killed_by_roi++;
                    }
                    #endif

                } // endfor: looped over all events

                // Debug option: provide  statisics
                #if defined(G_DUMP_MC)
                std::cout << " Killed by deadtime=";
                std::cout << n_killed_by_deadtime << std::endl;
                std::cout << " Killed by ROI=";
                std::cout << n_killed_by_roi << std::endl;
                #endif

            } // endfor: looped over all GTIs

        } // endfor: looped over all energy boundaries

    } // endif: model was valid

    // Return
    return list;
}


/***********************************************************************//**
 * @brief Read CTA effective area background model from XML element
 *
 * @param[in] xml XML element.
 *
 * Set up CTA effective area background model from the information provided by
 * an XML element. The XML element is expected to have the following
 * structure
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * Optionally, a temporal model may be provided using the following
 * syntax
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporalModel type="...">
 *         ...
 *       </temporalModel>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
void GCTAModelAeffBackground::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Initialise XML elements
    const GXmlElement* spectral = NULL;
    const GXmlElement* temporal = NULL;

    // Get pointer on spectrum
    spectral = xml.element("spectrum", 0);

    // Extract spectral model
    m_spectral = xml_spectral(*spectral);

    // Optionally get temporal model
    if (xml.elements("temporalModel")) {
        temporal   = xml.element("temporalModel", 0);
        m_temporal = xml_temporal(*temporal);
    }
    else {
        GModelTemporalConst constant;
        m_temporal = constant.clone();
    }

    // Read model attributes
    read_attributes(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA effective area background model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write CTA effective area background model information into an XML element.
 * The XML element will have the following structure
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * If the model contains a non-constant temporal model, the temporal
 * component will also be written following the syntax
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporalModel type="...">
 *         ...
 *       </temporalModel>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
void GCTAModelAeffBackground::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = xml.element("source", k);
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If we have a temporal model that is either not a constant, or a
    // constant with a normalization value that differs from 1.0 then
    // write the temporal component into the XML element. This logic
    // assures compatibility with the Fermi/LAT format as this format
    // does not handle temporal components.
    bool write_temporal = ((m_temporal != NULL) &&
                           (m_temporal->type() != "Constant" ||
                            (*m_temporal)[0].value() != 1.0));

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (write_temporal)     src->append(GXmlElement("temporalModel"));
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Optionally write temporal model
    if (write_temporal) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("temporalModel", 0);
            temporal()->write(*temp);
        }
    }

    // Write model attributes
    write_attributes(*src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA effective area background model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing CTA effective area background model information.
 ***************************************************************************/
std::string GCTAModelAeffBackground::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelAeffBackground ===");

        // Determine number of parameters per type
        int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
        int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type"));
        if (n_spectral > 0) {
            result.append("\""+spectral()->type()+"\"");
            if (n_temporal > 0) {
                result.append(" * ");
            }
        }
        if (n_temporal > 0) {
            result.append("\""+temporal()->type()+"\"");
        }

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of spectral par's") +
                      gammalib::str(n_spectral));
        for (int i = 0; i < n_spectral; ++i) {
            result.append("\n"+(*spectral())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of temporal par's") +
                      gammalib::str(n_temporal));
        for (int i = 0; i < n_temporal; ++i) {
            result.append("\n"+(*temporal())[i].print());
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelAeffBackground::init_members(void)
{
    // Initialise members
    m_spectral      = NULL;
    m_temporal      = NULL;
    m_n_mc_energies = 100;

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
 * @param[in] bgd CTA effective area background model.
 ***************************************************************************/
void GCTAModelAeffBackground::copy_members(const GCTAModelAeffBackground& bgd)
{
    // Copy cache
    m_npred_names    = bgd.m_npred_names;
    m_npred_energies = bgd.m_npred_energies;
    m_npred_times    = bgd.m_npred_times;
    m_npred_values   = bgd.m_npred_values;
    m_n_mc_energies  = bgd.m_n_mc_energies;

    // Clone spectral and temporal model components
    m_spectral = (bgd.m_spectral != NULL) ? bgd.m_spectral->clone() : NULL;
    m_temporal = (bgd.m_temporal != NULL) ? bgd.m_temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelAeffBackground::free_members(void)
{
    // Free memory
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModelData base class.
 ***************************************************************************/
void GCTAModelAeffBackground::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Determine the number of parameters
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i) {
            m_pars.push_back(&((*spectral())[i]));
        }

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i) {
            m_pars.push_back(&((*temporal())[i]));
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 *
 * Returns 'true' if models has a spectral and a temporal component.
 * Otherwise returns 'false'.
 ***************************************************************************/
bool GCTAModelAeffBackground::valid_model(void) const
{
    // Set result
    bool result = ((spectral() != NULL) && (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return pointer to spectral model from XML element
 *
 * @param[in] spectral XML element.
 * @return Pointer to spectral model.
 *
 * Returns pointer to spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelAeffBackground::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(spectral);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return pointer to temporal model from XML element
 *
 * @param[in] temporal XML element.
 * @return Pointer to temporal model.
 *
 * Returns pointer to temporal model that is defined in an XML element.
 ***************************************************************************/
GModelTemporal* GCTAModelAeffBackground::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(temporal);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Spatially integrate effective area for given energy
 *
 * @param[in] obs Observation.
 * @param[in] logE Log10 of reference energy in TeV.
 * @return Spatially integrated effective area for given energy.
 *
 * @exception GException::invalid_argument
 *            Invalid observation encountered.
 *
 * Spatially integrates the effective area for a given reference energy
 * over the region of interest.
 ***************************************************************************/
double GCTAModelAeffBackground::aeff_integral(const GObservation& obs,
                                              const double&       logE) const
{
    // Initialise result
    double value = 0.0;

    // Set number of iterations for Romberg integration.
    static const int iter_theta = 6;
    static const int iter_phi   = 6;

    // Get pointer on CTA observation
    const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);
    if (cta == NULL) {
        std::string msg = "Specified observation is not a CTA"
                          " observation.\n" + obs.print();
        throw GException::invalid_argument(G_NPRED, msg);
    }

    // Get pointer on CTA IRF response
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(cta->response());
    if (rsp == NULL) {
        std::string msg = "Specified observation does not contain"
                          " an IRF response.\n" + obs.print();
        throw GException::invalid_argument(G_NPRED, msg);
    }

    // Retrieve pointer to CTA effective area
    const GCTAAeff* aeff = rsp->aeff();
    if (aeff == NULL) {
        std::string msg = "Specified observation contains no effective area"
                          " information.\n" + obs.print();
        throw GException::invalid_argument(G_NPRED, msg);
    }

    // Get CTA event list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
    if (events == NULL) {
        std::string msg = "No CTA event list found in observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_NPRED, msg);
    }

    // Get ROI radius in radians
    double roi_radius = events->roi().radius() * gammalib::deg2rad;

    // Setup integration function
    GCTAModelAeffBackground::npred_roi_kern_theta integrand(aeff,
                                                            logE,
                                                            iter_phi);

    // Setup integration
    GIntegral integral(&integrand);

    // Set fixed number of iterations
    integral.fixed_iter(iter_theta);

    // Spatially integrate radial component
    value = integral.romberg(0.0, roi_radius);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::string origin  = "GCTAModelAeffBackground::aeff_integral";
        std::string message = " NaN/Inf encountered (value=" +
                              gammalib::str(value) + ", roi_radius=" +
                              gammalib::str(roi_radius) + ")";
        gammalib::warning(origin, message);
    }
    #endif

    // Return
    return value;

}


/***********************************************************************//**
 * @brief Kernel for offset angle integration of effective area background model
 *
 * @param[in] theta Offset angle from ROI centre (radians).
 * @return Integration kernel value.
 *
 * Computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \theta \times
 *                     \int_{0}^{2\pi}
 *                     B(\theta,\phi | E, t) d\phi
 * \f]
 *
 * where \f$B(\theta,\phi | E, t)\f$ is the background model for a specific
 * observed energy \f$E\f$ and time \f$t\f$.
 ***************************************************************************/
double GCTAModelAeffBackground::npred_roi_kern_theta::eval(const double& theta)
{
    // Initialise value
    double value = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Setup phi integration kernel
        GCTAModelAeffBackground::npred_roi_kern_phi integrand(m_aeff,
                                                              m_logE,
                                                              theta);

        // Setup integration
        GIntegral integral(&integrand);

        // Set fixed number of iterations
        integral.fixed_iter(m_iter);

        // Integrate over phi
        value = integral.romberg(0.0, gammalib::twopi) * std::sin(theta);

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::string origin  = "GCTAModelAeffBackground::npred_roi_kern_theta::eval"
                                  "(" + gammalib::str(theta) + ")";
            std::string message = " NaN/Inf encountered (value=" +
                                  gammalib::str(value) + ")";
            gammalib::warning(origin, message);
        }
        #endif

    } // endif: offset angle was positive

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle integration of effective area background model
 *
 * @param[in] phi Azimuth angle around ROI centre (radians).
 * @return Integration kernel value.
 *
 * Computes
 *
 * \f[
 *    B(\theta, \phi | E, t)
 * \f]
 ***************************************************************************/
double GCTAModelAeffBackground::npred_roi_kern_phi::eval(const double& phi)
{
    // Get background value
    double value = (*m_aeff)(m_logE, m_theta, phi);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::string origin  = "GCTAModelAeffBackground::npred_roi_kern_phi::eval"
                              "(" + gammalib::str(phi) + ")";
        std::string message = " NaN/Inf encountered (value=" +
                              gammalib::str(value) + ", theta=" +
                              gammalib::str(m_theta) + ")";
        gammalib::warning(origin, message);
    }
    #endif

    // Return Npred
    return value;
}
