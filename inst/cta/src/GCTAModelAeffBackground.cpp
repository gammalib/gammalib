/***************************************************************************
 *       GCTAModelAeffBackground.cpp - CTA Aeff background model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2020 by Michael Mayer                               *
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
#include "GCTASupport.hpp"

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
#define G_AEFF_INTEGRAL             "GCTAModelAeffBackground::aeff_integral("\
                                                    "GObservation&, double&)"

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
 * @brief Set spectral model component
 *
 * @param[in] spectral Pointer to spectral model component.
 *
 * Sets the spectral model component of the model.
 ***************************************************************************/
void GCTAModelAeffBackground::spectral(const GModelSpectral* spectral)
{
    // Free spectral model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_spectral != NULL) && (m_spectral != spectral)) {
        delete m_spectral;
    }

    // Clone spectral model component if it exists, otherwise set pointer
    // to NULL
    m_spectral = (spectral != NULL) ? spectral->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set temporal model component
 *
 * @param[in] temporal Pointer to temporal model component.
 *
 * Sets the temporal model component of the model.
 ***************************************************************************/
void GCTAModelAeffBackground::temporal(const GModelTemporal* temporal)
{
    // Free temporal model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_temporal != NULL) && (m_temporal != temporal)) {
        delete m_temporal;
    }

    // Clone temporal model component if it exists, otherwise set pointer
    // to NULL
    m_temporal = (temporal != NULL) ? temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return background rate in units of events MeV\f$^{-1}\f$
 *        s\f$^{-1}\f$ sr\f$^{-1}\f$
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Background rate (events MeV\f$^{-1}\f$ s\f$^{-1}\f$ sr\f$^{-1}\f$).
 *
 * The method returns a real rate, defined by the number of counts per MeV,
 * steradian and ontime.
 *
 * If the @p gradients flag is true the method will also set the parameter
 * gradients of the model parameters.
 ***************************************************************************/
double GCTAModelAeffBackground::eval(const GEvent&       event,
                                     const GObservation& obs,
                                     const bool&         gradients) const
{
    // Get reference on CTA pointing and effective area response from
    // observation and reference on CTA instrument direction from event
    const GCTAPointing& pnt  = gammalib::cta_pnt(G_EVAL, obs);
    const GCTAAeff&     aeff = gammalib::cta_rsp_aeff(G_EVAL, obs);
    const GCTAInstDir&  dir  = gammalib::cta_dir(G_EVAL, event);

    // Set instrument direction
    GCTAInstDir inst_dir = pnt.instdir(dir.dir());

    // Evaluate function
    double logE = event.energy().log10TeV();
    double spat = aeff(logE,
                       inst_dir.theta(),
                       inst_dir.phi(),
                       pnt.zenith(),
                       pnt.azimuth(), false);
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time(), gradients)
                  : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time(), gradients) : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // Optionally compute partial derivatives
    if (gradients) {

        // Multiply factors to spectral gradients
        if (spectral() != NULL) {
            double fact = spat * temp;
            if (fact != 1.0) {
                for (int i = 0; i < spectral()->size(); ++i)
                    (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact );
            }
        }

        // Multiply factors to temporal gradients
        if (temporal() != NULL) {
            double fact = spat * spec;
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
 * @brief Return spatially integrated background rate in units of
 *        events MeV\f$^{-1}\f$ s\f$^{-1}\f$
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated background rate
 *         (events MeV\f$^{-1}\f$ s\f$^{-1}\f$)
 *
 * Spatially integrates the effective area background model for a given
 * measured event energy and event time. The method returns a real rate,
 * defined as the number of counts per MeV and ontime.
 ***************************************************************************/
double GCTAModelAeffBackground::npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const
{
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
 * Draws a sample of events from the background model using a Monte
 * Carlo simulation. The region of interest, the energy boundaries and the
 * good time interval for the sampling will be extracted from the observation
 * argument that is passed to the method. The method also requires a random
 * number generator of type GRan which is passed by reference, hence the
 * state of the random number generator will be changed by the method.
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

        // Retrieve CTA pointing, effective area response and event list
        const GCTAPointing&  pnt    = gammalib::cta_pnt(G_MC, obs);
        const GCTAAeff&      aeff   = gammalib::cta_rsp_aeff(G_MC, obs);
        const GCTAEventList& events = gammalib::cta_event_list(G_MC, obs);

        // Get simulation region
        const GCTARoi&  roi     = events.roi();
        const GEbounds& ebounds = events.ebounds();
        const GGti&     gti     = events.gti();

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
        GEbounds spectral_ebounds(m_n_mc_energies,
                                  ebounds.emin(),
                                  ebounds.emax(),
                                  "LOG");
        GModelSpectralNodes spectral;
        for (int i = 0; i < spectral_ebounds.size(); ++i) {
            GEnergy energy    = spectral_ebounds.elogmean(i);
            double  intensity = aeff_integral(obs, energy.log10TeV());
            double  norm      = m_spectral->eval(energy, events.tstart());
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
                int n_killed_by_roi      = 0;
                #endif

                // Loop over events
                for (int i = 0; i < n_events; ++i) {

                    // Get Monte Carlo event energy from spectral model
                    GEnergy energy = spectral.mc(ebounds.emin(ieng),
                                                 ebounds.emax(ieng),
                                                 times[i],
                                                 ran);

                    // Get maximum effective area for rejection method
                    double max_aeff = aeff.max(energy.log10TeV(), pnt.zenith(),
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
                        double value = aeff(energy.log10TeV(), offset, phi,
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

                    // Rotate pointing direction by offset and azimuth angle
                    GSkyDir skydir = pnt.dir();
                    skydir.rotate_deg(phi    * gammalib::rad2deg,
                                      offset * gammalib::rad2deg);

                    // Convert rotated pointing direction in instrument system
                    GCTAInstDir mc_dir = pnt.instdir(skydir);

                    // Allocate event
                    GCTAEventAtom event;

                    // Set event attributes
                    event.dir(mc_dir);
                    event.energy(energy);
                    event.time(times[i]);

                    // Append event to list if it falls in ROI
                    if (events.roi().contains(event)) {
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
 * @param[in] chatter Chattiness.
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

    // Get reference on CTA pointing, effective area response and event list
    // from observation
    const GCTAPointing&  pnt    = gammalib::cta_pnt(G_AEFF_INTEGRAL, obs);
    const GCTAAeff&      aeff   = gammalib::cta_rsp_aeff(G_AEFF_INTEGRAL, obs);
    const GCTAEventList& events = gammalib::cta_event_list(G_AEFF_INTEGRAL, obs);

    // Get instrument direction of RoI centre
    GCTAInstDir roi_centre = pnt.instdir(events.roi().centre().dir());

    // Get ROI radius in radians
    double roi_radius = events.roi().radius() * gammalib::deg2rad;

    // Setup integration function
    GCTAModelAeffBackground::npred_roi_kern_theta integrand(&aeff,
                                                            logE,
                                                            roi_centre,
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
                                                              m_roi_centre,
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
 * @brief Kernel for azimuth angle integration of effective area background
 *        model
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
    // Compute detx and dety in radians
    double detx = m_roi_centre.detx();
    double dety = m_roi_centre.dety();
    if (m_theta > 0.0 ) {
        detx += m_theta * std::cos(phi);
        dety += m_theta * std::sin(phi);
    }

    // Convert into theta and phi
    GCTAInstDir dir(detx, dety);
    double theta_prime = dir.theta();
    double phi_prime   = dir.phi();

    // Get background value
    double value = (*m_aeff)(m_logE, theta_prime, phi_prime);

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
