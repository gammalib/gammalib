/***************************************************************************
 *       GCTAModelIrfBackground.cpp - CTA IRF background model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTAModelIrfBackground.cpp
 * @brief CTA IRF background model class implementation
 * @author Juergen Knoedlseder
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
#include "GCTAModelIrfBackground.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTABackground.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAModelIrfBackground g_cta_inst_background_seed;
const GModelRegistry         g_cta_inst_background_registry(&g_cta_inst_background_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL         "GCTAModelIrfBackground::eval(GEvent&, GObservation&)"
#define G_EVAL_GRADIENTS    "GCTAModelIrfBackground::eval_gradients(GEvent&,"\
                                                            " GObservation&)"
#define G_NPRED             "GCTAModelIrfBackground::npred(GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC               "GCTAModelIrfBackground::mc(GObservation&, GRan&)"
#define G_XML_SPECTRAL   "GCTAModelIrfBackground::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL   "GCTAModelIrfBackground::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_NPRED_CACHE

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC
//#define G_DEBUG_NPRED

/* __ Constants __________________________________________________________ */
const double g_cta_inst_background_npred_theta_eps = 1.0e-4;
const double g_cta_inst_background_npred_phi_eps   = 1.0e-4;


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelIrfBackground::GCTAModelIrfBackground(void) : GModelData()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA instrumental background model from the information
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
GCTAModelIrfBackground::GCTAModelIrfBackground(const GXmlElement& xml) :
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
 * Constructs a CTA instrumental background model from a spectral
 * model component. The temporal component is assumed to be constant.
 * Please refer to the classe GModelSpectral to learn more about the
 * definition of the spectral components.
 ***************************************************************************/
GCTAModelIrfBackground::GCTAModelIrfBackground(const GModelSpectral& spectral) :
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
 * @brief Copy constructor
 *
 * @param[in] bgd CTA instrument background model.
 ***************************************************************************/
GCTAModelIrfBackground::GCTAModelIrfBackground(const GCTAModelIrfBackground& bgd) :
                        GModelData(bgd)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(bgd);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelIrfBackground::~GCTAModelIrfBackground(void)
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
 * @param[in] bgd CTA instrument background model.
 * @return CTA instrument background model.
 ***************************************************************************/
GCTAModelIrfBackground& GCTAModelIrfBackground::operator=(const GCTAModelIrfBackground& bgd)
{
    // Execute only if object is not identical
    if (this != &bgd) {

        // Copy base class members
        this->GModelData::operator=(bgd);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(bgd);

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
 * @brief Clear CTA instrument background model
 *
 * This method properly resets the CTA instrument background model to an
 * initial state.
 ***************************************************************************/
void GCTAModelIrfBackground::clear(void)
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
 * @brief Clone CTA instrument background model
 *
 * @return Pointer to deep copy of CTA instrument background model.
 ***************************************************************************/
GCTAModelIrfBackground* GCTAModelIrfBackground::clone(void) const
{
    return new GCTAModelIrfBackground(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            Specified observation is not of the expected type.
 *
 * @todo Make sure that DETX and DETY are always set in GCTAInstDir.
 ***************************************************************************/
double GCTAModelIrfBackground::eval(const GEvent& event,
                                    const GObservation& obs) const
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

    // Retrieve pointer to CTA background
    const GCTABackground* bgd = rsp->background();
    if (bgd == NULL) {
        std::string msg = "Specified observation contains no background"
                          " information.\n" + obs.print();
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Extract CTA instrument direction from event
    const GCTAInstDir* dir  = dynamic_cast<const GCTAInstDir*>(&(event.dir()));
    if (dir == NULL) {
        std::string msg = "No CTA instrument direction found in event.";
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Set DETX and DETY in instrument direction
    GCTAInstDir inst_dir = cta->pointing().instdir(dir->dir());
    
    // Evaluate function
    double logE = event.energy().log10TeV();
    double spat = (*bgd)(logE, inst_dir.detx(), inst_dir.dety());
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time()) : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time()) : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // Apply deadtime correction
    value *= obs.deadc(event.time());

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            Specified observation is not of the expected type.
 *
 * @todo Make sure that DETX and DETY are always set in GCTAInstDir.
 ***************************************************************************/
double GCTAModelIrfBackground::eval_gradients(const GEvent& event,
                                              const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);
    if (cta == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Get pointer on CTA IRF response
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(cta->response());
    if (rsp == NULL) {
        std::string msg = "Specified observation does not contain an"
                          " IRF response.\n" + obs.print();
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Retrieve pointer to CTA background
    const GCTABackground* bgd = rsp->background();
    if (bgd == NULL) {
        std::string msg = "Specified observation contains no background"
                          " information.\n" + obs.print();
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Extract CTA instrument direction from event
    const GCTAInstDir* dir  = dynamic_cast<const GCTAInstDir*>(&(event.dir()));
    if (dir == NULL) {
        std::string msg = "No CTA instrument direction found in event.";
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Set DETX and DETY in instrument direction
    GCTAInstDir inst_dir = cta->pointing().instdir(dir->dir());

    // Evaluate function
    double logE = event.energy().log10TeV();
    double spat = (*bgd)(logE, inst_dir.detx(), inst_dir.dety());
    double spec = (spectral() != NULL)
                  ? spectral()->eval_gradients(event.energy(), event.time())
                  : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval_gradients(event.time())
                  : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // Apply deadtime correction
    double deadc = obs.deadc(event.time());
    value       *= deadc;

    // Multiply factors to spectral gradients
    if (spectral() != NULL) {
        double fact = spat * temp * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < spectral()->size(); ++i)
                (*spectral())[i].factor_gradient( (*spectral())[i].factor_gradient() * fact );
        }
    }

    // Multiply factors to temporal gradients
    if (temporal() != NULL) {
        double fact = spat * spec * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < temporal()->size(); ++i)
                (*temporal())[i].factor_gradient( (*temporal())[i].factor_gradient() * fact );
        }
    }

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
 * Spatially integrates the instrumental background model for a given
 * measured event energy and event time. This method also applies a deadtime
 * correction factor, so that the normalization of the model is a real rate
 * (counts/MeV/s).
 ***************************************************************************/
double GCTAModelIrfBackground::npred(const GEnergy&      obsEng,
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
				std::cout << "GCTAModelIrfBackground::npred:";
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
            
            // Retrieve pointer to CTA background
            const GCTABackground* bgd = rsp->background();
            if (bgd == NULL) {
                std::string msg = "Specified observation contains no background"
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

            // Get reference to ROI centre
            const GSkyDir& roi_centre = events->roi().centre().dir();

			// Get ROI radius in radians
			double roi_radius = events->roi().radius() * gammalib::deg2rad;

            // Get log10 of energy in TeV
            double logE = obsEng.log10TeV();

			// Setup integration function
			GCTAModelIrfBackground::npred_roi_kern_theta integrand(bgd, logE);

			// Setup integrator
			GIntegral integral(&integrand);
			integral.eps(g_cta_inst_background_npred_theta_eps);

			// Spatially integrate radial component
			npred = integral.romberg(0.0, roi_radius);

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
                std::string origin  = "GCTAModelIrfBackground::npred";
                std::string message = " NaN/Inf encountered (npred=" +
                                      gammalib::str(npred) + ", roi_radius=" +
                                      gammalib::str(roi_radius) + ")";
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
GCTAEventList* GCTAModelIrfBackground::mc(const GObservation& obs, GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Retrieve CTA observation
        const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);
        if (cta == NULL) {
            std::string msg = "Specified observation is not a CTA observation.\n" +
                              obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Get pointer on CTA IRF response
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(cta->response());
        if (rsp == NULL) {
            std::string msg = "Specified observation does not contain"
                              " an IRF response.\n" + obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Retrieve CTA response and pointing
        const GCTAPointing& pnt = cta->pointing();
        
        // Get pointer to CTA background
        const GCTABackground* bgd = rsp->background();
        if (bgd == NULL) {
            std::string msg = "Specified observation contains no background"
                              " information.\n" + obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Retrieve event list to access the ROI, energy boundaries and GTIs
        const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
        if (events == NULL) {
            std::string msg = "No CTA event list found in observation.\n" +
                              obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Get simulation region
        const GCTARoi&  roi     = events->roi();
        const GEbounds& ebounds = events->ebounds();
        const GGti&     gti     = events->gti();

        // Set simulation region for result event list
        list->roi(roi);
        list->ebounds(ebounds);
        list->gti(gti);

        // Create a spectral model that combines the information from the
        // background information and the spectrum provided by the model
        GModelSpectralNodes spectral(bgd->spectrum());
        for (int i = 0; i < spectral.nodes(); ++i) {
            GEnergy energy    = spectral.energy(i);
            double  intensity = spectral.intensity(i);
            double  norm      = m_spectral->eval(energy, events->tstart());
            spectral.intensity(i, norm*intensity);
        }

        // Loop over all energy boundaries
        for (int ieng = 0; ieng < ebounds.size(); ++ieng) {

            // Compute the background rate in model within the energy
            // boundaries from spectral component (units: cts/s).
            // Note that the time here is ontime. Deadtime correction will
            // be done later.
            double rate = spectral.flux(ebounds.emin(ieng), ebounds.emax(ieng));

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GCTAModelIrfBackground::mc(\"" << name() << "\": ";
            std::cout << "rate=" << rate << " cts/s)" << std::endl;
            #endif

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

                    // Get Monte Carlo event direction from spatial model.
                    // This only will set the DETX and DETY coordinates.
                    GCTAInstDir instdir = bgd->mc(energy, times[i], ran);

                    // Derive sky direction from instrument coordinates
                    GSkyDir skydir = pnt.skydir(instdir);

                    // Set sky direction in GCTAInstDir object
                    instdir.dir(skydir);

                    // Allocate event
                    GCTAEventAtom event;

                    // Set event attributes
                    event.dir(instdir);
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
 * @brief Read CTA instrument background model from XML element
 *
 * @param[in] xml XML element.
 *
 * Set up CTA instrument background model from the information provided by
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
void GCTAModelIrfBackground::read(const GXmlElement& xml)
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

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set observation identifiers
    ids(xml.attribute("id"));

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA instrument background model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write CTA instrument background model information into an XML element.
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
void GCTAModelIrfBackground::write(GXmlElement& xml) const
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

    // Set model type, name and optionally instruments
    src->attribute("name", name());
    src->attribute("type", type());
    if (instruments().length() > 0) {
        src->attribute("instrument", instruments());
    }
    std::string identifiers = ids();
    if (identifiers.length() > 0) {
        src->attribute("id", identifiers);
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA instrument background model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing CTA instrument background model information.
 ***************************************************************************/
std::string GCTAModelIrfBackground::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelIrfBackground ===");

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
void GCTAModelIrfBackground::init_members(void)
{
    // Initialise members
    m_spectral = NULL;
    m_temporal = NULL;

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
 * @param[in] bgd CTA background model.
 ***************************************************************************/
void GCTAModelIrfBackground::copy_members(const GCTAModelIrfBackground& bgd)
{
    // Copy cache
    m_npred_names    = bgd.m_npred_names;
    m_npred_energies = bgd.m_npred_energies;
    m_npred_times    = bgd.m_npred_times;
    m_npred_values   = bgd.m_npred_values;

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
void GCTAModelIrfBackground::free_members(void)
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
void GCTAModelIrfBackground::set_pointers(void)
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
bool GCTAModelIrfBackground::valid_model(void) const
{
    // Set result
    bool result = ((spectral() != NULL) && (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Construct spectral model from XML element
 *
 * @param[in] spectral XML element containing spectral model information.
 *
 * @exception GException::model_invalid_spectral
 *            Invalid spectral model type encountered.
 *
 * Returns pointer to a spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelIrfBackground::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model type
    std::string type = spectral.attribute("type");

    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(spectral);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct temporal model from XML element
 *
 * @param[in] temporal XML element containing temporal model information.
 *
 * @exception GException::model_invalid_temporal
 *            Invalid temporal model type encountered.
 *
 * Returns pointer to a temporal model that is defined in an XML element.
 ***************************************************************************/
GModelTemporal* GCTAModelIrfBackground::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model type
    std::string type = temporal.attribute("type");

    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(temporal);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_temporal(G_XML_TEMPORAL, type);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Kernel for offset angle integration of background model
 *
 * @param[in] theta Offset angle from ROI centre (radians).
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
double GCTAModelIrfBackground::npred_roi_kern_theta::eval(const double& theta)
{
    // Initialise value
    double value = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Setup phi integration kernel
        GCTAModelIrfBackground::npred_roi_kern_phi integrand(m_bgd,
                                                              m_logE,
                                                              theta);

        // Integrate over phi
        GIntegral integral(&integrand);
        integral.eps(g_cta_inst_background_npred_phi_eps);
        value = integral.romberg(0.0, gammalib::twopi) * std::sin(theta);

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::string origin  = "GCTAModelIrfBackground::npred_roi_kern_theta::eval"
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
 * @brief Kernel for azimuth angle integration of background model
 *
 * @param[in] phi Azimuth angle around ROI centre (radians).
 *
 * Computes
 *
 * \f[
 *    B(\theta, \phi | E, t)
 * \f]
 *
 * using
 *
 * \f[ {\rm detx} = \theta \cos \phi \f]
 * \f[ {\rm dety} = \theta \sin \phi \f]
 *
 * @todo Verify correct orientation of detx and dety with respect to phi
 ***************************************************************************/
double GCTAModelIrfBackground::npred_roi_kern_phi::eval(const double& phi)
{
	// Compute detx and dety
    double detx(0.0);
    double dety(0.0);
	if (m_theta > 0.0 ) {
		detx = m_theta * std::cos(phi);
		dety = m_theta * std::sin(phi);
	}

    // Get background value
    double value = (*m_bgd)(m_logE, detx, dety);

	// Debug: Check for NaN
	#if defined(G_NAN_CHECK)
	if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::string origin  = "GCTAModelIrfBackground::npred_roi_kern_phi::eval"
                              "(" + gammalib::str(phi) + ")";
        std::string message = " NaN/Inf encountered (value=" +
                              gammalib::str(value) + ", detx=" +
                              gammalib::str(detx) + ", dety=" +
                              gammalib::str(dety) + ")";
        gammalib::warning(origin, message);
	}
	#endif

    // Return Npred
    return value;
}
