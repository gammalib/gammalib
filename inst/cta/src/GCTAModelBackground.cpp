/***************************************************************************
 *           GCTAModelBackground.cpp - Background model class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelBackground.cpp
 * @brief Background model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GIntegral.hpp"
#include "GCTAModelSpatialRegistry.hpp"
#include "GCTAModelBackground.hpp"
#include "GCTAObservation.hpp"
#include "GCTAPointing.hpp"
#include "GCTAInstDir.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelBackground g_cta_background_seed;
const GModelRegistry      g_cta_background_registry(&g_cta_background_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL     "GCTAModelBackground::eval(GEvent&, GObservation&, bool&)"
#define G_NPRED "GCTAModelBackground::npred(GEnergy&, GTime&, GObservation&)"
#define G_MC                  "GCTAModelBackground::mc(GObservation&, GRan&)"
#define G_XML_SPATIAL        "GCTAModelBackground::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL      "GCTAModelBackground::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL      "GCTAModelBackground::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC                                  //!< Dump MC information


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty background model.
 ***************************************************************************/
GCTAModelBackground::GCTAModelBackground(void) : GModelData()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a background model from the information that is found in a XML
 * element. Please refer to the read() method to learn about the expected
 * structure of the XML element.
 ***************************************************************************/
GCTAModelBackground::GCTAModelBackground(const GXmlElement& xml) :
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
 * @brief Construct from spatial and spectral components
 *
 * @param[in] spatial Spatial model component.
 * @param[in] spectral Spectral model component.
 *
 * Constructs a background model from a @p spatial and a @p spectral model
 * component. The temporal component is assumed to be constant.
 ***************************************************************************/
GCTAModelBackground::GCTAModelBackground(const GCTAModelSpatial& spatial,
                                         const GModelSpectral&   spectral) :
                     GModelData()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_spatial  = spatial.clone();
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
 * @param[in] spatial Spatial model component.
 * @param[in] spectral Spectral model component.
 * @param[in] temporal Temporal model component.
 *
 * Constructs a background model from a @p spatial, a @p spectral and a
 * @p temporal model component.
 ***************************************************************************/
GCTAModelBackground::GCTAModelBackground(const GCTAModelSpatial& spatial,
                                         const GModelSpectral&   spectral,
                                         const GModelTemporal&   temporal) :
                     GModelData()
{
    // Initialise members
    init_members();

    // Clone model components
    m_spatial  = spatial.clone();
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
 * @param[in] model Background model.
 *
 * Constructs a background model by copying information from an existing
 * model. Note that the copy is a deep copy, so the original object can be
 * destroyed after the copy without any loss of information.
 ***************************************************************************/
GCTAModelBackground::GCTAModelBackground(const GCTAModelBackground& model) :
                     GModelData(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys a background model.
 ***************************************************************************/
GCTAModelBackground::~GCTAModelBackground(void)
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
 * @param[in] model Background model.
 *
 * Assigns the information from a background model to instance. Note that a
 * deep copy of the information is performed, so the original instance can be
 * destroyed after the assignment without any loss of information.
 ***************************************************************************/
GCTAModelBackground& GCTAModelBackground::operator=(const GCTAModelBackground& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelData::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members (this method also sets the parameter pointers)
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * Resets the instance to a clean initial state. All information that resided
 * in the object will be lost.
 ***************************************************************************/
void GCTAModelBackground::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelData::free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    this->GModelData::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of background model.
 ***************************************************************************/
GCTAModelBackground* GCTAModelBackground::clone(void) const
{
    return new GCTAModelBackground(*this);
}


/***********************************************************************//**
 * @brief Set spatial model component
 *
 * @param[in] spatial Pointer to spatial model component.
 *
 * Sets the spatial model component of the model.
 ***************************************************************************/
void GCTAModelBackground::spatial(const GCTAModelSpatial* spatial)
{
    // Free spatial model component only if it differs from current
    // component. This prevents unintentional deallocation of the argument
    if ((m_spatial != NULL) && (m_spatial != spatial)) {
        delete m_spatial;
    }

    // Clone spatial model component if it exists, otherwise set pointer
    // to NULL
    m_spatial = (spatial != NULL) ? spatial->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set spectral model component
 *
 * @param[in] spectral Pointer to spectral model component.
 *
 * Sets the spectral model component of the model.
 ***************************************************************************/
void GCTAModelBackground::spectral(const GModelSpectral* spectral)
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
void GCTAModelBackground::temporal(const GModelTemporal* temporal)
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
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Function value.
 *
 * Evaluates the background model which is a factorization of a spatial,
 * spectral and temporal model component. This method also applies a deadtime
 * correction factor, so that the normalization of the model is a real rate
 * (counts/exposure time).
 *
 * If the @p gradients flag is true the method will also set the parameter
 * gradients of the model parameters.
 *
 * @todo Add bookkeeping of last value and evaluate only if argument
 *       changed
 * @todo Verify that CTA instrument direction pointer is valid, or better,
 *       add an offset method to GCTAPointing. Ideally, we should precompute
 *       all offset angles (for an event cube this may only be done easily
 *       if the pointing has been fixed; otherwise we need a structure
 *       similar to the Fermi/LAT livetime cube that provides the effective
 *       sky exposure as function of offset angle).
 ***************************************************************************/
double GCTAModelBackground::eval(const GEvent&       event,
                                 const GObservation& obs,
                                 const bool&         gradients) const
{
    // Get reference on CTA instrument direction from event
    const GCTAInstDir& dir = gammalib::cta_dir(G_EVAL, event);

    // Evaluate function and gradients
    double spat = (spatial() != NULL)
                  ? spatial()->eval(dir, event.energy(), event.time(), gradients) : 1.0;
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

        // Multiply factors to radial gradients
        if (spatial() != NULL) {
            double fact = spec * temp * deadc;
            if (fact != 1.0) {
                for (int i = 0; i < spatial()->size(); ++i)
                    (*spatial())[i].factor_gradient((*spatial())[i].factor_gradient() * fact);
            }
        }

        // Multiply factors to spectral gradients
        if (spectral() != NULL) {
            double fact = spat * temp * deadc;
            if (fact != 1.0) {
                for (int i = 0; i < spectral()->size(); ++i)
                    (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact);
            }
        }

        // Multiply factors to temporal gradients
        if (temporal() != NULL) {
            double fact = spat * spec * deadc;
            if (fact != 1.0) {
                for (int i = 0; i < temporal()->size(); ++i)
                    (*temporal())[i].factor_gradient((*temporal())[i].factor_gradient() * fact);
            }
        }

    } // endif: computed partial derivatives

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated background model
 *
 * @param[in] energy Measured event energy.
 * @param[in] time Measured event time.
 * @param[in] obs Observation.
 * @return Integral of spatial model component.
 *
 * @exception GException::invalid_value
 *            Pointing direction differs from RoI centre.
 *
 * Spatially integrates the background model for a given measured event
 * energy and event time. This method also applies a deadtime correction
 * factor, so that the normalization of the model is a real rate
 * (counts/exposure time).
 ***************************************************************************/
double GCTAModelBackground::npred(const GEnergy&      energy,
                                  const GTime&        time,
                                  const GObservation& obs) const
{
    // Get reference on CTA pointing and event list from observation
    const GCTAPointing&  pnt    = gammalib::cta_pnt(G_NPRED, obs);
    const GCTAEventList& events = gammalib::cta_event_list(G_NPRED, obs);

    // Get reference to pointing direction and RoI centre
    const GSkyDir& pointing   = pnt.dir();
    const GSkyDir& roi_centre = events.roi().centre().dir();

    // Throw an exception if both differ significantly
    if (pointing.dist(roi_centre) > 1.0e-4) {
        std::string msg = "Pointing direction ("+
                          gammalib::str(pointing.ra_deg())+","+
                          gammalib::str(pointing.dec_deg())+") differs "
                          "significantly from RoI centre ("+
                          gammalib::str(roi_centre.ra_deg())+","+
                          gammalib::str(roi_centre.dec_deg())+"). "
                          "Method is only valid for RoI centres that are "
                          "identical to the pointing direction.";
        throw GException::invalid_value(G_NPRED, msg);
    }

    // Get spatially integrated model component
    double npred = spatial()->npred(energy, time, obs);

    // Multiply in spectral and temporal components
    npred *= spectral()->eval(energy, time);
    npred *= temporal()->eval(time);

    // Apply deadtime correction
    npred *= obs.deadc(time);

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
 * The method also applies a deadtime correction using a Monte Carlo process,
 * taking into account temporal deadtime variations. For this purpose, the
 * method makes use of the time dependent GObservation::deadc method.
 *
 * For each event in the returned event list, the sky direction, the nominal
 * coordinates (DETX and DETY), the energy and the time will be set.
 *
 * @todo Handle special case of cube spatial model
 ***************************************************************************/
GCTAEventList* GCTAModelBackground::mc(const GObservation& obs, 
                                       GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Get reference on CTA pointing, background response and event list
        // from observation
        const GCTAPointing&  pnt    = gammalib::cta_pnt(G_MC, obs);
        const GCTAEventList& events = gammalib::cta_event_list(G_MC, obs);

        // Get simulation region
        const GCTARoi&  roi     = events.roi();
        const GEbounds& ebounds = events.ebounds();
        const GGti&     gti     = events.gti();

        // Set simulation region for result event list
        list->roi(roi);
        list->ebounds(ebounds);
        list->gti(gti);

        // Set spectral model
        const GModelSpectral* spectral = m_spectral;

        // TODO: Handle special case of cube spatial model. This will
        // replace the spectral model by a spectral nodes model.

        // Get solid angle of spatial model. This only works for an energy
        // independent spatial model!!!
        double solidangle = m_spatial->npred(GEnergy(), GTime(), obs);

        // Loop over all energy boundaries
        for (int ieng = 0; ieng < ebounds.size(); ++ieng) {

            // Compute the background rate in model within the energy
            // boundaries from spectral component (units: cts/s/sr).
            double flux = spectral->flux(ebounds.emin(ieng), ebounds.emax(ieng));

            // Derive expecting rate (units: cts/s). Note that the time here
            // is good time. Deadtime correction will be done later.
            double rate = flux * solidangle;

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GCTAModelBackground::mc(\"" << name() << "\": ";
            std::cout << "flux=" << flux << " cts/s/sr, ";
            std::cout << "solidangle=" << solidangle << " sr, ";
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
                std::cout << " events=" << n_events << std::endl;
                int n_removed_by_deadtime = 0;
                int n_trial_outside_roi   = 0;
                #endif

                // Loop over events
                for (int i = 0; i < n_events; ++i) {

                    // Apply deadtime correction
                    double deadc = obs.deadc(times[i]);
                    if (deadc < 1.0) {
                        if (ran.uniform() > deadc) {
                            #if defined(G_DUMP_MC)
                            n_removed_by_deadtime++;
                            #endif
                            continue;
                        }
                    }

                    // Get Monte Carlo event energy from spectral model
                    GEnergy energy = spectral->mc(ebounds.emin(ieng),
                                                  ebounds.emax(ieng),
                                                  times[i],
                                                  ran);

                    // Get an instrument direction within the RoI. This is
                    // potentially tried 100 times so that if we really can't
                    // a valid instrument direction the code is not locked up
                    for (int k = 0; k < 100; ++k) {

                        // Get Monte Carlo event direction from spatial model.
                        // This only will set the DETX and DETY coordinates.
                        GCTAInstDir instdir = m_spatial->mc(energy, times[i], ran);

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

                        // Append event to list if it falls in RoI and break
                        // the look
                        if (events.roi().contains(event)) {
                            list->append(event);
                            break;
                        }

                        // ... otherwise optionally bookkeep the trial
                        #if defined(G_DUMP_MC)
                        else {
                            n_trial_outside_roi++;
                        }
                        #endif

                    } // endfor: trial look for instrument direction within RoI

                } // endfor: looped over all events

                // Debug option: provide statisics
                #if defined(G_DUMP_MC)
                std::cout << " Events removed due to deadtime=";
                std::cout << n_removed_by_deadtime << std::endl;
                std::cout << " Trials outside RoI=";
                std::cout << n_trial_outside_roi << std::endl;
                #endif

            } // endfor: looped over all GTIs

        } // endfor: looped over all energy boundaries

    } // endif: model was valid

    // Return
    return list;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the sky model from an XML element. The XML element is expected to
 * respect the following format:
 *
 *     <source name=".." type=".." instrument=".." id="..">
 *       <spectrum type="..">
 *         ..
 *       </spectrum>
 *       <spatialModel type="..">
 *         ..
 *       </spatialModel>
 *       <temporal type="..">
 *         ..
 *       </temporal>
 *     </source>
 *
 * The temporal element is optional. In no temporal element is specified a
 * constant component with unity normalization will be assumed.
 ***************************************************************************/
void GCTAModelBackground::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Get pointers on spatial and spectral model components
    const GXmlElement* spat = xml.element("spatialModel", 0);
    const GXmlElement* spec = xml.element("spectrum", 0);

    // Set spatial and spectral model components
    m_spatial  = xml_spatial(*spat);
    m_spectral = xml_spectral(*spec);

    // Handle optional temporal model
    if (xml.elements("temporal") > 0) {
        const GXmlElement* temp = xml.element("temporal", 0);
        m_temporal = xml_temporal(*temp);
    }
    else {
        GModelTemporalConst temporal;
        m_temporal = temporal.clone();
    }

    // Read model attributes
    read_attributes(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes the sky model into an XML source library. The format written to
 * the @p xml element is as follows:
 *
 *     <source name=".." type=".." instrument=".." id="..">
 *       <spectrum type="..">
 *         ..
 *       </spectrum>
 *       <spatialModel type="..">
 *         ..
 *       </spatialModel>
 *       <temporal type="..">
 *         ..
 *       </temporal>
 *     </source>
 *
 * For compatibility reasons the temporal element will only be written if it
 * is a non-constant component or a constant component with a normalization
 * that differs from unity.
 ***************************************************************************/
void GCTAModelBackground::write(GXmlElement& xml) const
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

    // If the temporal model is not a constant with unit normalization then
    // set cons to a NULL pointer
    GModelTemporalConst* cons = dynamic_cast<GModelTemporalConst*>(temporal());
    if (cons != NULL) {
        if (cons->norm() != 1.0) {
            cons = NULL;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spatial()  != NULL) src->append(GXmlElement("spatialModel"));
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (temporal() != NULL && cons == NULL) src->append(GXmlElement("temporal"));
    }

    // Write spatial model
    if (spatial() != NULL) {
        GXmlElement* spat = src->element("spatialModel", 0);
        spatial()->write(*spat);
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write temporal model (only if not a constant with unit normalization
    // factor)
    if (temporal() != NULL && cons == NULL) {
        GXmlElement* temp = src->element("temporal", 0);
        temporal()->write(*temp);
    }

    // Write model attributes
    write_attributes(*src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GCTAModelBackground::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelBackground ===");

        // Determine number of parameters per type
        int n_spatial  = (spatial()  != NULL) ? spatial()->size()  : 0;
        int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
        int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type")+type());

        // Append model type
        result.append("\n"+gammalib::parformat("Model components"));
        if (n_spatial > 0) {
            result.append("\""+spatial()->type()+"\"");
            if (n_spectral > 0 || n_temporal > 0) {
                result.append(" * ");
            }
        }
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
        result.append("\n"+gammalib::parformat("Number of spatial par's") +
                      gammalib::str(n_spatial));
        for (int i = 0; i < n_spatial; ++i) {
            result.append("\n"+(*spatial())[i].print());
        }
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
void GCTAModelBackground::init_members(void)
{
    // Initialise members
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Background model.
 ***************************************************************************/
void GCTAModelBackground::copy_members(const GCTAModelBackground& model)
{
    // Clone radial, spectral and temporal model components
    m_spatial  = (model.m_spatial  != NULL) ? model.m_spatial->clone()  : NULL;
    m_spectral = (model.m_spectral != NULL) ? model.m_spectral->clone() : NULL;
    m_temporal = (model.m_temporal != NULL) ? model.m_temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelBackground::free_members(void)
{
    // Free memory
    if (m_spatial  != NULL) delete m_spatial;
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_spatial  = NULL;
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
void GCTAModelBackground::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Determine the number of parameters
    int n_spatial  = (spatial()  != NULL) ? spatial()->size()  : 0;
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spatial + n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather spatial parameter pointers
        for (int i = 0; i < n_spatial; ++i) {
            m_pars.push_back(&((*spatial())[i]));
        }

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
 * @brief Return pointer to spatial model from XML element
 *
 * @param[in] spatial XML element.
 * @return Pointer to spatial model.
 *
 * Returns pointer to a spatial model that is defined in an XML element.
 ***************************************************************************/
GCTAModelSpatial* GCTAModelBackground::xml_spatial(const GXmlElement& spatial) const
{
    // Get radial model
    GCTAModelSpatialRegistry registry;
    GCTAModelSpatial*        ptr = registry.alloc(spatial);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return pointer to spectral model from XML element
 *
 * @param[in] spectral XML element.
 * @return Pointer to spectral model.
 *
 * Returns pointer to spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelBackground::xml_spectral(const GXmlElement& spectral) const
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
GModelTemporal* GCTAModelBackground::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(temporal);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 *
 * Returns 'true' if models has a spatial, a spectral and a temporal
 * component. Otherwise returns 'false'.
 ***************************************************************************/
bool GCTAModelBackground::valid_model(void) const
{
    // Set result
    bool result = ((spatial()  != NULL) &&
                   (spectral() != NULL) &&
                   (temporal() != NULL));

    // Return result
    return result;
}
