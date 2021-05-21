/***************************************************************************
 *      GCTAModelRadialAcceptance.cpp - Radial acceptance model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialAcceptance.cpp
 * @brief Radial acceptance model class implementation
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
#include "GCTAModelRadialRegistry.hpp"
#include "GCTAModelRadialAcceptance.hpp"
#include "GCTAObservation.hpp"
#include "GCTAPointing.hpp"
#include "GCTAInstDir.hpp"
#include "GCTARoi.hpp"
#include "GCTASupport.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRadialAcceptance g_cta_radial_acceptance_seed;
const GModelRegistry            g_cta_radial_acceptance_registry(&g_cta_radial_acceptance_seed);

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                "GCTAModelRadialAcceptance::operator() (int)"
#define G_EVAL     "GCTAModelRadialAcceptance::eval(GEvent&, GObservation&, "\
                                                                     "bool&)"
#define G_NPRED          "GCTAModelRadialAcceptance::npred(GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC            "GCTAModelRadialAcceptance::mc(GObservation&, GRan&)"
#define G_XML_RADIAL    "GCTAModelRadialAcceptance::xml_radial(GXmlElement&)"
#define G_XML_SPECTRAL "GCTAModelRadialAcceptance::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL "GCTAModelRadialAcceptance::xml_temporal(GXmlElement&)"

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
 * Constructs an empty CTA radial acceptance model.
 ***************************************************************************/
GCTAModelRadialAcceptance::GCTAModelRadialAcceptance(void) : GModelData()
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
 * Constructs a CTA radial acceptance model from the information that is
 * found in a XML element. Please refer to the read() method to learn about
 * the expected structure of the XML element.
 ***************************************************************************/
GCTAModelRadialAcceptance::GCTAModelRadialAcceptance(const GXmlElement& xml) :
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
 * @brief Construct from radial and spectral components
 *
 * @param[in] radial Radial CTA model component.
 * @param[in] spectral Spectral model component.
 *
 * Constructs a CTA radial acceptance model from a @p radial and a
 * @p spectral model component. The temporal component is assumed to be
 * constant.
 ***************************************************************************/
GCTAModelRadialAcceptance::GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                                                     const GModelSpectral&  spectral) :
                           GModelData()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_radial   = radial.clone();
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
 * @param[in] radial Radial CTA model component.
 * @param[in] spectral Spectral model component.
 * @param[in] temporal Temporal model component.
 *
 * Constructs a CTA radial acceptance model from a @p radial, a @p spectral
 * and a @p temporal model component.
 ***************************************************************************/
GCTAModelRadialAcceptance::GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                                                     const GModelSpectral&  spectral,
                                                     const GModelTemporal&  temporal) :
                           GModelData()
{
    // Initialise members
    init_members();

    // Clone model components
    m_radial   = radial.clone();
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
 * @param[in] model Radial acceptance model.
 *
 * Constructs a CTA radial acceptance model by copying information from an
 * existing model. Note that the copy is a deep copy, so the original object
 * can be destroyed after the copy without any loss of information.
 ***************************************************************************/
GCTAModelRadialAcceptance::GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model) :
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
 * Destroys a CTA radial acceptance model.
 ***************************************************************************/
GCTAModelRadialAcceptance::~GCTAModelRadialAcceptance(void)
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
 * @param[in] model Radial acceptance model.
 *
 * Assigns the information from a CTA radial acceptance model to the actual
 * object. Note that a deep copy of the information is performed, so the
 * original object can be destroyed after the assignment without any loss of
 * information.
 ***************************************************************************/
GCTAModelRadialAcceptance& GCTAModelRadialAcceptance::operator=(const GCTAModelRadialAcceptance& model)
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
 * Resets the object to a clean initial state. All information that resided
 * in the object will be lost.
 ***************************************************************************/
void GCTAModelRadialAcceptance::clear(void)
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
 * @return Pointer to deep copy of CTA radial acceptance model.
 ***************************************************************************/
GCTAModelRadialAcceptance* GCTAModelRadialAcceptance::clone(void) const
{
    return new GCTAModelRadialAcceptance(*this);
}


/***********************************************************************//**
 * @brief Set radial model component
 *
 * @param[in] radial Pointer to radial model component.
 *
 * Sets the radial model component of the model.
 ***************************************************************************/
void GCTAModelRadialAcceptance::radial(const GCTAModelRadial* radial)
{
    // Free spatial model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_radial != NULL) && (m_radial != radial)) {
        delete m_radial;
    }

    // Clone spatial model component if it exists, otherwise set pointer
    // to NULL
    m_radial = (radial != NULL) ? radial->clone() : NULL;

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
void GCTAModelRadialAcceptance::spectral(const GModelSpectral* spectral)
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
void GCTAModelRadialAcceptance::temporal(const GModelTemporal* temporal)
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
 * Evaluates tha CTA radial acceptance model which is a factorization of a
 * spatial, spectral and temporal model component. The method returns a
 * real rate, defined by the number of counts per MeV, steradian and ontime.
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
double GCTAModelRadialAcceptance::eval(const GEvent&       event,
                                       const GObservation& obs,
                                       const bool&         gradients) const
{
    // Get reference on CTA pointing from observation and reference on CTA
    // instrument direction from event
    const GCTAPointing& pnt = gammalib::cta_pnt(G_EVAL, obs);
    const GCTAInstDir&  dir = gammalib::cta_dir(G_EVAL, event);

    // Compute offset angle (in degrees)
    double offset = dir.dir().dist_deg(pnt.dir());

    // Evaluate function and gradients
    double rad  = (radial()   != NULL)
                  ? radial()->eval(offset, gradients) : 1.0;
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time(), gradients)
                  : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time(), gradients) : 1.0;

    // Compute value
    double value = rad * spec * temp;

    // Optionally compute partial derivatives
    if (gradients) {

        // Multiply factors to radial gradients
        if (radial() != NULL) {
            double fact = spec * temp;
            if (fact != 1.0) {
                for (int i = 0; i < radial()->size(); ++i)
                    (*radial())[i].factor_gradient((*radial())[i].factor_gradient() * fact);
            }
        }

        // Multiply factors to spectral gradients
        if (spectral() != NULL) {
            double fact = rad * temp;
            if (fact != 1.0) {
                for (int i = 0; i < spectral()->size(); ++i)
                    (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact);
            }
        }

        // Multiply factors to temporal gradients
        if (temporal() != NULL) {
            double fact = rad * spec;
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
 * @brief Return spatially integrated background rate in units of
 *        events MeV\f$^{-1}\f$ s\f$^{-1}\f$
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated background rate
 *         (events MeV\f$^{-1}\f$ s\f$^{-1}\f$)
 *
 * Spatially integrates the data model for a given measured event energy and
 * event time. The method returns a real rate, defined as the number of
 * counts per MeV and ontime.
 ***************************************************************************/
double GCTAModelRadialAcceptance::npred(const GEnergy&      obsEng,
                                        const GTime&        obsTime,
                                        const GObservation& obs) const
{
    // Initialise result
    double npred = 0.0;

    // Evaluate only if model is valid
    if (valid_model()) {

        // Retrieve CTA pointing and event list
        const GCTAPointing&  pnt    = gammalib::cta_pnt(G_NPRED, obs);
        const GCTAEventList& events = gammalib::cta_event_list(G_NPRED, obs);

        // Get ROI radius in radians
        double roi_radius = events.roi().radius() * gammalib::deg2rad;

        // Get distance from ROI centre in radians
        double roi_distance = events.roi().centre().dir().dist(pnt.dir());

        // Setup integration function
        GCTAModelRadialAcceptance::roi_kern integrand(radial(), roi_radius, roi_distance);

        // Setup integrator
        GIntegral integral(&integrand);

        // Setup integration boundaries
        double rmin = (roi_distance > roi_radius) ? roi_distance-roi_radius : 0.0;
        double rmax = roi_radius + roi_distance;

        // Spatially integrate radial component
        npred = integral.romberg(rmin, rmax);

        // Multiply in spectral and temporal components
        npred *= spectral()->eval(obsEng, obsTime);
        npred *= temporal()->eval(obsTime);

    } // endif: model was valid

    // Return
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 *
 * Draws a sample of events from the radial acceptance model using a Monte
 * Carlo simulation. The pointing information, the energy boundaries and the
 * good time interval for the sampling will be extracted from the observation
 * argument that is passed to the method. The method also requires a random
 * number generator of type GRan which is passed by reference, hence the
 * state of the random number generator will be changed by the method.
 ***************************************************************************/
GCTAEventList* GCTAModelRadialAcceptance::mc(const GObservation& obs, 
                                             GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Retrieve CTA pointing
        const GCTAPointing& pnt = gammalib::cta_pnt(G_MC, obs);

        // Convert CTA pointing direction in instrument system
        GCTAInstDir pnt_dir(pnt.dir());

        // Loop over all energy boundaries
        for (int ieng = 0; ieng < obs.events()->ebounds().size(); ++ieng) {

            // Compute the on-axis background rate in model within the
            // energy boundaries from spectral component (units: cts/s/sr)
            double flux = spectral()->flux(obs.events()->ebounds().emin(ieng),
                                           obs.events()->ebounds().emax(ieng));

            // Compute solid angle used for normalization
            double area = radial()->omega();

            // Derive expecting rate (units: cts/s). Note that the time here
            // is good time. Deadtime correction will be done later.
            double rate = flux * area;

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GCTAModelRadialAcceptance::mc(\"" << name() << "\": ";
            std::cout << "flux=" << flux << " cts/s/sr, ";
            std::cout << "area=" << area << " sr, ";
            std::cout << "rate=" << rate << " cts/s)" << std::endl;
            #endif

            // Loop over all good time intervals
            for (int itime = 0; itime < obs.events()->gti().size(); ++itime) {

                // Get event arrival times from temporal model
                GTimes times = m_temporal->mc(rate,
                                              obs.events()->gti().tstart(itime),
                                              obs.events()->gti().tstop(itime),
                                              ran);

                // Get number of events
                int n_events = times.size();

                // Reserve space for events
                if (n_events > 0) {
                    list->reserve(n_events);
                }

                // Loop over events
                for (int i = 0; i < n_events; ++i) {

                    // Set event energy
                    GEnergy energy = spectral()->mc(obs.events()->ebounds().emin(ieng),
                                                    obs.events()->ebounds().emax(ieng),
                                                    times[i],
                                                    ran);

                    // Get Monte Carlo event direction from radial model.
                    // This only will set the DETX and DETY coordinates.
                    GCTAInstDir instdir = radial()->mc(ran);

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

                    // Append event to list
                    list->append(event);

                } // endfor: looped over all events

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
 * The model is composed of a spectrum component ('spectral'), a radial
 * component ('radialModel'), and, optionally, of a temporal component
 * ('temporal'). If no temporal component is found a constant model is
 * assumed.
 ***************************************************************************/
void GCTAModelRadialAcceptance::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Initialise XML elements
    const GXmlElement* rad  = NULL;
    const GXmlElement* spec = NULL;

    // Get pointers on spectrum and radial model
    rad  = xml.element("radialModel", 0);
    spec = xml.element("spectrum", 0);

    // Clone radial and spectral models
    m_radial   = xml_radial(*rad);
    m_spectral = xml_spectral(*spec);

    // Optionally get temporal model
    try {
        const GXmlElement* temp = xml.element("temporal", 0);
        m_temporal = xml_temporal(*temp);
    }
    catch (GException::xml_name_not_found &e) {
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
 * @todo Document method.
 ***************************************************************************/
void GCTAModelRadialAcceptance::write(GXmlElement& xml) const
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
        if (radial()   != NULL) src->append(GXmlElement("radialModel"));
        if (write_temporal)     src->append(GXmlElement("temporal"));
    }

    // Write radial model
    if (radial()) {
        GXmlElement* rad = src->element("radialModel", 0);
        radial()->write(*rad);
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Optionally write temporal model
    if (write_temporal) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("temporal", 0);
            temporal()->write(*temp);
        }
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
std::string GCTAModelRadialAcceptance::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelRadialAcceptance ===");

        // Determine number of parameters per type
        int n_radial   = (radial()   != NULL) ? radial()->size()   : 0;
        int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
        int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type"));
        if (n_radial > 0) {
            result.append("\""+radial()->type()+"\"");
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
        result.append("\n"+gammalib::parformat("Number of radial par's") +
                      gammalib::str(n_radial));
        for (int i = 0; i < n_radial; ++i) {
            result.append("\n"+(*radial())[i].print());
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
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelRadialAcceptance::init_members(void)
{
    // Initialise members
    m_radial   = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelRadialAcceptance::copy_members(const GCTAModelRadialAcceptance& model)
{
    // Clone radial, spectral and temporal model components
    m_radial   = (model.m_radial   != NULL) ? model.m_radial->clone()   : NULL;
    m_spectral = (model.m_spectral != NULL) ? model.m_spectral->clone() : NULL;
    m_temporal = (model.m_temporal != NULL) ? model.m_temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelRadialAcceptance::free_members(void)
{
    // Free memory
    if (m_radial   != NULL) delete m_radial;
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_radial   = NULL;
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
void GCTAModelRadialAcceptance::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Determine the number of parameters
    int n_radial   = (radial()   != NULL) ? radial()->size()   : 0;
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_radial + n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather radial parameter pointers
        for (int i = 0; i < n_radial; ++i) {
            m_pars.push_back(&((*radial())[i]));
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
 * @brief Verifies if model has all components
 *
 * Returns 'true' if models has a radial, a spectral and a temporal
 * component. Otherwise returns 'false'.
 ***************************************************************************/
bool GCTAModelRadialAcceptance::valid_model(void) const
{
    // Set result
    bool result = ((radial()   != NULL) &&
                   (spectral() != NULL) &&
                   (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Construct radial model from XML element
 *
 * @param[in] radial XML element containing radial model information.
 *
 * @exception GException::invalid_argument
 *            Invalid radial model type encountered.
 *
 * Returns pointer to a radial model that is defined in an XML element.
 ***************************************************************************/
GCTAModelRadial* GCTAModelRadialAcceptance::xml_radial(const GXmlElement& radial) const
{
    // Get radial model type
    std::string type = radial.attribute("type");

    // Get radial model
    GCTAModelRadialRegistry registry;
    GCTAModelRadial*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(radial);
    }

    // ... otherwise throw an exception
    else {
        std::string msg = "Invalid radial model type \""+type+"\" encountered. "
                          "No such model exists in the registry of radial "
                          "models. Please specify valid radial model type.";
        throw GException::invalid_argument(G_XML_RADIAL, msg);
    }

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
GModelSpectral* GCTAModelRadialAcceptance::xml_spectral(const GXmlElement& spectral) const
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
GModelTemporal* GCTAModelRadialAcceptance::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(temporal);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Integration kernel for the Npred method
 *
 * @param[in] offset Offset angle (radians).
 *
 * Computes
 * \f[N_{\rm pred} = \int_{\rm ROI} r(\theta) \sin \theta {\rm d}\phi\f]
 * where
 * \f$r(\theta)\f$ is the radial acceptance function,
 * \f$\theta\f$ is the measured offset angle from the pointing direction
 * in radians, and
 * \f$\phi\f$ is the measured azimuth angle. The integration is done over
 * the arc of the azimuth angle that lies within the ROI. This integration
 * is done analytically using the "roi_arclength" support function.
 ***************************************************************************/
double GCTAModelRadialAcceptance::roi_kern::eval(const double& offset)
{
    // Circumvent const correctness
    GCTAModelRadial* radial = const_cast<GCTAModelRadial*>(m_parent);

    // Initialise Npred value
    double value = 0.0;

    // Continue only if offset > 0
    if (offset > 0.0) {

        // Get arclength for given radius in radians.
        double phi = gammalib::roi_arclength(offset,
                                             m_dist,
                                             m_cosdist,
                                             m_sindist,
                                             m_roi,
                                             m_cosroi);

        // Get kernel value if phi > 0
        if (phi > 0.0) {
            value = radial->eval(offset*gammalib::rad2deg) * phi *
                    std::sin(offset);
        }

    } // endif: offset was positive

    // Return
    return value;
}
