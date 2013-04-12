/***************************************************************************
 *      GCTAModelRadialAcceptance.cpp - Radial acceptance model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
#include "GCTAException.hpp"
#include "GCTASupport.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRadialAcceptance g_cta_radial_acceptance_seed;
const GModelRegistry            g_cta_radial_acceptance_registry(&g_cta_radial_acceptance_seed);

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                "GCTAModelRadialAcceptance::operator() (int)"
#define G_EVAL                     "GCTAModelRadialAcceptance::eval(GEvent&,"\
                                                            " GObservation&)"
#define G_EVAL_GRADIENTS "GCTAModelRadialAcceptance::eval_gradients(GEvent&,"\
                                                            " GObservation&)"
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
 * @brief Constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA radial acceptance model from the information that is
 * found in a XML element. Please refer to the method
 * GCTAModelRadialAcceptance::read
 * to learn more about the information that is expected in the XML element.
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
 * Constructs a CTA radial acceptance model from a radial and a spectral
 * model component. The temporal component is assumed to be constant.
 * Please refer to the classes GCTAModelRadial and GModelSpectral to learn
 * more about the definition of the radial and spectral components.
 ***************************************************************************/
GCTAModelRadialAcceptance::GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                                                     const GModelSpectral&  spectral)
                                                               : GModelData()
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
GCTAModelRadialAcceptance& GCTAModelRadialAcceptance::operator= (const GCTAModelRadialAcceptance& model)
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
 * Clone a CTA radial acceptance model. Cloning performs a deep copy of the
 * information, so the original object can be destroyed after cloning without
 * any loss of information.
 ***************************************************************************/
GCTAModelRadialAcceptance* GCTAModelRadialAcceptance::clone(void) const
{
    return new GCTAModelRadialAcceptance(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * Evaluates tha CTA radial acceptance model which is a factorization of a
 * spatial, spectral and temporal model component. This method also applies
 * a deadtime correction factor, so that the normalization of the model is
 * a real rate (counts/exposure time).
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
double GCTAModelRadialAcceptance::eval(const GEvent& event,
                                       const GObservation& obs) const
{
    // Extract CTA pointing direction
    GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing());
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_EVAL);
    }

    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Compute offset angle (in degrees)
    double offset = cta_dir->dist_deg(pnt->dir());

    // Evaluate function and gradients
    double rad  = (radial()   != NULL)
                  ? radial()->eval(offset) : 1.0;
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time()) : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time()) : 1.0;

    // Compute value
    double value = rad * spec * temp;

    // Apply deadtime correction
    value *= obs.deadc(event.time());

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation (not used).
 *
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found in observation
 *
 * Evaluates tha CTA radial acceptance model and parameter gradients. The CTA
 * radial acceptance model is a factorization of a spatial, spectral and
 * temporal model component. This method also applies a deadtime correction
 * factor, so that the normalization of the model is a real rate
 * (counts/exposure time).
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
double GCTAModelRadialAcceptance::eval_gradients(const GEvent& event,
                                                 const GObservation& obs) const
{
    // Extract CTA pointing direction
    GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing());
    if (pnt == NULL) {
        throw GCTAException::no_pointing(G_EVAL_GRADIENTS);
    }

    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Compute offset angle (in degrees)
    double offset = cta_dir->dist_deg(pnt->dir());

    // Evaluate function and gradients
    double rad  = (radial()   != NULL)
                  ? radial()->eval_gradients(offset) : 1.0;
    double spec = (spectral() != NULL)
                  ? spectral()->eval_gradients(event.energy(), event.time()) : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval_gradients(event.time()) : 1.0;

    // Compute value
    double value = rad * spec * temp;

    // Apply deadtime correction
    double deadc = obs.deadc(event.time());
    value       *= deadc;

    // Multiply factors to radial gradients
    if (radial() != NULL) {
        double fact = spec * temp * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < radial()->size(); ++i)
                (*radial())[i].factor_gradient( (*radial())[i].factor_gradient() * fact );
        }
    }

    // Multiply factors to spectral gradients
    if (spectral() != NULL) {
        double fact = rad * temp * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < spectral()->size(); ++i)
                (*spectral())[i].factor_gradient( (*spectral())[i].factor_gradient() * fact );
        }
    }

    // Multiply factors to temporal gradients
    if (temporal() != NULL) {
        double fact = rad * spec * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < temporal()->size(); ++i)
                (*temporal())[i].factor_gradient( (*temporal())[i].factor_gradient() * fact );
        }
    }

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated data model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 *
 * @exception GException::no_list
 *            No valid CTA event list found in observation
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found in observation
 *
 * Spatially integrates the data model for a given measured event energy and
 * event time. This method also applies a deadtime correction factor, so that
 * the normalization of the model is a real rate (counts/exposure time).
 ***************************************************************************/
double GCTAModelRadialAcceptance::npred(const GEnergy&      obsEng,
                                        const GTime&        obsTime,
                                        const GObservation& obs) const
{
    // Initialise result
    double npred = 0.0;

    // Evaluate only if model is valid
    if (valid_model()) {

        // Get pointer on CTA events list
        const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
        if (events == NULL) {
            throw GException::no_list(G_NPRED);
        }

        // Get CTA pointing direction
        GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing());
        if (pnt == NULL) {
            throw GCTAException::no_pointing(G_NPRED);
        }

        // Get ROI radius in radians
        double roi_radius = events->roi().radius() * gammalib::deg2rad;

        // Get distance from ROI centre in radians
        double roi_distance = events->roi().centre().dist(pnt->dir());

        // Setup integration function
        GCTAModelRadialAcceptance::roi_kern integrand(radial(), roi_radius, roi_distance);

        // Setup integrator
        GIntegral integral(&integrand);

        // Setup integration boundaries
        double rmin = (roi_distance > roi_radius) ? roi_distance-roi_radius : 0.0;
        double rmax = roi_radius + roi_distance;

        // Spatially integrate radial component
        npred = integral.romb(rmin, rmax);

        // Multiply in spectral and temporal components
        npred *= spectral()->eval(obsEng, obsTime);
        npred *= temporal()->eval(obsTime);

        // Apply deadtime correction
        npred *= obs.deadc(obsTime);

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
 * @exception GCTAException::no_pointing
 *            No CTA pointing found in observation.
 *
 * Draws a sample of events from the radial acceptance model using a Monte
 * Carlo simulation. The pointing information, the energy boundaries and the
 * good time interval for the sampling will be extracted from the observation
 * argument that is passed to the method. The method also requires a random
 * number generator of type GRan which is passed by reference, hence the
 * state of the random number generator will be changed by the method.
 *
 * The method also applies a deadtime correction using a Monte Carlo process,
 * taking into account temporal deadtime variations. For this purpose, the
 * method makes use of the time dependent GObservation::deadc method.
 ***************************************************************************/
GCTAEventList* GCTAModelRadialAcceptance::mc(const GObservation& obs, 
                                             GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Extract CTA pointing direction at beginning of observation
        GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing());
        if (pnt == NULL) {
            throw GCTAException::no_pointing(G_MC);
        }

        // Convert CTA pointing direction in instrument system
        GCTAInstDir pnt_dir(pnt->dir());

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

                    // Apply deadtime correction
                    double deadc = obs.deadc(times[i]);
                    if (deadc < 1.0) {
                        if (ran.uniform() > deadc) {
                            continue;
                        }
                    }

                    // Set event direction
                    GCTAInstDir dir = radial()->mc(pnt_dir, ran);

                    // Set event energy
                    GEnergy energy = spectral()->mc(obs.events()->ebounds().emin(ieng),
                                                    obs.events()->ebounds().emax(ieng),
                                                    times[i],
                                                    ran);

                    // Allocate event
                    GCTAEventAtom event;

                    // Set event attributes
                    event.dir(dir);
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
 * ('lightcurve'). If no temporal component is found a constant model is
 * assumed.
 ***************************************************************************/
void GCTAModelRadialAcceptance::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Initialise XML elements
    const GXmlElement* rad  = NULL;
    const GXmlElement* spec = NULL;
    const GXmlElement* temp = NULL;

    // Get pointers on spectrum and radial model
    rad  = xml.element("radialModel", 0);
    spec = xml.element("spectrum", 0);

    // Clone radial and spectral models
    m_radial   = xml_radial(*rad);
    m_spectral = xml_spectral(*spec);

    // Optionally get temporal model
    try {
        temp = xml.element("lightcurve", 0);
        m_temporal = xml_temporal(*temp);
    }
    catch (GException::xml_name_not_found &e) {
        GModelTemporalConst temporal;
        m_temporal = temporal.clone();
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

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (radial()   != NULL) src->append(GXmlElement("radialModel"));
    }

    // Set model type, name and optionally instruments
    src->attribute("name", name());
    src->attribute("type", type());
    if (instruments().length() > 0) {
        src->attribute("instrument", instruments());
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write radial model
    if (radial()) {
        GXmlElement* rad = src->element("radialModel", 0);
        radial()->write(*rad);
    }

    // Write temporal model
    if (temporal()) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("lightcurve", 0);
            temporal()->write(*temp);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
 * @exception GCTAException::model_invalid_radial
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
        throw GCTAException::model_invalid_radial(G_XML_RADIAL, type);
    }

    // Return pointer
    return ptr;
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
GModelSpectral* GCTAModelRadialAcceptance::xml_spectral(const GXmlElement& spectral) const
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
GModelTemporal* GCTAModelRadialAcceptance::xml_temporal(const GXmlElement& temporal) const
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
 * is done analytically using the "cta_roi_arclength" support function.
 ***************************************************************************/
double GCTAModelRadialAcceptance::roi_kern::eval(double offset)
{
    // Circumvent const correctness
    GCTAModelRadial* radial = const_cast<GCTAModelRadial*>(m_parent);

    // Initialise Npred value
    double value = 0.0;
    
    // Continue only if offset > 0
    if (offset > 0.0) {

        // Get arclength for given radius in radians.
        double phi = cta_roi_arclength(offset, m_dist, m_cosdist, m_sindist,
                                       m_roi, m_cosroi);

        // Get kernel value if phi > 0
        if (phi > 0.0) {
            value = radial->eval(offset*gammalib::rad2deg) * phi *
                    std::sin(offset);
        }

    } // endif: offset was positive

    // Return
    return value;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
