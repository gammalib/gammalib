/***************************************************************************
 *     GCTAModelRadialAcceptance.cpp  -  Radial acceptance model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAModelRadialAcceptance.cpp
 * @brief Radial acceptance model class implementation
 * @author J. Knodlseder
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
#define G_DUMP_MC 0                                 //!< Dump MC information


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
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
 * @brief Copy constructor
 *
 * @param[in] model Radial acceptance model.
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
    GTime time; // not used
    GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing(time));
    if (pnt == NULL)
        throw GCTAException::no_pointing(G_EVAL);

    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Compute offset angle (in degrees)
    double offset = cta_dir->dist_deg(pnt->dir());

    // Initialise value
    double value = 1.0;

    // Evaluate function
    if (radial()   != NULL) value *= radial()->eval(offset);
    if (spectral() != NULL) value *= spectral()->eval(event.energy());
    if (temporal() != NULL) value *= temporal()->eval(event.time());

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
    GTime time; // not used
    GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing(time));
    if (pnt == NULL)
        throw GCTAException::no_pointing(G_EVAL);

    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Compute offset angle (in degrees)
    double offset = cta_dir->dist_deg(pnt->dir());

    // Evaluate function and gradients
    double rad  = (radial()   != NULL) ? radial()->eval_gradients(offset) : 1.0;
    double spec = (spectral() != NULL) ? spectral()->eval_gradients(event.energy()) : 1.0;
    double temp = (temporal() != NULL) ? temporal()->eval_gradients(event.time()) : 1.0;

    // Compute value
    double value = rad * spec * temp;

    // Multiply factors to radial gradients
    if (radial() != NULL) {
        double fact = spec * temp;
        if (fact != 1.0) {
            for (int i = 0; i < radial()->size(); ++i)
                (*radial())[i].gradient( (*radial())[i].gradient() * fact );
        }
    }

    // Multiply factors to spectral gradients
    if (spectral() != NULL) {
        double fact = rad * temp;
        if (fact != 1.0) {
            for (int i = 0; i < spectral()->size(); ++i)
                (*spectral())[i].gradient( (*spectral())[i].gradient() * fact );
        }
    }

    // Multiply factors to temporal gradients
    if (temporal() != NULL) {
        double fact = rad * spec;
        if (fact != 1.0) {
            for (int i = 0; i < temporal()->size(); ++i)
                (*temporal())[i].gradient( (*temporal())[i].gradient() * fact );
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
 * Spatially integrates the data model for a given measured event energy
 * and event time.
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
        if (events == NULL)
            throw GException::no_list(G_NPRED);

        // Get CTA pointing direction
        GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing(obsTime));
        if (pnt == NULL)
            throw GCTAException::no_pointing(G_NPRED);

        // Get ROI radius in radians
        double roi_radius = events->roi().radius() * deg2rad;

        // Get distance from ROI centre in radians
        double roi_distance = events->roi().centre().dist(pnt->dir());

        // Setup integration function
        GCTAModelRadialAcceptance::roi_kern integrand(radial(), roi_radius, roi_distance);

        // Setup integrator
        GIntegral integral(&integrand);

        // Spatially integrate radial component
        npred = integral.romb(0.0, roi_radius);

        // Multiply in spectral and temporal components
        npred *= spectral()->eval(obsEng);
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
 * @exception GCTAException::no_pointing
 *            No pointing found
 *
 * This method requires that 
 *   the pointing,
 *   the energy boundaries, and
 *   the good time interval
 * of the observation have been set up previously.
 ***************************************************************************/
GCTAEventList* GCTAModelRadialAcceptance::mc(const GObservation& obs, 
                                             GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Extract CTA pointing direction at beginning of observation
        GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing(obs.events()->tstart()));
        if (pnt == NULL)
            throw GCTAException::no_pointing(G_MC);

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

            // Derive expecting rate (units: cts/s)
            double rate = flux * area;

            // Debug option: dump rate
            #if G_DUMP_MC
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
                                              obs.events()->gti().tstop(itime), ran);

                // Reserve space for events
                if (times.size() > 0)
                    list->reserve(times.size());

                // Loop over events
                for (int i = 0; i < times.size(); ++i) {

                    // Set event direction
                    GCTAInstDir dir = radial()->mc(pnt_dir, ran);

                    // Set event energy
                    GEnergy energy = spectral()->mc(obs.events()->ebounds().emin(ieng),
                                                    obs.events()->ebounds().emax(ieng),
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
    GXmlElement* rad  = NULL;
    GXmlElement* spec = NULL;
    GXmlElement* temp = NULL;

    // Get pointers on spectrum and radial model
    rad  = (GXmlElement*)xml.element("radialModel", 0);
    spec = (GXmlElement*)xml.element("spectrum", 0);

    // Clone radial and spectral models
    m_radial   = xml_radial(*rad);
    m_spectral = xml_spectral(*spec);

    // Optionally get temporal model
    try {
        temp = (GXmlElement*)xml.element("lightcurve", 0);
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

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GCTAModelRadialAcceptance::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = (GXmlElement*)xml.element("source", k);
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = new GXmlElement("source");
        src->attribute("name") = name();
        if (spectral() != NULL) src->append(new GXmlElement("spectrum"));
        if (radial()   != NULL) src->append(new GXmlElement("radialModel"));
        xml.append(src);
    }

    // Set model type, name and optionally instruments
    src->attribute("name", name());
    src->attribute("type", type());
    if (instruments().length() > 0)
        src->attribute("instrument", instruments());

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = (GXmlElement*)src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write radial model
    if (radial()) {
        GXmlElement* rad = (GXmlElement*)src->element("radialModel", 0);
        radial()->write(*rad);
    }

    // Write temporal model
    if (temporal()) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = (GXmlElement*)src->element("lightcurve", 0);
            temporal()->write(*temp);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GCTAModelRadialAcceptance::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAModelRadialAcceptance ===");

    // Determine number of parameters per type
    int n_radial   = (radial()   != NULL) ? radial()->size()   : 0;
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

    // Append name and instruments
    result.append("\n"+parformat("Name")+name());
    result.append("\n"+parformat("Instruments"));
    if (m_instruments.size() > 0) {
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (i > 0)
                result.append(", ");
            result.append(m_instruments[i]);
        }
    }
    else
        result.append("all");

    // Append model type
    result.append("\n"+parformat("Model type"));
    if (n_radial > 0) {
        result.append("\""+radial()->type()+"\"");
        if (n_spectral > 0 || n_temporal > 0)
            result.append(" * ");
    }
    if (n_spectral > 0) {
        result.append("\""+spectral()->type()+"\"");
        if (n_temporal > 0)
            result.append(" * ");
    }
    if (n_temporal > 0)
        result.append("\""+temporal()->type()+"\"");

    // Append parameters
    result.append("\n"+parformat("Number of parameters")+str(size()));
    result.append("\n"+parformat("Number of radial par's")+str(n_radial));
    for (int i = 0; i < n_radial; ++i)
        result.append("\n"+(*radial())[i].print());
    result.append("\n"+parformat("Number of spectral par's")+str(n_spectral));
    for (int i = 0; i < n_spectral; ++i)
        result.append("\n"+(*spectral())[i].print());
    result.append("\n"+parformat("Number of temporal par's")+str(n_temporal));
    for (int i = 0; i < n_temporal; ++i)
        result.append("\n"+(*temporal())[i].print());

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
        for (int i = 0; i < n_radial; ++i)
            m_pars.push_back(&((*radial())[i]));

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i)
            m_pars.push_back(&((*spectral())[i]));

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i)
            m_pars.push_back(&((*temporal())[i]));

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
    if (ptr != NULL)
        ptr->read(radial);

    // ... otherwise throw an exception
    else
        throw GCTAException::model_invalid_radial(G_XML_RADIAL, type);

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
    if (ptr != NULL)
        ptr->read(spectral);

    // ... otherwise throw an exception
    else
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);

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
    if (ptr != NULL)
        ptr->read(temporal);

    // ... otherwise throw an exception
    else
        throw GException::model_invalid_temporal(G_XML_TEMPORAL, type);

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
    GCTAModelRadial* radial = (GCTAModelRadial*)m_parent;

    // Get arclength for given radius in radians.
    double phi = cta_roi_arclength(offset, m_dist, m_cosdist, m_sindist,
                                   m_roi, m_cosroi);

    // Get PSF value
    double value = radial->eval(offset*rad2deg) * phi * sin(offset);

    // Return
    return value;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
