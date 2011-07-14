/***************************************************************************
 *            GModelSky.hpp  -  Abstract virtual sky model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GModelSky.cpp
 * @brief Abstract sky model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GEvent.hpp"
#include "GPointing.hpp"
#include "GResponse.hpp"
#include "GObservation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NPRED           "GModelSky::npred(GEnergy&, GTime&, GObservation&)"
#define G_XML_SPATIAL                  "GModelSky::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL                "GModelSky::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL                "GModelSky::xml_temporal(GXmlElement&)"
#define G_SPATIAL             "GModelSky::spatial(GEvent&, GEnergy&, GTime&,"\
                                                       " GObservation, bool)"
#define G_SPECTRAL                     "GModelSky::spectral(GEvent&, GTime&," \
                                                      " GObservation&, bool)"
#define G_TEMPORAL        "GModelSky::temporal(GEvent&, GObservation&, bool)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DUMP_MC 0                                 //!< Dump MC information


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSky::GModelSky(void) : GModel()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GModelSky::GModelSky(const GXmlElement& xml) : GModel(xml)
{
    // Initialise private members for clean destruction
    init_members();

    // Get pointers on spectrum and spatial model
    GXmlElement* spec = (GXmlElement*)xml.element("spectrum", 0);
    GXmlElement* spat = (GXmlElement*)xml.element("spatialModel", 0);

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(*spat);
    m_spectral = xml_spectral(*spec);
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Sky model.
 ***************************************************************************/
GModelSky::GModelSky(const GModelSky& model) : GModel(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct sky model from spatial and spectral XML elements
 *
 * @param[in] spatial Spatial XML element.
 * @param[in] spectral Spectral XML element.
 ***************************************************************************/
GModelSky::GModelSky(const GXmlElement& spatial, const GXmlElement& spectral)
                                                                   : GModel()
{
    // Initialise private members for clean destruction
    init_members();

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(spatial);
    m_spectral = xml_spectral(spectral);
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSky::~GModelSky(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Sky model.
 ***************************************************************************/
GModelSky& GModelSky::operator= (const GModelSky& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModel::operator=(model);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns value of source model
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * This method evaluates the factorized source model for a given true photon
 * arrival direction, true photon energy and true photon arrival time.
 ***************************************************************************/
double GModelSky::value(const GSkyDir& srcDir, const GEnergy& srcEng,
                        const GTime& srcTime)
{
    // Initialise source model
    double source = 1.0;

    // Evaluate source model
    if (m_spatial  != NULL) source *= m_spatial->eval(srcDir);
    if (m_spectral != NULL) source *= m_spectral->eval(srcEng);
    if (m_temporal != NULL) source *= m_temporal->eval(srcTime);

    // Return
    return source;
}


/***********************************************************************//**
 * @brief Returns parameter gradients of source model
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * This method returns the parameter gradients of the factorized source model
 * for a given true photon arrival direction, true photon energy and true
 * photon arrival time.
 ***************************************************************************/
GVector GModelSky::gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                             const GTime& srcTime)
{
    // Evaluate source model gradients
    if (m_spatial  != NULL) m_spatial->eval_gradients(srcDir);
    if (m_spectral != NULL) m_spectral->eval_gradients(srcEng);
    if (m_temporal != NULL) m_temporal->eval_gradients(srcTime);

    // Set vector of gradients
    GVector gradients;
    if (size() > 0) {
        gradients = GVector(size());
        for (int i = 0; i < size(); ++i)
            gradients[i] = m_pars[i]->gradient();
    }

    // Return gradients
    return gradients;
}


/***********************************************************************//**
 * @brief Evaluate model
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * This method evaluates the source model for a given event within a given
 * observation.
 ***************************************************************************/
double GModelSky::eval(const GEvent& event, const GObservation& obs) const
{
    // Evaluate function
    double value = temporal(event, obs, false);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate model and parameter gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * This method evaluates the source model and model gradients for a given
 * event within a given observation.
 ***************************************************************************/
double GModelSky::eval_gradients(const GEvent& event, 
                                 const GObservation& obs) const
{
    // Evaluate function
    double value = temporal(event, obs, true);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated sky model
 *
 * @param[in] obsEng Measured photon energy.
 * @param[in] obsTime Measured photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GException::no_response
 *            No valid instrument response function defined.
 *
 * Computes
 * \f[N"_{\rm pred} = \int_{\rm ROI}
 *    S(\vec{p}, E, t) PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *    {\rm d}\vec{p'}\f]
 * where
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t)\f$ is the point
 * spread function,
 * \f$\vec{p'}\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy,
 * \f$t\f$ is the true photon arrival time, and
 * \f$d\f$ is the instrument pointing.
 *
 * \f${\rm ROI}\f$ is the region of interest that is stored in the
 * GObservation::m_roi member. The integration over the ROI is performed
 * by the GResponse::npred() method.
 *
 * @todo The actual method is only correct if no energy and time dispersion
 *       exists. For the moment we set srcEng=obsEng and srcTime=obsTime.
 *       Formally, Equation (2) of the instrument document has to be
 *       computed, which is an integration over source energy, time
 *       and arrival direction. For the moment, only the integration over
 *       arrival direction is performed by GResponse::npred().
 ***************************************************************************/
double GModelSky::npred(const GEnergy& obsEng, const GTime& obsTime,
                        const GObservation& obs) const
{
    // Initialise result
    double npred = 0.0;

    // Continue only if model is valid)
    if (valid_model()) {

        // Get response function
        GResponse* rsp = obs.response();
        if (rsp == NULL)
            throw GException::no_response(G_NPRED);

        // Here we make the simplifying approximations
        // srcEng=obsEng and srcTime=obsTime. To be fully correct we should
        // integrate over true energy and true time here ... at least true
        // time if we want to consider energy dispersion ...
        GEnergy srcEng  = obsEng;
        GTime   srcTime = obsTime;

        // Compute response components
        double npred_spatial  = rsp->npred(*this, srcEng, srcTime, obs);
        double npred_spectral = spectral()->eval(srcEng);
        double npred_temporal = temporal()->eval(srcTime);
                
        // Compute response
        npred = npred_spatial * npred_spectral * npred_temporal;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnan(npred) || isinfinite(npred)) {
            std::cout << "*** ERROR: GModelSky::npred:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", npred_spatial=" << npred_spatial;
            std::cout << ", npred_spectral=" << npred_spectral;
            std::cout << ", npred_temporal=" << npred_temporal;
            std::cout << ", srcEng=" << srcEng;
            std::cout << ", srcTime=" << srcTime;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: model was valid

    // Return npred
    return npred;
}


/***********************************************************************//**
 * @brief Read sky model from XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GModelSky::read(const GXmlElement& xml)
{
    // Clear sky model
    clear();

    // Get pointers on spectrum and spatial model
    GXmlElement* spec = (GXmlElement*)xml.element("spectrum", 0);
    GXmlElement* spat = (GXmlElement*)xml.element("spatialModel", 0);

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(*spat);
    m_spectral = xml_spectral(*spec);
    m_temporal = temporal.clone();

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
 * @param[in] xml Source library.
 ***************************************************************************/
void GModelSky::write(GXmlElement& xml) const
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
        if (spatial()  != NULL) src->append(new GXmlElement("spatialModel"));
        xml.append(src);
    }

    // Set model attributes
    src->attribute("name", name());
    if (spatial() != NULL) {
        if (spatial()->type() == "PointSource")
            src->attribute("type", "PointSource");
        else
            src->attribute("type", "DiffuseSource");
    }
    else
        src->attribute("type", "unknown");
    std::string instruments = this->instruments();
    if (instruments.length() > 0)
        src->attribute("instrument", instruments);

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = (GXmlElement*)src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write spatial model
    if (spatial() != NULL) {
        GXmlElement* spat = (GXmlElement*)src->element("spatialModel", 0);
        spatial()->write(*spat);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return simulated list of photons
 *
 * @param[in] area Simulation surface area (cm2).
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (deg).
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] tmin Minimum photon arrival time.
 * @param[in] tmax Maximum photon arrival time.
 * @param[in] ran Random number generator. 
 *
 * This method returns a list of photons that has been derived by Monte Carlo
 * simulation from the model. A simulation region is define by specification
 * of a simulation cone (a circular region on the sky),
 * of an energy range [emin, emax], and
 * of a time interval [tmin, tmax].
 * The simulation cone may eventually cover the entire sky (by setting
 * the radius to 180 degrees), yet simulations will be more efficient if
 * only the sky region will be simulated that is actually observed by the
 * telescope.
 *
 * @todo Check usage for diffuse models
 * @todo Implement photon arrival direction simulation for diffuse models
 * @todo Implement unique model ID to assign as Monte Carlo ID
 ***************************************************************************/
GPhotons GModelSky::mc(const double& area,
                       const GSkyDir& dir,  const double&  radius,
                       const GEnergy& emin, const GEnergy& emax,
                       const GTime&   tmin, const GTime&   tmax,
                       GRan& ran) const
{
    // Allocate photons
    GPhotons photons;

    // Continue only if model is valid)
    if (valid_model()) {

        // Get point source pointer
        GModelSpatialPtsrc* ptsrc = dynamic_cast<GModelSpatialPtsrc*>(m_spatial);

        // Check if model will produce any photons in the specified
        // simulation region. If the model is a point source we check if the
        // source is located within the simulation code. If the model is a
        // diffuse source we check if the source overlaps with the simulation
        // code
        bool use_model = true;
        if (ptsrc != NULL) {
            if (dir.dist(ptsrc->dir()) > radius)
                use_model = false;
        }
        else {
            //TODO
        }

        // Continue only if model overlaps with simulation region
        if (use_model) {

            // Compute flux within [emin, emax] in model from spectral
            // component (units: ph/cm2/s)
            double flux = m_spectral->flux(emin, emax);

            // Derive expecting counting rate within simulation surface
            // (units: ph/s)
            double rate = flux * area;

            // Debug option: dump rate
            #if G_DUMP_MC
            std::cout << "GModelSky::mc(\"" << name() << "\": ";
            std::cout << "flux=" << flux << " ph/cm2/s, ";
            std::cout << "rate=" << rate << " ph/s)" << std::endl;
            #endif

            // Get photon arrival times from temporal model
            GTimes times = m_temporal->mc(rate, tmin, tmax, ran);

            // Reserve space for photons
            if (times.size() > 0)
                photons.reserve(times.size());

            // Loop over photons
            for (int i = 0; i < times.size(); ++i) {

                // Allocate photon
                GPhoton photon;

                // Set photon arrival time
                photon.time(times[i]);

                // Set incident photon direction
                photon.dir(m_spatial->mc(ran));

                // Set photon energy
                photon.energy(m_spectral->mc(emin, emax, ran));

                // Append photon
                photons.push_back(photon);

            } // endfor: looped over photons

        } // endif: model was used
    } // endif: model was valid

    // Return photon list
    return photons;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSky::init_members(void)
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
 * @param[in] model Sky model.
 ***************************************************************************/
void GModelSky::copy_members(const GModelSky& model)
{
    // Clone spatial and spectral models
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
void GModelSky::free_members(void)
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
 * @brief Set parameter pointers
 *
 * Gathers all parameter pointers from the model.
 ***************************************************************************/
void GModelSky::set_pointers(void)
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
        for (int i = 0; i < n_spatial; ++i)
            m_pars.push_back(&((*spatial())[i]));

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
 * @brief Construct spatial model from XML element
 *
 * @param[in] spatial XML element containing spatial model information.
 *
 * @exception GException::model_invalid_spatial
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelSpatial* GModelSky::xml_spatial(const GXmlElement& spatial) const
{
    // Get spatial model type
    std::string type = spatial.attribute("type");

    // Get spatial model
    GModelSpatialRegistry registry;
    GModelSpatial*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL)
        ptr->read(spatial);

    // ... otherwise throw an exception
    else
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct spectral model from XML element
 *
 * @param[in] spectral XML element containing spectral model information.
 *
 * @exception GException::model_invalid_spectral
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelSpectral* GModelSky::xml_spectral(const GXmlElement& spectral) const
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
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelTemporal* GModelSky::xml_temporal(const GXmlElement& temporal) const
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
 * @brief Returns spatial model component
 *
 * @param[in] event Observed event.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @param[in] grad Evaluate gradients.
 *
 * @exception GException::no_response
 *            Observation has no valid instrument response
 *
 * This method computes the spatial model component for a given true photon
 * energy and arrival time.
 ***************************************************************************/
double GModelSky::spatial(const GEvent& event,
                          const GEnergy& srcEng, const GTime& srcTime,
                          const GObservation& obs, bool grad) const
{
    // Initialise result
    double value = 0.0;

    // Continue only if the model has a spatial component
    if (m_spatial != NULL) {

        // Get response function
        GResponse* rsp = obs.response();
        if (rsp == NULL)
            throw GException::no_response(G_SPATIAL);

        // Get IRF value. This method returns the spatial component of the
        // source model.
        double irf = rsp->irf(event, *this, srcEng, srcTime, obs);

        // Case A: evaluate gradients
        if (grad) {

            // Evaluate source model
            double spec = (spectral() != NULL) ? spectral()->eval_gradients(srcEng)  : 1.0;
            double temp = (temporal() != NULL) ? temporal()->eval_gradients(srcTime) : 1.0;

            // Set value
            value = spec * temp * irf;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (isnan(value) || isinfinite(value)) {
                std::cout << "*** ERROR: GModelSky::spatial:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (value=" << value;
                std::cout << ", spec=" << spec;
                std::cout << ", temp=" << temp;
                std::cout << ", irf=" << irf;
                std::cout << ")" << std::endl;
            }
            #endif

            // Multiply factors to spectral gradients
            if (spectral() != NULL) {
                double fact = temp * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < spectral()->size(); ++i)
                        (*spectral())[i].gradient((*spectral())[i].gradient() * fact);
                }
            }

            // Multiply factors to temporal gradients
            if (temporal() != NULL) {
                double fact = spec * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < temporal()->size(); ++i)
                        (*temporal())[i].gradient((*temporal())[i].gradient() * fact);
                }
            }

        } // endif: gradient evaluation has been requested

        // Case B: evaluate no gradients
        else {

            // Evaluate source model
            double spec = (m_spectral != NULL) ? m_spectral->eval(srcEng) : 1.0;
            double temp = (m_temporal != NULL) ? m_temporal->eval(srcTime) : 1.0;

            // Set value
            value = spec * temp * irf;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (isnan(value) || isinfinite(value)) {
                std::cout << "*** ERROR: GModelSky::spatial:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (value=" << value;
                std::cout << ", spec=" << spec;
                std::cout << ", temp=" << temp;
                std::cout << ", irf=" << irf;
                std::cout << ")" << std::endl;
            }
            #endif

        }

    } // endif: Gamma-ray source model had a spatial component

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over spectral component
 *
 * @param[in] event Observed event.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @param[in] grad Evaluate gradients.
 *
 * @exception GException::no_response
 *            Observation has no valid instrument response
 * @exception GException::feature_not_implemented
 *            Energy integration not yet implemented
 *
 * This method integrates the source model over the spectral component. If
 * the response function has no energy dispersion then no spectral
 * integration is needed and the observed photon energy is identical to the
 * true photon energy.
 *
 * @todo Needs implementation of spectral integration to handle energy
 *       dispersion.
 ***************************************************************************/
double GModelSky::spectral(const GEvent& event, const GTime& srcTime,
                           const GObservation& obs, bool grad) const
{
    // Initialise result
    double value = 0.0;

    // Get response function
    GResponse* rsp = obs.response();
    if (rsp == NULL)
        throw GException::no_response(G_SPECTRAL);

    // Determine if energy integration is needed
    bool integrate = rsp->hasedisp();

    // Case A: Integraion
    if (integrate) {
        throw GException::feature_not_implemented(G_SPECTRAL);
    }

    // Case B: No integration (assume no energy dispersion)
    else
        value = spatial(event, event.energy(), srcTime, obs, grad);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnan(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSky::spectral:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", event=" << event;
        std::cout << ", srcTime=" << srcTime;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over temporal component
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] grad Evaluate gradients.
 *
 * @exception GException::no_response
 *            Observation has no valid instrument response
 * @exception GException::feature_not_implemented
 *            Temporal integration not yet implemented
 *
 * This method integrates the source model over the temporal component. If
 * the response function has no time dispersion then no temporal integration
 * is needed and the observed photon arrival time is identical to the true
 * photon arrival time.
 *
 * @todo Needs implementation of temporal integration to handle time
 *       dispersion.
 ***************************************************************************/
double GModelSky::temporal(const GEvent& event, const GObservation& obs,
                           bool grad) const
{
    // Initialise result
    double value = 0.0;

    // Get response function
    GResponse* rsp = obs.response();
    if (rsp == NULL)
        throw GException::no_response(G_TEMPORAL);

    // Determine if time integration is needed
    bool integrate = rsp->hastdisp();

    // Case A: Integraion
    if (integrate) {
        throw GException::feature_not_implemented(G_TEMPORAL);
    }

    // Case B: No integration (assume no time dispersion)
    else
        value = spectral(event, event.time(), obs, grad);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnan(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSky::temporal:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", event=" << event;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 ***************************************************************************/
bool GModelSky::valid_model(void) const
{
    // Set result
    bool result = ((m_spatial  != NULL) &&
                   (m_spectral != NULL) &&
                   (m_temporal != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GModelSky::print_model(void) const
{
    // Initialise result string
    std::string result;

    // Determine number of parameters per type
    int n_spatial  = (m_spatial  != NULL) ? m_spatial->size()  : 0;
    int n_spectral = (m_spectral != NULL) ? m_spectral->size() : 0;
    int n_temporal = (m_temporal != NULL) ? m_temporal->size() : 0;

    // Append header
    result.append(parformat("Name")+name());
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
    if (n_spatial > 0) {
        result.append("\""+spatial()->type()+"\"");
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
    result.append("\n"+parformat("Number of spatial par's")+str(n_spatial));
    for (int i = 0; i < n_spatial; ++i)
        result.append("\n"+(*spatial())[i].print());
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
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
