/***************************************************************************
 *                    GModelSky.cpp - Sky model class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @brief Sky model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GIntegral.hpp"
#include "GModelRegistry.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GSource.hpp"
#include "GResponse.hpp"

/* __ Globals ____________________________________________________________ */
const GModelSky         g_pointsource_seed("PointSource");
const GModelSky         g_extendedsource_seed("ExtendedSource");
const GModelSky         g_diffusesource_seed("DiffuseSource");
const GModelRegistry    g_pointsource_registry(&g_pointsource_seed);
const GModelRegistry    g_extendedsource_registry(&g_extendedsource_seed);
const GModelRegistry    g_diffusesource_registry(&g_diffusesource_seed);

/* __ Method name definitions ____________________________________________ */
#define G_NPRED           "GModelSky::npred(GEnergy&, GTime&, GObservation&)"
#define G_XML_SPATIAL                  "GModelSky::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL                "GModelSky::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL                "GModelSky::xml_temporal(GXmlElement&)"
#define G_INTEGRATE_TIME  "GModelSky::integrate_time(GEvent&, GObservation&,"\
                                                                     " bool)"
#define G_INTEGRATE_ENERGY     "GModelSky::integrate_energy(GEvent&, GTime&," \
                                                      " GObservation&, bool)"
#define G_INTEGRATE_DIR "GModelSky::integrate_dir(GEvent&, GEnergy&, GTime&,"\
                                                       " GObservation, bool)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC                                  //!< Dump MC information
//#define G_DUMP_MC_DETAIL                  //!< Dump detailed MC information


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Construct an empty sky model. An empty sky model has no type.
 ***************************************************************************/
GModelSky::GModelSky(void) : GModel()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Type constructor
 *
 * @param[in] type Model type.
 *
 * Construct an empty sky model of the specified @p type. This constructor
 * does basically the same than the void constructor with the only difference
 * that the m_type member of the class is set to the specified @p type
 * string.
 ***************************************************************************/
GModelSky::GModelSky(const std::string& type) : GModel()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a sky model from the model information that is found in the
 * @p xml element. Only the spatial and spectral information is handled by
 * this constructor. For the temporal component, a constant model of type
 * GModelTemporalConst will be allocated.
 ***************************************************************************/
GModelSky::GModelSky(const GXmlElement& xml) : GModel(xml)
{
    // Initialise members
    init_members();

    // Get pointers on spectrum and spatial model
    const GXmlElement* spec = xml.element("spectrum", 0);
    const GXmlElement* spat = xml.element("spatialModel", 0);

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(*spat);
    m_spectral = xml_spectral(*spec);
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Set model type dependent on spatial model type
    set_type();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct sky model from spatial and spectral XML elements
 *
 * @param[in] spatial Spatial XML element.
 * @param[in] spectral Spectral XML element.
 *
 * Constructs a sky model from the @p spatial and @p spectral information
 * that is found in the respective XML elements. For the temporal component,
 * a constant model of type GModelTemporalConst will be allocated.
 ***************************************************************************/
GModelSky::GModelSky(const GXmlElement& spatial, const GXmlElement& spectral) :
           GModel()
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

    // Set model type dependent on spatial model type
    set_type();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct sky model from spatial, spectral and temporal XML elements
 *
 * @param[in] spatial Spatial XML element.
 * @param[in] spectral Spectral XML element.
 * @param[in] temporal Temporal XML element.
 *
 * Constructs a sky model from the @p spatial, @p spectral and @p temporal
 * information that is found in the respective XML elements.
 ***************************************************************************/
GModelSky::GModelSky(const GXmlElement& spatial,
                     const GXmlElement& spectral,
                     const GXmlElement& temporal) :
           GModel()
{
    // Initialise private members for clean destruction
    init_members();

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(spatial);
    m_spectral = xml_spectral(spectral);
    m_temporal = xml_temporal(temporal);

    // Set parameter pointers
    set_pointers();

    // Set model type dependent on spatial model type
    set_type();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct sky model from spatial and spectral components
 *
 * @param[in] spatial Spatial model component.
 * @param[in] spectral Spectral model component.
 *
 * Constructs a sky model from @p spatial and @p spectral model components.
 * For the temporal component, a constant model of type GModelTemporalConst
 * will be allocated.
 ***************************************************************************/
GModelSky::GModelSky(const GModelSpatial&  spatial,
                     const GModelSpectral& spectral) : 
           GModel()
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

    // Set model type dependent on spatial model type
    set_type();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct sky model from spatial, spectral and temporal components
 *
 * @param[in] spatial Spatial model component.
 * @param[in] spectral Spectral model component.
 * @param[in] temporal Temporal model component.
 *
 * Constructs a sky model from @p spatial, @p spectral and @p temporal model
 * components.
 ***************************************************************************/
GModelSky::GModelSky(const GModelSpatial& spatial,
                     const GModelSpectral& spectral,
                     const GModelTemporal& temporal) : 
           GModel()
{
    // Initialise members
    init_members();

    // Clone model components
    m_spatial  = spatial.clone();
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Set model type dependent on spatial model type
    set_type();

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
 * @return Sky model.
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
 * @brief Clear sky model
 *
 * This method properly resets the sky model to an initial state.
 ***************************************************************************/
void GModelSky::clear(void)
{
    // Free class members
    free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone sky model
 *
 * @return Pointer to deep copy of sky model.
 ***************************************************************************/
GModelSky* GModelSky::clone(void) const
{
    // Clone sky model
    return new GModelSky(*this);
}


/***********************************************************************//**
 * @brief Return value of sky model for a given photon
 *
 * @param[in] photon Photon.
 * @return Value of sky model.
 *
 * Returns the value of the sky model for a given @p photon. If no model
 * components exist the model will return a value of 1.
 *
 * @todo We probably should return a value of 0 is no model components exist.
 * But we may also make sure that the model never has NULL pointers, which
 * would avoid all the pointer checks.
 ***************************************************************************/
double GModelSky::value(const GPhoton& photon)
{
    // Initialise source model
    double source = 1.0;

    // Evaluate source model
    if (m_spatial  != NULL) source *= m_spatial->eval(photon);
    if (m_spectral != NULL) source *= m_spectral->eval(photon.energy(),
                                                       photon.time());
    if (m_temporal != NULL) source *= m_temporal->eval(photon.time());

    // Return
    return source;
}


/***********************************************************************//**
 * @brief Return parameter gradients of sky model for a given photon
 *
 * @param[in] photon Photon.
 * @return Vector of parameter gradients
 *
 * Returns a vector of parameter gradients for the sky model for a given
 * @p photon. If there are no parameters in the sky model, an empty vector
 * will be returned.
 ***************************************************************************/
GVector GModelSky::gradients(const GPhoton& photon)
{
    // Evaluate source model gradients
    if (m_spatial  != NULL) m_spatial->eval_gradients(photon);
    if (m_spectral != NULL) m_spectral->eval_gradients(photon.energy(),
                                                       photon.time());
    if (m_temporal != NULL) m_temporal->eval_gradients(photon.time());

    // Set vector of gradients
    GVector gradients;
    if (size() > 0) {
        gradients = GVector(size());
        for (int i = 0; i < size(); ++i) {
            gradients[i] = m_pars[i]->factor_gradient();
        }
    }

    // Return gradients
    return gradients;
}


/***********************************************************************//**
 * @brief Evaluate sky model for a given event of an observation
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Value of sky model
 *
 * Evalutes the value of the sky model for an @p event of a specific
 * observation @p obs. This method will subsequently call the following
 * methods:
 * - integrate_time() (performs integral over time dispersion)
 * - spectral() (performs integral over energy dispersion)
 * - spatial() (performs integral over spatial dispersion)
 ***************************************************************************/
double GModelSky::eval(const GEvent& event, const GObservation& obs) const
{
    // Evaluate function
    double value = integrate_time(event, obs, false);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate sky model and parameter gradients for a given event of an
 *        observation
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Value of sky model
 *
 * Evalutes the value of the sky model and of the parameter for an @p event
 * of a specific observation @p obs. This method will subsequently call the
 * following methods:
 * - integrate_time() (performs integral over time dispersion)
 * - spectral() (performs integral over energy dispersion)
 * - spatial() (performs integral over spatial dispersion)
 *
 * While the value of the sky model is returned by the method, the parameter
 * gradients are set as GModelPar members.
 ***************************************************************************/
double GModelSky::eval_gradients(const GEvent& event, 
                                 const GObservation& obs) const
{
    // Evaluate function
    double value = integrate_time(event, obs, true);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated sky model
 *
 * @param[in] obsEng Measured photon energy.
 * @param[in] obsTime Measured photon arrival time.
 * @param[in] obs Observation.
 * @return Spatially integrated sky model.
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
 * The method takes care of any instrument dependent scale factors. These
 * scale factors will be applied to the predicted number of model counts.
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
        const GResponse& rsp = obs.response();

        // Here we make the simplifying approximations
        // srcEng=obsEng and srcTime=obsTime. To be fully correct we should
        // integrate over true energy and true time here ... at least true
        // time if we want to consider energy dispersion ...
        GEnergy srcEng  = obsEng;
        GTime   srcTime = obsTime;

        // Set source
        GSource source(this->name(), m_spatial, srcEng, srcTime);

        // Compute response components
        double npred_spatial  = rsp.npred(source, obs);
        double npred_spectral = spectral()->eval(srcEng, srcTime);
        double npred_temporal = temporal()->eval(srcTime);

        // Compute response
        npred = npred_spatial * npred_spectral * npred_temporal;

        // If required, apply instrument specific model scaling
        if (!m_scales.empty()) {
            npred *= scale(obs.instrument()).value();
        }

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
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
 *     </source>
 *
 * Optionally, the model may also contain scale parameters following the
 * format:
 *
 *     <source name=".." type=".." instrument=".." id="..">
 *       <spectrum type="..">
 *         ..
 *       </spectrum>
 *       <spatialModel type="..">
 *         ..
 *       </spatialModel>
 *       <scaling>
 *         <instrument name=".." scale="1.0" min="0.1" max="10.0" value="1.0" free="0"/>
 *         <instrument name=".." scale="1.0" min="0.1" max="10.0" value="0.5" free="0"/>
 *       </scaling>
 *     </source>
 *
 * (see GModel::read_scales() for more information on instrument scales).
 ***************************************************************************/
void GModelSky::read(const GXmlElement& xml)
{
    // Clear sky model
    clear();

    // Get pointers on spectrum and spatial model
    const GXmlElement* spec = xml.element("spectrum", 0);
    const GXmlElement* spat = xml.element("spatialModel", 0);

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

    // Read instrument scales
    read_scales(xml);

    // Set observation identifiers
    ids(xml.attribute("id"));

    // Set parameter pointers
    set_pointers();

    // Set model type dependent on spatial model type
    set_type();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml Source library.
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
 *       <scaling>
 *         <instrument name=".." scale="1.0" min="0.1" max="10.0" value="1.0" free="0"/>
 *         <instrument name=".." scale="1.0" min="0.1" max="10.0" value="0.5" free="0"/>
 *       </scaling>
 *     </source>
 *
 * The scaling element will only be written optionally in case that instrument
 * dependent scaling factors exist (see GModel::write_scales() for more
 * information on instrument scales).
 ***************************************************************************/
void GModelSky::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = static_cast<GXmlElement*>(xml.element("source", k));
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (spatial()  != NULL) src->append(GXmlElement("spatialModel"));
    }

    // Set model attributes
    src->attribute("name", name());
    src->attribute("type", type());
    std::string instruments = this->instruments();
    if (instruments.length() > 0) {
        src->attribute("instrument", instruments);
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

    // Write spatial model
    if (spatial() != NULL) {
        GXmlElement* spat = src->element("spatialModel", 0);
        spatial()->write(*spat);
    }

    // Write instrument scales
    write_scales(*src);

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
 * @param[in,out] ran Random number generator.
 * @return List of photons
 *
 * Returns a list of photons that has been derived by Monte Carlo simulation
 * from the model. A simulation region is define by specification of
 * - a simulation cone, which is a circular region on the sky defined by
 *   a centre direction @p dir and a @p radius,
 * - an energy range [@p emin, @p emax], and
 * - a time interval [@p tmin, @p tmax].
 *
 * Only photons with parameters in the simulation region will be returned
 * by the method.
 *
 * The simulation cone may eventually cover the entire sky (by setting
 * the radius to 180 degrees), yet simulations will be more efficient if
 * only the sky region will be simulated that is actually observed by the
 * telescope.
 *
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

        // Determine the spatial model normalization within the simulation
        // cone and check whether the model will produce any photons in that
        // cone.
        double norm      = m_spatial->norm(dir, radius);
        bool   use_model = (norm > 0.0) ? true : false;

        // Continue only if model overlaps with simulation region
        if (use_model) {

            // Initialise de-allocation flag
            bool free_spectral = false;

            // Set pointer to spectral model
            GModelSpectral* spectral = m_spectral;

            // If the spectral model is a diffuse cube then create a node
            // function spectral model that is the product of the diffuse
            // cube node function and the spectral model evaluated at the
            // energies of the node function
            GModelSpatialDiffuseCube* cube = dynamic_cast<GModelSpatialDiffuseCube*>(m_spatial);
            if (cube != NULL) {

                // Set MC cone
                cube->set_mc_cone(dir, radius);

                // Allocate node function to replace the spectral component
                GModelSpectralNodes* nodes = new GModelSpectralNodes(cube->spectrum());
                for (int i = 0; i < nodes->nodes(); ++i) {
                    GEnergy energy    = nodes->energy(i);
                    double  intensity = nodes->intensity(i);
                    double  norm      = m_spectral->eval(energy, tmin);
                    nodes->intensity(i, norm*intensity);
                }
                
                // Signal that node function needs to be de-allocated later
                free_spectral = true;

                // Set the spectral model pointer to the node function
                spectral = nodes;

            } // endif: spatial model was a diffuse cube

            // Compute flux within [emin, emax] in model from spectral
            // component (units: ph/cm2/s)
            double flux = spectral->flux(emin, emax);

            // Derive expecting counting rate within simulation surface
            // (units: ph/s)
            double rate = flux * area * norm;

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GModelSky::mc(\"" << name() << "\": ";
            std::cout << "flux=" << flux << " ph/cm2/s, ";
            std::cout << "rate=" << rate << " ph/s)" << std::endl;
            #endif

            // Get photon arrival times from temporal model
            GTimes times = m_temporal->mc(rate, tmin, tmax, ran);

            // Debug option: dump number of times
            #if defined(G_DUMP_MC_DETAIL)
            std::cout << "  Times=" << times.size() << std::endl;
            #endif

            // Reserve space for photons
            if (times.size() > 0) {
                photons.reserve(times.size());
            }

            // Loop over photons
            for (int i = 0; i < times.size(); ++i) {

                // Debug option: dump photon index
                #if defined(G_DUMP_MC_DETAIL)
                std::cout << "  Photon=" << i << std::endl;
                #endif

                // Allocate photon
                GPhoton photon;

                // Set photon arrival time
                photon.time(times[i]);

                // Debug option: dump time
                #if defined(G_DUMP_MC_DETAIL)
                std::cout << "    Time=" << times[i] << std::endl;
                #endif

                // Set photon energy
                photon.energy(spectral->mc(emin, emax, photon.time(), ran));

                // Debug option: dump energy
                #if defined(G_DUMP_MC_DETAIL)
                std::cout << "    Energy=" << photon.energy() << std::endl;
                #endif

                // Set incident photon direction
                photon.dir(m_spatial->mc(photon.energy(), photon.time(), ran));

                // Debug option: dump direction
                #if defined(G_DUMP_MC_DETAIL)
                std::cout << "    Direction=" << photon.dir() << std::endl;
                #endif

                // Append photon
                if (dir.dist_deg(photon.dir()) <= radius) {
                    photons.append(photon);
                }

            } // endfor: looped over photons

            // Free spectral model if required
            if (free_spectral) delete spectral;

        } // endif: model was used
    } // endif: model was valid

    // Return photon list
    return photons;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSky::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSky ===");

        // Append model
        //result.append("\n"+print_model(chatter));

        // Determine number of parameters per type
        int n_spatial  = (m_spatial  != NULL) ? m_spatial->size()  : 0;
        int n_spectral = (m_spectral != NULL) ? m_spectral->size() : 0;
        int n_temporal = (m_temporal != NULL) ? m_temporal->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type")+type());

        // Append model components
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
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of spatial par's"));
        result.append(gammalib::str(n_spatial));
        for (int i = 0; i < n_spatial; ++i) {
            result.append("\n"+(*spatial())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of spectral par's"));
        result.append(gammalib::str(n_spectral));
        for (int i = 0; i < n_spectral; ++i) {
            result.append("\n"+(*spectral())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of temporal par's"));
        result.append(gammalib::str(n_temporal));
        for (int i = 0; i < n_temporal; ++i) {
            result.append("\n"+(*temporal())[i].print());
        }

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
void GModelSky::init_members(void)
{
    // Initialise members
    m_type.clear();
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
    // Copy attributes
    m_type = model.m_type;
    
    // Clone model components
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
 * Gathers all parameter pointers from the model into a flat array of model
 * pointers. 
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
 * @brief Set model type based on spatial model component
 *
 * Set the model type depending on the class type of the spatial model
 * component.
 *
 * @todo A method could be implemented in the GModelSpatial class that
 *       determines the model type. This is however not very critical.
 ***************************************************************************/
void GModelSky::set_type(void)
{
    // Clear model type
    m_type.clear();

    // Continue only if we have a spatial model component
    if (m_spatial != NULL) {
    
        // Is spatial model a point source?
        if (dynamic_cast<const GModelSpatialPointSource*>(m_spatial) != NULL) {
            m_type = "PointSource";
        }
        
        // ... otherwise, is spatial model a radial source?
        else if (dynamic_cast<const GModelSpatialRadial*>(m_spatial) != NULL) {
            m_type = "ExtendedSource";
        }

        // ... otherwise, is spatial model an elliptical source?
        else if (dynamic_cast<const GModelSpatialElliptical*>(m_spatial) != NULL) {
            m_type = "ExtendedSource";
        }

        // ... otherwise we have a diffuse model
        else {
            m_type = "DiffuseSource";
        }
    
    } // endif: there was a spatial model component

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
    if (ptr != NULL) {
        ptr->read(spatial);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);
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
 * @brief Compute event probability for the sky model
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] grad Evaluate gradients.
 * @return Probability of measuring the event
 *         (\f${\rm s}^{-1}\f$ \f${\rm MeV}^{-1}\f$ \f${\rm sr}^{-1}\f$).
 *
 * @exception GException::feature_not_implemented
 *            Temporal integration not yet implemented
 *
 * Computes the probability \f$P(\vec{p'}, E', t')\f$ of measuring an
 * @p event with instrument direction \f$\vec{p'}\f$, energy \f$E'\f$ at
 * time \f$t'\f$ from a source \f$S(\vec{p}, E, t)\f$ using
 *
 * \f[
 *    P(\vec{p'}, E', t') = \int_{0}^{t'+\Delta t}
 *                          \int_{E'-\Delta E}^{\infty}
 *                          \int_{\Omega} 
 *                          S(\vec{p}, E, t) \,
 *                          R(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *                          {\rm d}\vec{p} \, {\rm d}E \,{\rm d}t
 * \f]
 *
 * (in units of \f${\rm s}^{-1}\f$ \f${\rm MeV}^{-1}\f$ \f${\rm sr}^{-1}\f$).
 * The terms \f$\Delta t\f$ and \f$\Delta E\f$ account for the statistical
 * jitter related to the measurement process and are of the order of a few
 * times the rms in the time and energy measurements.
 *
 * The method performs the integration of the true time. The integration
 * over the true photon energy and photon arrival direction is performed
 * by integrate_energy() and integrate_dir(), respectively.
 *
 * @todo Needs implementation of temporal integration to handle time
 *       dispersion.
 ***************************************************************************/
double GModelSky::integrate_time(const GEvent& event,
                                 const GObservation& obs,
                                 bool grad) const
{
    // Initialise result
    double value = 0.0;

    // Get response function
    const GResponse& rsp = obs.response();

    // Determine if time integration is needed
    bool integrate = rsp.use_tdisp();

    // Case A: Integration
    if (integrate) {
        throw GException::feature_not_implemented(G_INTEGRATE_TIME);
    }

    // Case B: No integration (assume no time dispersion)
    else {
        value = integrate_energy(event, event.time(), obs, grad);
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSky::integrate_time:";
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
 * @brief Integrate sky model over true photon energy and arrival direction
 *
 * @param[in] event Observed event.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @param[in] grad Evaluate gradients.
 *
 * @exception GException::feature_not_implemented
 *            Energy integration not yet implemented
 *
 * Integrates the sky model over the true photon energy and arrival
 * direction using
 *
 * \f[
 *    \int_{E'-\Delta E}^{\infty}
 *    \int_{\Omega} 
 *    S(\vec{p}, E, t) \,
 *    R(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *    {\rm d}\vec{p} \, {\rm d}E
 * \f]
 *
 * The method performs the integration of the true photon energy. The
 * integration over the true photon arrival direction is performed by
 * integrate_dir().
 *
 * @todo Needs implementation of energy integration to handle energy
 *       dispersion.
 ***************************************************************************/
double GModelSky::integrate_energy(const GEvent& event,
                                   const GTime& srcTime,
                                   const GObservation& obs,
                                   bool grad) const
{
    // Initialise result
    double value = 0.0;

    // Get response function
    const GResponse& rsp = obs.response();

    // Determine if energy integration is needed
    bool integrate = rsp.use_edisp();

    // Case A: Integration
    if (integrate) {
    
        // Retrieve true energy boundaries
        GEbounds ebounds = rsp.ebounds_src(event.energy());
    
        // Loop over all boundaries
        for (int i = 0; i < ebounds.size(); ++i) {

            // Get boundaries in MeV
            double emin = ebounds.emin(i).MeV();
            double emax = ebounds.emax(i).MeV();

            // Continue only if valid
            if (emax > emin) {

                // Setup integration function
                GModelSky::edisp_kern integrand(this, event, srcTime, obs, grad);
                GIntegral integral(&integrand);

                // Set integration precision
                integral.eps(1.0e-3);

                // Do Romberg integration
                emin   = std::log(emin);
                emax   = std::log(emax);
                value += integral.romb(emin, emax);
    
            } // endif: interval was valid
        } // endfor: looped over intervals

    }

    // Case B: No integration (assume no energy dispersion)
    else {
        value = integrate_dir(event, event.energy(), srcTime, obs, grad);
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSky::integrate_energy:";
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
 * @brief Integrate sky model over true photon arrival direction
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
 * Integrates the sky model over the true photon arrival direction using
 *
 * \f[
 *    \int_{\Omega} 
 *    S(\vec{p}, E, t) \,
 *    R(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *    {\rm d}\vec{p}
 * \f]
 *
 * The method takes care of any instrument dependent scale factors. These
 * scale factors will be applied to the IRF so that they are correctly
 * taken into account in the spectral and temporal model gradient
 * computations.
 ***************************************************************************/
double GModelSky::integrate_dir(const GEvent&       event,
                                const GEnergy&      srcEng,
                                const GTime&        srcTime,
                                const GObservation& obs,
                                bool grad) const
{
    // Initialise result
    double value = 0.0;

    // Continue only if the model has a spatial component
    if (m_spatial != NULL) {

        // Get response function
        const GResponse& rsp = obs.response();

        // Set source
        GSource source(this->name(), m_spatial, srcEng, srcTime);
        
        // Get IRF value. This method returns the spatial component of the
        // source model.
        double irf = rsp.irf(event, source, obs);

        // If required, apply instrument specific model scaling
        if (!m_scales.empty()) {
            irf *= scale(obs.instrument()).value();
        }

        // Case A: evaluate gradients
        if (grad) {

            // Evaluate source model
            double spec = (spectral() != NULL) ? spectral()->eval_gradients(srcEng, srcTime) : 1.0;
            double temp = (temporal() != NULL) ? temporal()->eval_gradients(srcTime) : 1.0;

            // Set value
            value = spec * temp * irf;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
                std::cout << "*** ERROR: GModelSky::integrate_dir:";
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
                    for (int i = 0; i < spectral()->size(); ++i) {
                        (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact);
                    }
                }
            }

            // Multiply factors to temporal gradients
            if (temporal() != NULL) {
                double fact = spec * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < temporal()->size(); ++i) {
                        (*temporal())[i].factor_gradient((*temporal())[i].factor_gradient() * fact);
                    }
                }
            }

        } // endif: gradient evaluation has been requested

        // Case B: evaluate no gradients
        else {

            // Evaluate source model
            double spec = (m_spectral != NULL) ? m_spectral->eval(srcEng, srcTime) : 1.0;
            double temp = (m_temporal != NULL) ? m_temporal->eval(srcTime) : 1.0;

            // Set value
            value = spec * temp * irf;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
                std::cout << "*** ERROR: GModelSky::integrate_dir:";
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
 * @brief Integration kernel for edisp_kern() method
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel needed for the edisp_kern()
 * method.
 ***************************************************************************/
double GModelSky::edisp_kern::eval(const double& x)
{
    // Set energy
    GEnergy eng;
    double expx = std::exp(x);
    eng.MeV(expx);

    // Get function value
    double value = m_parent->integrate_dir(m_event, eng, m_srcTime, m_obs, m_grad);

    // Save value if needed
    #if defined(G_NAN_CHECK)
    double value_out = value;
    #endif

    // Correct for variable substitution
    value *= expx;

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSky::edisp_kern::eval";
        std::cout << "(x=" << x << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << " (value_out=" << value_out;
        std::cout << " exp(x)=" << expx;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}
