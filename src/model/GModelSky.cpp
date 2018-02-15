/***************************************************************************
 *                    GModelSky.cpp - Sky model class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialComposite.hpp"
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
const GModelSky         g_compositesource_seed("CompositeSource");
const GModelRegistry    g_pointsource_registry(&g_pointsource_seed);
const GModelRegistry    g_extendedsource_registry(&g_extendedsource_seed);
const GModelRegistry    g_diffusesource_registry(&g_diffusesource_seed);
const GModelRegistry    g_compositesource_registry(&g_compositesource_seed);

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
 * @brief Set spatial model component
 *
 * @param[in] spatial Pointer to spatial model component.
 *
 * Sets the spatial model component of the model.
 ***************************************************************************/
void GModelSky::spatial(const GModelSpatial* spatial)
{
    // Free spatial model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_spatial != NULL) && (m_spatial != spatial)) {
        delete m_spatial;
    }

    // Clone spatial model component if it exists, otherwise set pointer
    // to NULL
    m_spatial = (spatial != NULL) ? spatial->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Set model type dependent on spatial model type
    set_type();

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
void GModelSky::spectral(const GModelSpectral* spectral)
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
void GModelSky::temporal(const GModelTemporal* temporal)
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
    if (m_spatial  != NULL) m_spatial->eval(photon, true);
    if (m_spectral != NULL) m_spectral->eval(photon.energy(), photon.time(),
                                             true);
    if (m_temporal != NULL) m_temporal->eval(photon.time(), true);

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
 * @param[in] gradients Compute gradients?
 * @return Value of sky model.
 *
 * Evalutes the value of the sky model for an @p event of a specific
 * observation @p obs.
 *
 * If the @p gradients flag is true the method will also compute the
 * parameter gradients for all model parameters.
 ***************************************************************************/
double GModelSky::eval(const GEvent&       event,
                       const GObservation& obs,
                       const bool&         gradients) const
{
    // Evaluate function
    double value = obs.response()->convolve(*this, event, obs, gradients);

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
 * Computes
 *
 * \f[
 *    N_{\rm pred}(E',t') = \int_{\rm ROI}
 *                          P(p',E',t') \, dp' \, dE'
 * \f]
 *
 * of the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * where                         
 * \f$S(p,E,t)\f$ is the source model,
 * \f$R(p',E',t'|p,E,t)\f$ is the instrument response function,
 * \f$p'\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$p\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy, and
 * \f$t\f$ is the true photon arrival time.
 *
 * The method calls the GResponse::nroi() method that does the relevant
 * integration.
 ***************************************************************************/
double GModelSky::npred(const GEnergy&      obsEng,
                        const GTime&        obsTime,
                        const GObservation& obs) const
{
    // Initialise result
    double npred = 0.0;

    // Continue only if model is valid)
    if (valid_model()) {

        // Compute Nroi
        npred = obs.response()->nroi(*this, obsEng, obsTime, obs);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
            std::cout << "*** ERROR: GModelSky::npred:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", obsEng=" << obsEng;
            std::cout << ", obsTime=" << obsTime;
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
 *     <source name=".." type=".." instrument=".." id=".." ts="..">
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
 *
 * Optionally, the model may also contain scale parameters following the
 * format:
 *
 *     <source name=".." type=".." instrument=".." id=".." ts="..">
 *       <spectrum type="..">
 *         ..
 *       </spectrum>
 *       <spatialModel type="..">
 *         ..
 *       </spatialModel>
 *       <temporal type="..">
 *         ..
 *       </temporal>
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

    // Set spatial and spectral models
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
 *     <source name=".." type=".." instrument=".." id=".." ts="..">
 *       <spectrum type="..">
 *         ..
 *       </spectrum>
 *       <spatialModel type="..">
 *         ..
 *       </spatialModel>
 *       <temporal type="..">
 *         ..
 *       </temporal>
 *       <scaling>
 *         <instrument name=".." scale="1.0" min="0.1" max="10.0" value="1.0" free="0"/>
 *         <instrument name=".." scale="1.0" min="0.1" max="10.0" value="0.5" free="0"/>
 *       </scaling>
 *     </source>
 *
 * For compatibility reasons the temporal element will only be written if it
 * is a non-constant component or a constant component with a normalization
 * that differs from unity.
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
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (spatial()  != NULL) src->append(GXmlElement("spatialModel"));
        if (temporal() != NULL && cons == NULL) src->append(GXmlElement("temporal"));
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
        double norm      = m_spatial->mc_norm(dir, radius);
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
            GModelSpatialDiffuseCube* cube =
                         dynamic_cast<GModelSpatialDiffuseCube*>(m_spatial);
            if (cube != NULL) {

                // Set MC cone
                cube->set_mc_cone(dir, radius);

                // Allocate node function to replace the spectral component
                GModelSpectralNodes* nodes =
                                     new GModelSpectralNodes(cube->spectrum());
                for (int i = 0; i < nodes->nodes(); ++i) {
                    GEnergy energy    = nodes->energy(i);
                    double  intensity = nodes->intensity(i);
                    double  value     = m_spectral->eval(energy, tmin);
                    nodes->intensity(i, value*intensity);
                }

                // Signal that node function needs to be de-allocated later
                free_spectral = true;

                // Set the spectral model pointer to the node function
                spectral = nodes;

                // Kluge: if there are no nodes then the spectral->flux method
                // will throw an exception. We therefore set here the use_model
                // flag to false in case that there are no spectral nodes
                if (nodes->nodes() == 0) {
                    use_model = false;
                }

            } // endif: spatial model was a diffuse cube

            // Compute flux within [emin, emax] in model from spectral
            // component (units: ph/cm2/s)
            double flux = (use_model) ? spectral->flux(emin, emax) : 0.0;

            // Derive expecting counting rate within simulation surface
            // (units: ph/s)
            double rate = flux * area * norm;

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GModelSky::mc(\"" << name() << "\": ";
            std::cout << "flux=" << flux << " ph/cm2/s, ";
            std::cout << "rate=" << rate << " ph/s, ";
            std::cout << "norm=" << norm << ")" << std::endl;
            #endif

            // Continue only if rate is positive
            if (rate > 0.0) {

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

                    // Set photon energy. If an invalid_return_value exception
                    // occurs the energy returned by the spectral Monte Carlo
                    // method is invalid and the photon is skipped.
                    try {
                        photon.energy(spectral->mc(emin, emax, photon.time(),
                                                   ran));
                    }
                    catch (GException::invalid_return_value) {
                        continue;
                    }

                    // Debug option: dump energy
                    #if defined(G_DUMP_MC_DETAIL)
                    std::cout << "    Energy=" << photon.energy() << std::endl;
                    #endif

                    // Set incident photon direction. If an invalid_return_value
                    // exception occurs the sky direction returned by the
                    // spatial Monte Carlo method is invalid and the photon is
                    // skipped.
                    try {
                        photon.dir(m_spatial->mc(photon.energy(), photon.time(),
                                                 ran));
                    }
                    catch (GException::invalid_return_value) {
                        continue;
                    }

                    // Debug option: dump direction
                    #if defined(G_DUMP_MC_DETAIL)
                    std::cout << "    Direction=" << photon.dir() << std::endl;
                    #endif

                    // Append photon
                    if (dir.dist_deg(photon.dir()) <= radius) {
                        photons.append(photon);
                    }

                } // endfor: looped over photons

            } // endif: rate was positive

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
 * @param[in] chatter Chattiness.
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

        // ... otherwise, is spatial model a composite source?
        else if (dynamic_cast<const GModelSpatialComposite*>(m_spatial) != NULL) {
            m_type = "CompositeSource";
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
 * @brief Return pointer to spatial model from XML element
 *
 * @param[in] spatial XML element.
 * @return Pointer to spatial model.
 *
 * Returns pointer to spatial model that is defined in an XML element.
 ***************************************************************************/
GModelSpatial* GModelSky::xml_spatial(const GXmlElement& spatial) const
{
    // Get spatial model
    GModelSpatialRegistry registry;
    GModelSpatial*        ptr = registry.alloc(spatial);

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
GModelSpectral* GModelSky::xml_spectral(const GXmlElement& spectral) const
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
GModelTemporal* GModelSky::xml_temporal(const GXmlElement& temporal) const
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
 * @return True is model is valid, false otherwise.
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
