/***************************************************************************
 *       GModelSpatialDiffuseConst.cpp - Spatial isotropic model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseConst.cpp
 * @brief Isotropic spatial model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialDiffuseConst.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDiffuseConst g_spatial_const_seed;
const GModelSpatialRegistry     g_spatial_const_registry(&g_spatial_const_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialDiffuseConst g_spatial_const_legacy_seed(true, "ConstantValue");
const GModelSpatialRegistry     g_spatial_const_legacy_registry(&g_spatial_const_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_MC_NORM     "GModelSpatialDiffuseConst::mc_norm(GSkyDir&, double&)"
#define G_READ                "GModelSpatialDiffuseConst::read(GXmlElement&)"
#define G_WRITE              "GModelSpatialDiffuseConst::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(void) :
                           GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model type constructor
 *
 * @param[in] dummy Dummy flag.
 * @param[in] type Model type.
 *
 * Constructs empty isotropic spatial model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const bool&        dummy,
                                                     const std::string& type) :
                           GModelSpatialDiffuse()
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
 * Constructs isotropic spatial model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const GXmlElement& xml) :
                           GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Value constructor
 *
 * @param[in] value Isotropic intensity value (\f$sr^{-1}\f$).
 *
 * Constructs isotropic spatial model by assigning an intensity of @p value
 * in units of \f$sr^{-1}\f$ to the diffuse emission. This constructor
 * explicitly sets the @a m_value member of the model.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const double& value) :
                           GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set parameter
    m_value.value(value);

    // Perform autoscaling of parameter
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Isotropic spatial model.
 *
 * Constructs isotropic spatial model by copying another isotropic spatial
 * model.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const GModelSpatialDiffuseConst& model) :
                           GModelSpatialDiffuse(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destructs isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst::~GModelSpatialDiffuseConst(void)
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
 * @param[in] model Isotropic spatial model.
 * @return Isotropic spatial model.
 *
 * Assigns an isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst& GModelSpatialDiffuseConst::operator=(const GModelSpatialDiffuseConst& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialDiffuse::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
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
 * @brief Clear isotropic spatial model
 *
 * Clears the isotropic spatial model. This method is equivalent to creating
 * an isotropic spatial model using the void constructor.
 ***************************************************************************/
void GModelSpatialDiffuseConst::clear(void)
{
    // Free class members
    free_members();
    this->GModelSpatialDiffuse::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialDiffuse::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone isotropic spatial model
 *
 * @return Pointer to deep copy of isotropic spatial model.
 *
 * Returns a pointer to a deep copy of the isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst* GModelSpatialDiffuseConst::clone(void) const
{
    // Clone isotropic spatial model
    return new GModelSpatialDiffuseConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate isotropic spatial model value
 *
 * @param[in] photon Incident photon.
 * @param[in] gradients Compute gradients?
 * @return Model value (\f$sr^{-1}\f$).
 *
 * Evaluates the isotropic spatial model for a given @p photon, characterised
 * by a sky direction, energy and arrival time \f$(\vec{p},E,t)\f$. By
 * definition, the isotropic spatial model and gradient are independent of
 * sky direction, energy and time. The model value is given by
 *
 * \f[
 *    M_\mathrm{S}(\vec{p}|E,t) = \frac{v}{4\pi}
 * \f]
 *
 * where \f$v\f$ is the value() parameter, which is divided by the solid
 * angle \f$4\pi\f$ of the celestial sphere. 
 *
 * If the @p gradients flag is true the method will also evaluate the partial
 * derivatives of the model. The value() gradient is given by
 *
 * \f[
 *    \frac{\partial M_\mathrm{S}(\vec{p}|E,t)}{\partial v_v} =
 *    \frac{v_s}{4\pi}
 * \f]
 *
 * where \f$v=v_v v_s\f$ is the factorisation of the value() parameter into
 * a value \f$v_v\f$ and scale \f$v_s\f$ term.
 ***************************************************************************/
double GModelSpatialDiffuseConst::eval(const GPhoton& photon,
                                       const bool&    gradients) const
{
    // Set normalization constant
    const double norm = 1.0 / gammalib::fourpi;

    // Compute model value
    double value = norm * m_value.value();

    // Optionally compute partial derivatives
    if (gradients) {

        // Compute partial derivatives of the parameter values
        double g_norm = (m_value.is_free()) ? norm * m_value.scale() : 0.0;

        // Set gradient
        m_value.factor_gradient(g_norm);

    } // endif: computed partial derivatives

    // Return model value
    return (value);
}


/***********************************************************************//**
 * @brief Return MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Returns an arbitrary direction on the celestial sphere within a simulation
 * cone region that has been defined by a previous call to the mc_norm()
 * method (if no such call was made the entire sky will be assumed as the
 * simulation cone). The sky direction is independent of event @p energy and
 * arrival @p time.
 ***************************************************************************/
GSkyDir GModelSpatialDiffuseConst::mc(const GEnergy& energy,
                                      const GTime&   time,
                                      GRan&          ran) const
{
    // Simulate offset from simulation cone centre
    double theta = std::acos(1.0 - ran.uniform() * (1.0 - m_mc_cos_radius)) *
                   gammalib::rad2deg;
    double phi   = 360.0 * ran.uniform();

    // Initialise sky direction with simulation cone centre
    GSkyDir dir(m_mc_centre);

    // Rotate sky direction by offset angles (phi, theta)
    dir.rotate_deg(phi, theta);

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Return normalization of constant diffuse source for Monte Carlo
 *        simulations
 *
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (deg).
 * @return Model normalization.
 *
 * @exception GException::invalid_argument
 *            Radius of simulation cone not coprised in [0,180] degrees.
 *
 * Returns the normalization of a constant diffuse source. The normalization
 * is given by the product between the fraction \f$f\f$ of the sky covered
 * by the simulation cone 
 *
 * \f[
 *    f = \frac{2\pi \left( 1-\cos(r) \right)}{4\pi}
 * \f]
 *
 * (where \f$r\f$ is the radius of the simulation cone) and the model
 * parameter value(). The normalization only depends of the radius of the
 * simulation cone and is invariant to its centre.
 ***************************************************************************/
double GModelSpatialDiffuseConst::mc_norm(const GSkyDir& dir,
                                          const double&  radius) const
{
    // Throw exception if the radius is invalid
    if (radius < 0.0 || radius > 180.0) {
        std::string msg = "Invalid simulation cone radius "+
                          gammalib::str(radius)+" specified. Please provide "
                          "a value comprised within 0 and 180 degrees.";
        throw GException::invalid_argument(G_MC_NORM, msg);
    }

    // Store cosine of simulation cone radius and cone centre
    m_mc_centre     = dir;
    m_mc_cos_radius = std::cos(radius * gammalib::deg2rad);

    // Compute fraction of sky covered by the simulation cone
    double fraction = 0.5 * (1.0 - m_mc_cos_radius);

    // Multiply by normalization value of model
    double norm = fraction * value();

    // Return normalization
    return (norm);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the isotropic source model information from an XML element. The XML
 * element is expected to have the following format:
 *
 *     <spatialModel type="DiffuseIsotropic">
 *       <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialDiffuseConst::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 1);

    // Get value parameter
    const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, m_value.name());

    // Read value parameter
    m_value.read(*par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write the isotropic source model information into an XML element. The XML
 * element will have the following format:
 *
 *     <spatialModel type="DiffuseIsotropic">
 *       <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialDiffuseConst::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get or create Value parameter
    GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, m_value.name());

    // Write parameter
    m_value.write(*par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns isotropic flux integrated in sky region
 *
 * @param[in] region Sky region.
 * @param[in] srcEng Energy.
 * @param[in] srcTime Time.
 * @return Flux (adimensional or ph/cm2/s).
 *
 * Returns isotropic flux within a sky region. The flux \f$F\f$ is computed
 * using
 *
 * \f[F = \frac{{\tt value}}{4 \pi} \times \Omega\f]
 *
 * where
 * - \f${\tt value}\f$ is the normalisation factor returned by the value()
 *   method, and
 * - \f$\Omega\f$ is the solid angle of the sky region.
 ***************************************************************************/
double GModelSpatialDiffuseConst::flux(const GSkyRegion& region,
                                       const GEnergy&    srcEng,
                                       const GTime&      srcTime) const
{
    // Set normalization constant
    const double norm = 1.0 / gammalib::fourpi;

    // Compute flux in sky region
    double flux = norm * m_value.value() * region.solidangle();

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Print isotropic source model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialDiffuseConst::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialDiffuseConst ===");

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
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
void GModelSpatialDiffuseConst::init_members(void)
{
    // Initialise model type
    m_type = "DiffuseIsotropic";

    // Initialise Value
    m_value.clear();
    m_value.name("Value");
    m_value.fix();
    m_value.value(1.0);
    m_value.scale(1.0);
    m_value.range(0.0, 10.0);
    m_value.gradient(0.0);
    m_value.has_grad(true);

    // Initialise other members
    m_mc_centre.clear();
    m_mc_cos_radius = -1.0; // Set simulation cone to allsky

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Isotropic spatial model.
 ***************************************************************************/
void GModelSpatialDiffuseConst::copy_members(const GModelSpatialDiffuseConst& model)
{
    // Copy members
    m_type          = model.m_type;   // Needed to conserve model type
    m_value         = model.m_value;
    m_mc_centre     = model.m_mc_centre;
    m_mc_cos_radius = model.m_mc_cos_radius;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialDiffuseConst::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialDiffuseConst::set_region(void) const
{
    // Set sky region circle (all sky)
    GSkyRegionCircle region(0.0, 0.0, 180.0);

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialDiffuseConst*>(this)->m_region = region;

    // Return
    return;
}
