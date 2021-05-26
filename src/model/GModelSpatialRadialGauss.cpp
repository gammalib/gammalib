/***************************************************************************
 *    GModelSpatialRadialGauss.cpp - Radial Gaussian source model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialRadialGauss.cpp
 * @brief Radial Gaussian model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialRadialGauss.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialGauss g_radial_gauss_seed;
const GModelSpatialRegistry    g_radial_gauss_registry(&g_radial_gauss_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialRadialGauss g_radial_gauss_legacy_seed(true, "GaussFunction");
const GModelSpatialRegistry    g_radial_gauss_legacy_registry(&g_radial_gauss_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_EVAL   "GModelSpatialRadialGauss::eval(double&, GEnergy&, GTime&, "\
                                                                     "bool&)"
#define G_READ                 "GModelSpatialRadialGauss::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialRadialGauss::write(GXmlElement&)"

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
 * Constructs empty radial Gaussian model.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(void) : GModelSpatialRadial()
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
 * Constructs empty radial Gaussian model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const bool&        dummy,
                                                   const std::string& type) :
                          GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] dir Sky position of Gaussian.
 * @param[in] sigma Width of Gaussian (in degrees).
 *
 * Constructs a Gaussian spatial model using a sky direction (@p dir) and
 * a Gaussian width parameter @p sigma in degrees.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const GSkyDir& dir,
                                                   const double&  sigma) :
                          GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Assign direction and sigma
    this->dir(dir);
    this->sigma(sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a Gaussian spatial model by extracting information from an XML
 * element. See the method read() for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const GXmlElement& xml) :
                          GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Radial Gaussian model.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const GModelSpatialRadialGauss& model) : 
                          GModelSpatialRadial(model)
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
 ***************************************************************************/
GModelSpatialRadialGauss::~GModelSpatialRadialGauss(void)
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
 * @param[in] model Radial Gaussian model.
 * @return Radial Gaussian model.
 ***************************************************************************/
GModelSpatialRadialGauss& GModelSpatialRadialGauss::operator=(const GModelSpatialRadialGauss& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialRadial::operator=(model);

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
 * @brief Clear radial Gauss model
 ***************************************************************************/
void GModelSpatialRadialGauss::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatialRadial::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialRadial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone radial Gauss model
 ***************************************************************************/
GModelSpatialRadialGauss* GModelSpatialRadialGauss::clone(void) const
{
    // Clone radial Gauss model
    return new GModelSpatialRadialGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate Gaussian source model
 *
 * @param[in] theta Angular distance from Gaussian centre (radians).
 * @param[in] energy Photon energy (not used).
 * @param[in] time Photon arrival time (not used).
 * @param[in] gradients Compute gradients?
 * @return Model value (\f${\rm sr}^{-1}\f$).
 *
 * Evaluates the spatial component for a Gaussian source model using
 *
 * \f[
 *    M_{\rm S}({\bf p} | E, t) = M_{\rm S}(\theta) =
 *    \frac{1}{2 \pi \sigma^2} \exp 
 *    \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)
 * \f]
 *
 * where
 * - \f$\theta\f$ is the angular separation from the centre of the model, and
 * - \f$\sigma\f$ is the Gaussian width.
 *
 * If @p gradients is `true` the method will also compute analytical parameter
 * gradients
 *
 * \f[
 *    \frac{\partial M_{\rm S}}{\partial \alpha_0} =
 *    \frac{\partial M_{\rm S}}{\partial \theta}
 *    \frac{\partial \theta}{\partial \alpha_0}
 * \f]
 *
 * \f[
 *    \frac{\partial M_{\rm S}}{\partial \beta_0} =
 *    \frac{\partial M_{\rm S}}{\partial \theta}
 *    \frac{\partial \theta}{\partial \beta_0}
 * \f]
 *
 * with
 *
 * \f[
 *    \frac{\partial M_{\rm S}}{\partial \theta} =
 *    -\frac{\theta}{\sigma^2} M_{\rm S}(\theta)
 * \f]
 *
 * and
 *
 * \f[
 *    \frac{\partial M_{\rm S}}{\partial \sigma} =
 *    \left( \frac{\theta^2}{\sigma^3} - \frac{2}{\sigma} \right)
 *    M_{\rm S}(\theta)
 * \f]
 *
 * where
 * - \f$\alpha_0\f$ is the Right Ascension of the model centre, and
 * - \f$\beta_0\f$ is the Declination of the model centre.
 *
 * The computation of the analytical gradients with respect to the model
 * centre relies on the computation of the partial derivatives
 * \f$\frac{\partial \theta}{\partial \alpha_0}\f$ and
 * \f$\frac{\partial \theta}{\partial \beta_0}\f$ which are assumed to
 * be stored in `m_ra.factor_gradient()` and `m_dec.factor_gradient()`
 * when the method is entered. No check is performed whether these values
 * are actually available and valid.
 *
 * Note that the model normalisation is only correct in the small angle
 * approximation.
 ***************************************************************************/
double GModelSpatialRadialGauss::eval(const double&  theta,
                                      const GEnergy& energy,
                                      const GTime&   time,
                                      const bool&    gradients) const
{
    // Update
    update(gradients);

    // Compute value
    double theta2    = theta * theta;
    double arg2      = m_inv_sigma2_rad * theta2;
    double value     = m_value_norm * std::exp(-0.5 * arg2);

    // Optionally compute gradients
    if (gradients) {

        // Evaluate position gradients. Note that this requires that the
        // partial derivatives of theta with respect to Right Ascension
        // and Declination are preset in the m_ra.factor_gradient() and
        // m_dec.factor_gradient() members.
        double g_ra  = 0.0;
        double g_dec = 0.0;
        if (m_ra.is_free() || m_dec.is_free()) {
            double g_theta = -theta * m_g_theta_norm * value;
            if (m_ra.is_free()) {
                g_ra = g_theta * m_ra.factor_gradient() * m_ra.scale();
            }
            if (m_dec.is_free()) {
                g_dec = g_theta * m_dec.factor_gradient() * m_dec.scale();
            }
        }
        m_ra.factor_gradient(g_ra);
        m_dec.factor_gradient(g_dec);

        // Evaluate sigma gradient
        double g_sigma = 0.0;
        if (m_sigma.is_free()) {
            g_sigma = (arg2 - 2.0) * value * m_g_sigma_norm;
        }
        m_sigma.factor_gradient(g_sigma);

    } // endif: gradient computation was requested

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::string msg = "Model value not a number:";
        for (int i = 0; i < m_pars.size(); ++i) {
            msg += " " + m_pars[i]->name() + "=";
            msg += gammalib::str(m_pars[i]->value());
        }
        msg += " energy=" + energy.print();
        msg += " time=" + time.print();
        gammalib::warning(G_EVAL, msg);
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction from the 2D Gaussian distribution as
 * function of the photon @p energy and arrival @p time.
 *
 * @todo This method is only valid in the small angle approximation.
 ***************************************************************************/
GSkyDir GModelSpatialRadialGauss::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Simulate offset from photon arrival direction
    double theta = sigma() * ran.chisq2();
    double phi   = 360.0   * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Checks where model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the radial gauss model.
 ***************************************************************************/
bool GModelSpatialRadialGauss::contains(const GSkyDir& dir,
                                        const double&  margin) const
{
    // Compute distance to centre (radians)
    double distance = dir.dist(this->dir());

    // Return flag
    return (distance <= theta_max() + margin*gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * Returns \f$5 \sigma\f$ as approximate edge of the Gaussian. This limit
 * is of course arbitrary, but allows to limit the integration region for
 * response computation.
 ***************************************************************************/
double GModelSpatialRadialGauss::theta_max(void) const
{
    // Return value
    return (sigma() * gammalib::deg2rad * 5.0);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the radial Gauss model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="RadialGaussian">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="RadialGaussian">
 *       <parameter name="GLON"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"  scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement a test of the sigma and sigma boundary. The sigma
 *       and sigma minimum should be >0.
 ***************************************************************************/
void GModelSpatialRadialGauss::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 3);

    // Read Gaussian location
    GModelSpatialRadial::read(xml);

    // Get parameters
    const GXmlElement* sigma = gammalib::xml_get_par(G_READ, xml, m_sigma.name());

    // Read parameters
    m_sigma.read(*sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Writes the radial disk model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="RadialGaussian">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadialGauss::write(GXmlElement& xml) const
{
    // Write Gaussian location
    GModelSpatialRadial::write(xml);

    // Get or create parameters
    GXmlElement* sigma = gammalib::xml_need_par(G_WRITE, xml, m_sigma.name());

    // Write parameters
    m_sigma.write(*sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Gaussian source information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialGauss ===");

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
 *
 * Note that the minimum Gaussian width is set to 1 arcsec.
 ***************************************************************************/
void GModelSpatialRadialGauss::init_members(void)
{
    // Initialise model type
    m_type = "RadialGaussian";

    // Initialise Gaussian sigma
    m_sigma.clear();
    m_sigma.name("Sigma");
    m_sigma.unit("deg");
    m_sigma.value(2.778e-4); // 1 arcsec
    m_sigma.min(2.778e-4);   // 1 arcsec
    m_sigma.free();
    m_sigma.scale(1.0);
    m_sigma.gradient(0.0);

    // Signal that this model provides analytical parameter gradients
    m_ra.has_grad(true);
    m_dec.has_grad(true);
    m_sigma.has_grad(true);

    // Set parameter pointer(s)
    m_pars.push_back(&m_sigma);

    // Initialise pre-computation cache
    m_last_sigma     = 0.0;
    m_inv_sigma2_rad = 0.0;
    m_value_norm     = 0.0;
    m_g_theta_norm   = 0.0;
    m_g_sigma_norm   = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Gaussian spatial model.
 *
 * We do not have to push back the members on the parameter stack as this
 * should have been done by init_members() that was called before. Otherwise
 * we would have sigma twice on the stack.
 ***************************************************************************/
void GModelSpatialRadialGauss::copy_members(const GModelSpatialRadialGauss& model)
{
    // Copy members
    m_type  = model.m_type;   // Needed to conserve model type
    m_sigma = model.m_sigma;

    // Copy pre-computation cache
    m_last_sigma     = model.m_last_sigma;
    m_inv_sigma2_rad = model.m_inv_sigma2_rad;
    m_value_norm     = model.m_value_norm;
    m_g_theta_norm   = model.m_g_theta_norm;
    m_g_sigma_norm   = model.m_g_sigma_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialGauss::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * @param[in] gradients Gradient computation requested?
 ***************************************************************************/
void GModelSpatialRadialGauss::update(const bool& gradients) const
{
    // Update if radius has changed
    if (m_last_sigma != sigma()) {

        // Store last values
        m_last_sigma = sigma();

        // Compute Gaussian sigma in radians
        double sigma_rad  = sigma() * gammalib::deg2rad;
        double sigma2_rad = sigma_rad * sigma_rad;
        if (sigma2_rad > 0.0) {
            m_inv_sigma2_rad = 1.0 / sigma2_rad;
            m_value_norm     = 1.0 / (gammalib::twopi * sigma2_rad);
        }
        else {
            m_inv_sigma2_rad = 0.0;
            m_value_norm     = 0.0;
        }

        // Pre-compute gradient stuff
        if (gradients) {
            m_g_theta_norm = m_inv_sigma2_rad * gammalib::deg2rad;
            m_g_sigma_norm = m_sigma.scale() * gammalib::deg2rad / sigma_rad;
        }

    } // endif: update required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialRadialGauss::set_region(void) const
{
    // Set sky region circle (5 times Gaussian sigma)
    GSkyRegionCircle region(m_ra.value(), m_dec.value(), m_sigma.value() * 5.0);

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialRadialGauss*>(this)->m_region = region;

    // Return
    return;
}
