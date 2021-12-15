/***************************************************************************
 *    GModelSpatialRadialGeneralGauss.cpp - Radial Generalized Gaussian source model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021- by Luigi Tibaldo                                   *
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
 * @file GModelSpatialRadialGeneralGauss.cpp
 * @brief Radial Generallized Gaussian model class implementation
 * @author Luigi Tibaldo
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialRadialGeneralGauss.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialGeneralGauss g_radial_general_gauss_seed;
const GModelSpatialRegistry    g_radial_general_gauss_registry(&g_radial_general_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL   "GModelSpatialRadialGeneralGauss::eval(double&, GEnergy&, GTime&, "\
                                                                     "bool&)"
#define G_READ                 "GModelSpatialRadialGeneralGauss::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialRadialGeneralGauss::write(GXmlElement&)"
#define G_MC          "GModelSpatialRadialGeneralGauss::mc(GEnergy&, GTime&, GRan&)"

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
GModelSpatialRadialGeneralGauss::GModelSpatialRadialGeneralGauss(void) : GModelSpatialRadial()
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
GModelSpatialRadialGeneralGauss::GModelSpatialRadialGeneralGauss(const bool&        dummy,
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
 * @param[in] radius Width of Generalized Gaussian (in degrees).
 * @param[in] ridx Reciprocal of exponential index of radial profile.
 *
 * Constructs a Generalized Gaussian spatial model using a sky direction (@p dir),
 * a Gaussian width parameter @p radius in degrees and a reciprocal of the 
 * exponential index @p ridx.
 ***************************************************************************/
GModelSpatialRadialGeneralGauss::GModelSpatialRadialGeneralGauss(const GSkyDir& dir,
                                                   const double&  radius,
						   const double& ridx) :
                          GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Assign direction and radius
    this->dir(dir);
    this->radius(radius);
    this->ridx(ridx);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a Generalized Gaussian spatial model by extracting information from an XML
 * element. See the method read() for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialGeneralGauss::GModelSpatialRadialGeneralGauss(const GXmlElement& xml) :
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
GModelSpatialRadialGeneralGauss::GModelSpatialRadialGeneralGauss(const GModelSpatialRadialGeneralGauss& model) : 
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
GModelSpatialRadialGeneralGauss::~GModelSpatialRadialGeneralGauss(void)
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
GModelSpatialRadialGeneralGauss& GModelSpatialRadialGeneralGauss::operator=(const GModelSpatialRadialGeneralGauss& model)
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
void GModelSpatialRadialGeneralGauss::clear(void)
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
GModelSpatialRadialGeneralGauss* GModelSpatialRadialGeneralGauss::clone(void) const
{
    // Clone radial Gauss model
    return new GModelSpatialRadialGeneralGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate Generalized Gaussian source model
 *
 * @param[in] theta Angular distance from source centre (radians).
 * @param[in] energy Photon energy (not used).
 * @param[in] time Photon arrival time (not used).
 * @param[in] gradients Compute gradients?
 * @return Model value (\f${\rm sr}^{-1}\f$).
 *
 * Evaluates the spatial component for a Generalized Gaussian source model using
 *
 * \f[
 *    M_{\rm S}({\bf p} | E, t) = M_{\rm S}(\theta) =
 *    \frac{1}{2 \pi r^2 \eta \Gamma(\eta)} \exp 
 *    \left(-(\frac{\theta}{r})^1/\eta \right)
 * \f]
 *
 * where
 * - \f$\theta\f$ is the angular separation from the centre of the model, and
 * - \f$r\f$ is the Generalized Gaussian radius.
 * - \f$\eta\f$ is the reciprocal of the radial profile exponent.
 *
 * If @p gradients is `true` the method will also compute parameter
 * gradients
 *
 *
 * Note that the model normalisation is only correct in the small angle
 * approximation and for  \f$\eta\f$ order unity or smaller.
 ***************************************************************************/
double GModelSpatialRadialGeneralGauss::eval(const double&  theta,
                                      const GEnergy& energy,
                                      const GTime&   time,
                                      const bool&    gradients) const
{
    // Update
    update();

    // Compute value
    double arg       = -1.0 * std::pow(m_inv_radius_rad * theta,m_inv_ridx);
    double value     = m_value_norm * std::exp(arg);

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
 * Draws an arbitrary sky direction from the 2D Generalized  Gaussian distribution
 * as function of the photon @p energy and arrival @p time.
 *
 * @todo This method is only valid in the small angle approximation.
 ***************************************************************************/
GSkyDir GModelSpatialRadialGeneralGauss::mc(const GEnergy& energy,
                                            const GTime&   time,
                                            GRan&          ran) const
{

    // Update precomputation cache
    update();

    // Initialise simulation
    #if defined(G_DEBUG_MC)
    int    n_samples = 0;
    #endif
    double u_max     = m_value_norm * std::sin(theta_max());
    double value     = 0.0;
    double u         = 1.0;
    double theta     = 0.0;

    // Throw an exception if the maximum value is zero
    if (u_max <= 0.0) {
        std::string msg = "Non positive maximum radial genral Gauss model value.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Simulate offset from photon arrival direction using a rejection method
    do {
        theta = ran.uniform() * theta_max();
        value = eval(theta, energy, time) * std::sin(theta);
        u     = ran.uniform() * u_max;
        #if defined(G_DEBUG_MC)
        n_samples++;
        #endif
    } while (u > value);
    #if defined(G_DEBUG_MC)
    std::cout << "#=" << n_samples << " ";
    #endif

    // Simulate azimuth angle
    double phi = 360.0 * ran.uniform();

    // Rotate pointing direction by offset and azimuth angle
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta * gammalib::rad2deg);

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
bool GModelSpatialRadialGeneralGauss::contains(const GSkyDir& dir,
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
double GModelSpatialRadialGeneralGauss::theta_max(void) const
{
    // Return value
  return (radius() * gammalib::deg2rad * 5.0);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the radial Generalized Gauss model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="RadialGaussian">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
*       <parameter name="R_Index" scale="1.0" value="0.5"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="RadialGaussian">
 *       <parameter name="GLON"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"  scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
*       <parameter name="R_Index" scale="1.0" value="0.5"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement a test of the radius and radius boundary. The sigma
 *       and sigma minimum should be >0.
 ***************************************************************************/
void GModelSpatialRadialGeneralGauss::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 4);

    // Read location
    GModelSpatialRadial::read(xml);

    // Get parameters
    const GXmlElement* radius = gammalib::xml_get_par(G_READ, xml, m_radius.name());
    const GXmlElement* ridx   = gammalib::xml_get_par(G_READ, xml, m_ridx.name());

    // Read parameters
    m_radius.read(*radius);
    m_ridx.read(*ridx);

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
 *       <parameter name="R_Index" scale="1.0" value="0.5"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadialGeneralGauss::write(GXmlElement& xml) const
{
    // Write Gaussian location
    GModelSpatialRadial::write(xml);

    // Get or create parameters
    GXmlElement* radius = gammalib::xml_need_par(G_WRITE, xml, m_radius.name());
    GXmlElement* ridx   = gammalib::xml_need_par(G_WRITE, xml, m_ridx.name());

    // Write parameters
    m_radius.write(*radius);
    m_ridx.write(*ridx);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Generalized Gaussian source information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialGeneralGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialGeneralGauss ===");

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
 * Note that the minimum Gaussian width is set to 1 arcsec
 * and w set the minimum reciprocal index to 1.e-15
 ***************************************************************************/
void GModelSpatialRadialGeneralGauss::init_members(void)
{
    // Initialise model type
    m_type = "RadialGeneralGaussian";

    // Initialise Radius
    m_radius.clear();
    m_radius.name("Radius");
    m_radius.unit("deg");
    m_radius.value(2.778e-4); // 1 arcsec
    m_radius.min(2.778e-4); // 1 arcsec
    m_radius.free();
    m_radius.scale(1.0);
    m_radius.gradient(0.0);
    m_radius.has_grad(false);

    // Initialise Reciprocal index
    m_ridx.clear();
    m_ridx.name("R_Index");
    m_ridx.value(0.5); // Gaussian
    m_ridx.min(1.e-15);   // need to be > 0
    m_ridx.free();
    m_ridx.scale(1.0);
    m_ridx.gradient(0.0);
    m_ridx.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_radius);
    m_pars.push_back(&m_ridx);

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius and width shell is not meaningful
    m_last_radius = 0.0;
    m_inv_radius_rad = 0.0;
    m_last_ridx = 0.0;
    m_inv_ridx = 0.0;
    m_value_norm = 0.0;

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
void GModelSpatialRadialGeneralGauss::copy_members(const GModelSpatialRadialGeneralGauss& model)
{
    // Copy members
    m_type  = model.m_type;   // Needed to conserve model type
    m_radius = model.m_radius;
    m_ridx = model.m_ridx;

    // Copy pre-computation cache
    m_last_radius     = model.m_last_radius;
    m_inv_radius_rad = model.m_inv_radius_rad;
    m_last_ridx      = model.m_last_ridx;
    m_inv_ridx       = model.m_inv_ridx;
    m_value_norm     = model.m_value_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialGeneralGauss::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * @param[in] gradients Gradient computation requested?
 ***************************************************************************/
void GModelSpatialRadialGeneralGauss::update() const
{
    // Update if radius has changed
  if (m_last_radius != radius() or m_last_ridx != ridx()) {

        // Store last values
        m_last_radius = radius();
	m_last_ridx = ridx();

        // Compute radius in radians
        double radius_rad  = radius() * gammalib::deg2rad;
        if (radius_rad > 0.0 and ridx() > 0) {
            m_inv_radius_rad = 1.0 / radius_rad;
	    m_inv_ridx = 1.0 / ridx() ;
            m_value_norm     = 1.0 / (gammalib::twopi * radius_rad * radius_rad * ridx() * std::exp(gammalib::gammln(2 * ridx())) );
        }
        else {
            m_inv_radius_rad = 0.0;
	    m_inv_ridx = 1.0;
            m_value_norm     = 0.0;
        }

    } // endif: update required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialRadialGeneralGauss::set_region(void) const
{
    // Set sky region circle (5 times radius)
    GSkyRegionCircle region(m_ra.value(), m_dec.value(), m_radius.value() * 5.0);

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialRadialGeneralGauss*>(this)->m_region = region;

    // Return
    return;
}
