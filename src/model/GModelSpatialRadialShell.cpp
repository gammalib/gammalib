/***************************************************************************
 *      GModelSpatialRadialShell.cpp - Radial shell source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Christoph Deil                              *
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
 * @file GModelSpatialRadialShell.cpp
 * @brief Radial shell model class implementation
 * @author Christoph Deil
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialShell g_radial_shell_seed;
const GModelSpatialRegistry    g_radial_shell_registry(&g_radial_shell_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialRadialShell g_radial_shell_legacy_seed(true, "ShellFunction");
const GModelSpatialRegistry    g_radial_shell_legacy_registry(&g_radial_shell_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_MC          "GModelSpatialRadialShell::mc(GEnergy&, GTime&, GRan&)"
#define G_READ                 "GModelSpatialRadialShell::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialRadialShell::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_MC                                     //!< Debug MC method


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty radial shell model.
 ***************************************************************************/
GModelSpatialRadialShell::GModelSpatialRadialShell(void) : GModelSpatialRadial()
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
 * Constructs empty radial shell model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialRadialShell::GModelSpatialRadialShell(const bool&        dummy,
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
 * @brief Shell constructor
 *
 * @param[in] dir Sky position of shell centre.
 * @param[in] radius Inner shell radius (degrees).
 * @param[in] width Shell width (degrees).
 *
 * Constructs the shell model from the shell centre (@p dir), the inner
 * shell @p radius, and the shell @p width.
 ***************************************************************************/
GModelSpatialRadialShell::GModelSpatialRadialShell(const GSkyDir& dir,
                                                   const double&  radius,
                                                   const double&  width) :
                          GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Assign parameters
    this->dir(dir);
    this->radius(radius);
    this->width(width);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs shell spatial model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialShell::GModelSpatialRadialShell(const GXmlElement& xml) :
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
 * @param[in] model Radial shell source model.
 ***************************************************************************/
GModelSpatialRadialShell::GModelSpatialRadialShell(const GModelSpatialRadialShell& model) :
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
GModelSpatialRadialShell::~GModelSpatialRadialShell(void)
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
 * @param[in] model Radial shell source model.
 * @return Radial shell source model.
 ***************************************************************************/
GModelSpatialRadialShell& GModelSpatialRadialShell::operator=(const GModelSpatialRadialShell& model)
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
 * @brief Clear radial shell model
 ***************************************************************************/
void GModelSpatialRadialShell::clear(void)
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
 * @brief Clone radial shell model
 *
 * @return Pointer to deep copy of radial shell model.
 ***************************************************************************/
GModelSpatialRadialShell* GModelSpatialRadialShell::clone(void) const
{
    // Clone radial shell model
    return new GModelSpatialRadialShell(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from shell centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * Evaluates the spatial part for a shell source model. The shell source
 * model is a radial function \f$f(\theta)\f$, where \f$\theta\f$ is the
 * angular separation between shell centre and the actual location.
 *
 * The function is given by
 * \f[
 * S_{\rm p}(\vec{p} | E, t) = {\tt m\_norm} \left \{
 *  \begin{array}{l l}
 *     \displaystyle
 *     \sqrt{ \sin^2 \theta_{\rm out} - \sin^2 \theta } -
 *     \sqrt{ \sin^2 \theta_{\rm in}  - \sin^2 \theta }
 *     & \mbox{if $\theta \le \theta_{\rm in}$} \\
 *     \\
 *    \displaystyle
 *     \sqrt{ \sin^2 \theta_{\rm out} - \sin^2 \theta }
 *     & \mbox{if $\theta_{\rm in} < \theta \le \theta_{\rm out}$} \\
 *     \\
 *    \displaystyle
 *    0 & \mbox{if $\theta > \theta_{\rm out}$}
 *  \end{array}
 *  \right .
 * \f]
 * 
 * Here, \f$\theta_{\rm in}\f$ and \f$\theta_{\rm out}\f$ are the shell inner
 * and outer radius.
 *
 * The method will not compute analytical parameter gradients, even if the
 * @p gradients argument is set to true. Radial disk parameter gradients
 * need to be computed numerically.
 ***************************************************************************/
double GModelSpatialRadialShell::eval(const double&  theta,
                                      const GEnergy& energy,
                                      const GTime&   time,
                                      const bool&    gradients) const
{
    // Update precomputation cache
    update();

    // Set x appropriately
    double x = std::sin(theta);
    x *= x;

    // Compute result
    double result = 0.0;
    if (x < m_x_out) {
        result = std::sqrt(m_x_out - x);
        if (x < m_x_in) {
            result -= std::sqrt(m_x_in - x);
        }
        result *= m_norm;
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(result) || gammalib::is_infinite(result)) {
        std::cout << "*** ERROR: GModelSpatialRadialShell::eval";
        std::cout << "(theta=" << theta << "): NaN/Inf encountered";
        std::cout << " (result=" << result;
        std::cout << ", m_norm=" << m_norm;
        std::cout << ", x=" << x;
        std::cout << ", m_x_out=" << m_x_out;
        std::cout << ", m_x_in=" << m_x_in;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return normalised value
    return result;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * @exception GException::invalid_return_value
 *            Invalid u_max value
 *
 * Draws an arbitrary sky position from the 2D shell distribution.
 ***************************************************************************/
GSkyDir GModelSpatialRadialShell::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Update precomputation cache
    update();

    // Initialise simulation
    #if defined(G_DEBUG_MC)
    int    n_samples = 0;
    #endif
    double u_max     = m_norm * m_x_out;
    double value     = 0.0;
    double u         = 1.0;
    double theta     = 0.0;

    // Throw an exception if the maximum value is zero
    if (u_max <= 0.0) {
        std::string msg = "Non positive maximum radial shell model value.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Simulate offset from photon arrival direction using a rejection method
    do {
        theta = ran.uniform() * m_theta_out;
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
 * Signals whether a sky direction is contained in the radial shell model.
 ***************************************************************************/
bool GModelSpatialRadialShell::contains(const GSkyDir& dir,
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
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialShell::theta_max(void) const
{
    // Return value
    return ((radius()+width()) * gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Reads the radial shell model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="RadialShell">
 *       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.30"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Width"  scale="1.0" value="0.10"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="RadialShell">
 *       <parameter name="GLON"   scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.30"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Width"  scale="1.0" value="0.10"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement tests of the radius and radius boundary and the width and
 *       width boundary. The radius and radius boundary should be >=0, the
 *       width and width boundary should be >0.
 ***************************************************************************/
void GModelSpatialRadialShell::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || npars != 4) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Shell model requires exactly 4 parameters.");
    }

    // Read shell location
    GModelSpatialRadial::read(xml);

    // Get parameters
    const GXmlElement* radius = gammalib::xml_get_par(G_READ, xml, m_radius.name());
    const GXmlElement* width  = gammalib::xml_get_par(G_READ, xml, m_width.name());

    // Read parameters
    m_radius.read(*radius);
    m_width.read(*width);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type 'GaussFunction'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the radial shell model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="RadialShell">
 *       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.30"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Width"  scale="1.0" value="0.10"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadialShell::write(GXmlElement& xml) const
{
    // Write shell location
    GModelSpatialRadial::write(xml);

    // Get or create parameters
    GXmlElement* radius = gammalib::xml_need_par(G_WRITE, xml, m_radius.name());
    GXmlElement* width  = gammalib::xml_need_par(G_WRITE, xml, m_width.name());

    // Write parameters
    m_radius.write(*radius);
    m_width.write(*width);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialShell::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialShell ===");

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
void GModelSpatialRadialShell::init_members(void)
{
    // Initialise model type
    m_type = "RadialShell";

    // Initialise Radius
    m_radius.clear();
    m_radius.name("Radius");
    m_radius.unit("deg");
    m_radius.value(0.0);
    m_radius.min(0.0);
    m_radius.free();
    m_radius.scale(1.0);
    m_radius.gradient(0.0);
    m_radius.has_grad(false);  // Radial components never have gradients

    // Initialise Width
    m_width.clear();
    m_width.name("Width");
    m_width.unit("deg");
    m_width.value(2.778e-4);
    m_width.min(2.778e-4);   // 1 arcsec
    m_width.free();
    m_width.scale(1.0);
    m_width.gradient(0.0);
    m_width.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_radius);
    m_pars.push_back(&m_width);

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius and width shell is not meaningful
    m_last_radius = 0.0;
    m_last_width  = 0.0;
    m_theta_in    = 0.0;
    m_x_in        = 0.0;
    m_theta_out   = 0.0;
    m_x_out       = 0.0;
    m_norm        = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial shell source model.
 *
 * We do not have to push back the members on the parameter stack as this
 * should have been done by init_members() that was called before. Otherwise
 * we would have the radius and width twice on the stack.
 ***************************************************************************/
void GModelSpatialRadialShell::copy_members(const GModelSpatialRadialShell& model)
{
    // Copy members
    m_type   = model.m_type;   // Needed to conserve model type
    m_radius = model.m_radius;
    m_width  = model.m_width;

    // Copy precomputation cache
    m_last_radius = model.m_last_radius;
    m_last_width  = model.m_last_width;
    m_theta_in    = model.m_theta_in;
    m_x_in        = model.m_x_in;
    m_theta_out   = model.m_theta_out;
    m_x_out       = model.m_x_out;
    m_norm        = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialShell::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Performs precomputations that are valid for a given pair of radius and
 * width values. The following members are set by this method:
 * m_last_radius, m_last_width, m_theta_in, m_theta_out, m_x_in, m_x_out
 * and m_norm.
 *
 * m_theta_in contains the inner shell radius in radians, while m_theta_out
 * contains the outer shell radius in radians.
 *
 * \f${\tt m\_x\_in} = \sin^2 {\tt m\_theta\_in}\f$,
 * \f${\tt m\_x\_out} = \sin^2 {\tt m\_theta\_out}\f$, and
 * \f[{\tt m\_norm} = \frac{1}{2 \pi}
 *    \frac{\sqrt{1-\cos 2 {\tt m\_theta\_out}} -
 *          \sqrt{1-\cos 2 {\tt m\_theta\_in}}}{2 \sqrt{2}} +
 *    \frac{1+\cos 2 {\tt m\_theta\_out}}{4} \ln \left(
 *          \frac{\sqrt{2} \cos {\tt m\_theta\_out}}
 *               {\sqrt{2} + \sqrt{1 - \cos 2 {\tt m\_theta\_out}}} \right) -
 *    \frac{1+\cos 2 {\tt m\_theta\_in}}{4} \ln \left(
 *          \frac{\sqrt{2} \cos {\tt m\_theta\_in}}
 *               {\sqrt{2} + \sqrt{1 - \cos 2 {\tt m\_theta\_in}}} \right)\f]
 ***************************************************************************/
void GModelSpatialRadialShell::update() const
{
    // Set constants
    //const double c1 = gammalib::twopi / 3.0;
    const double c2 = 1.0 / (2.0 * gammalib::sqrt_two);

    // Update if radius or width have changed
    if (m_last_radius != radius() || m_last_width != width()) {

        // Store last values
        m_last_radius = radius();
        m_last_width  = width();

        // Perform precomputations
        m_theta_in           = radius() * gammalib::deg2rad;
        m_theta_out          = (radius() + width()) * gammalib::deg2rad;
        double sin_theta_in  = std::sin(m_theta_in);
        double sin_theta_out = std::sin(m_theta_out);
        double term1         = (f1(m_theta_out) - f1(m_theta_in)) * c2;
        double term2         = f2(m_theta_out);
        double term3         = f2(m_theta_in);
        double denom         = gammalib::twopi * (term1 + term2 - term3);
        m_norm               = (denom > 0.0) ? 1.0 / denom : 0.0;
        m_x_in               = sin_theta_in  * sin_theta_in;
        m_x_out              = sin_theta_out * sin_theta_out;

    } // endif: update required

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(m_norm) || gammalib::is_infinite(m_norm)) {
        std::cout << "*** ERROR: GModelSpatialRadialShell::update:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (m_norm=" << m_norm;
        std::cout << ", radius=" << radius();
        std::cout << ", width=" << width();
        std::cout << ", m_theta_in=" << m_theta_in;
        std::cout << ", m_theta_out=" << m_theta_out;
        std::cout << ", term1=" << (f1(m_theta_out) - f1(m_theta_in)) * c2;
        std::cout << ", term2=" << f2(m_theta_out);
        std::cout << ", term3=" << f2(m_theta_in);
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return function 1 value needed for precomputation
 *
 * Computes \f$f1(x) = \sqrt{1 - \cos 2 x}\f$.
 ***************************************************************************/
double GModelSpatialRadialShell::f1(double x)
{
    // Compute value
    double f1 = std::sqrt(1.0 - std::cos(2.0 * x));

    // Return value
    return f1;
}


/***********************************************************************//**
 * @brief Return function 2 value needed for precomputation
 *
 * Compute
 * \f[f2(x) = \frac{1+\cos 2x}{4} 
 *    \ln \left( \frac{\sqrt{2} \cos x}{\sqrt{2} + \sqrt{ 1 - \cos 2 x}}
 *        \right)\f].
 ***************************************************************************/
double GModelSpatialRadialShell::f2(double x)
{
    // Compute value
    double t1 = (1.0 + std::cos(2.0*x)) / 4.0;
    double t2 = gammalib::sqrt_two * std::cos(x);
    double t3 = gammalib::sqrt_two + f1(x);
    double f2 = t1 * std::log(t2 / t3);

    // Return value
    return f2;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialRadialShell::set_region(void) const
{
    // Set sky region circle (maximum Gaussian sigma times a scaling
    // factor (actually 3))
    GSkyRegionCircle region(m_ra.value(), m_dec.value(), m_radius.value() + m_width.value());

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialRadialShell*>(this)->m_region = region;

    // Return
    return;
}
