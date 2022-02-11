/***************************************************************************
 *    GModelSpatialRadialGeneralGauss.cpp - Generalised radial Gaussian    *
 *                           source model class                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Luigi Tibaldo                                    *
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
 * @file GModelSpatialEllipticalGeneralGauss.cpp
 * @brief Generalised elliptical Gaussian model class implementation
 * @author Luigi Tibaldo
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialEllipticalGeneralGauss.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */
namespace {
    const double c_theta_max = 5.0; //!< max distance from center in units of major semiaxis
}

/* __ Globals ____________________________________________________________ */
const GModelSpatialEllipticalGeneralGauss g_elliptical_general_gauss_seed;
const GModelSpatialRegistry               g_elliptical_general_gauss_registry(&g_elliptical_general_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR                 "GModelSpatialEllipticalGeneralGauss::"\
           "GModelSpatialEllipticalGeneralGauss(GSkyDir&, double&, double&, "\
                                            "double&, double&, std::string&)"
#define G_READ      "GModelSpatialEllipticalGeneralGauss::read(GXmlElement&)"
#define G_WRITE    "GModelSpatialEllipticalGeneralGauss::write(GXmlElement&)"
#define G_MC     "GModelSpatialEllipticalGeneralGauss::mc(GEnergy&, GTime&, "\
                                                                     "GRan&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_SMALL_ANGLE_APPROXIMATION      //!< Use small angle approximation

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpatialEllipticalGeneralGauss::GModelSpatialEllipticalGeneralGauss(void) :
                                     GModelSpatialElliptical()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Elliptical Gaussian constructor
 *
 * @param[in] dir Centre of elliptical Gaussian.
 * @param[in] semimajor Semi-major axis (deg).
 * @param[in] semiminor Semi-minor axis (deg).
 * @param[in] posangle Position angle of semi-major axis (deg).
 * @param[in] ridx Reciprocal of radial profile index.
 * @param[in] coordsys Coordinate system (either "CEL" or "GAL")
 *
 * @exception GException::invalid_argument
 *            Invalid @p coordsys argument specified.
 *
 * Construct elliptical Gaussian model from the ellipse centre direction
 * (@p dir), the @p semimajor and @p semiminor axes, the position
 * angle (@p posangle), and the reciprocal of the radial profile index
 * (@p ridx). The @p coordsys parameter specifies whether the sky direction
 * should be interpreted in the celestial or Galactic coordinate system.
 ***************************************************************************/
GModelSpatialEllipticalGeneralGauss::GModelSpatialEllipticalGeneralGauss(const GSkyDir&     dir,
                                                                         const double&      semimajor,
                                                                         const double&      semiminor,
                                                                         const double&      posangle,
                                                                         const double&      ridx,
                                                                         const std::string& coordsys) :
                                     GModelSpatialElliptical()
{
    // Throw an exception if the coordinate system is invalid
    if ((coordsys != "CEL") && (coordsys != "GAL")) {
        std::string msg = "Invalid coordinate system \""+coordsys+"\" "
                          "specified. Please specify either \"CEL\" or "
                          "\"GAL\".";
        throw GException::invalid_argument(G_CONSTRUCTOR, msg);
    }

    // Initialise members
    init_members();

    // Set parameter names
    if (coordsys == "CEL") {
        m_lon.name("RA");
        m_lat.name("DEC");
    }
    else {
        m_lon.name("GLON");
        m_lat.name("GLAT");
    }

    // Assign parameters
    this->dir(dir);
    this->semiminor(semiminor);
    this->semimajor(semimajor);
    this->posangle(posangle);
    this->ridx(ridx);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs generalised elliptical Gaussian model by extracting information
 * from an XML element. See the read() method for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpatialEllipticalGeneralGauss::GModelSpatialEllipticalGeneralGauss(const GXmlElement& xml) :
                                     GModelSpatialElliptical()
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
 * @param[in] model Generalised elliptical Gaussian model.
 ***************************************************************************/
GModelSpatialEllipticalGeneralGauss::GModelSpatialEllipticalGeneralGauss(const GModelSpatialEllipticalGeneralGauss& model) :
                                     GModelSpatialElliptical(model)
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
GModelSpatialEllipticalGeneralGauss::~GModelSpatialEllipticalGeneralGauss(void)
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
 * @param[in] model Generalised elliptical Gaussian model.
 * @return Generalised elliptical Gaussian model.
 ***************************************************************************/
GModelSpatialEllipticalGeneralGauss& GModelSpatialEllipticalGeneralGauss::operator=(const GModelSpatialEllipticalGeneralGauss& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialElliptical::operator=(model);

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
 * @brief Clear generalised elliptical Gaussian model
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatialElliptical::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialElliptical::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone generalised elliptical Gaussian model
 *
 * @return Pointer to deep copy of generalised elliptical Gaussian model.
 ***************************************************************************/
GModelSpatialEllipticalGeneralGauss* GModelSpatialEllipticalGeneralGauss::clone(void) const
{
    // Clone generalised elliptical Gaussian model
    return new GModelSpatialEllipticalGeneralGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from disk centre (radians).
 * @param[in] posangle Position angle (counterclockwise from North) (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * Evaluates the spatial component for a generalised elliptical Gaussian
 * source model.
 *
 * The generalised elliptical Gaussian function is defined by
 * \f$S_{\rm p}(\theta, \phi | E, t)\f$, where
 * \f$\theta\f$ is the angular separation between elliptical Gaussian centre
 * and the actual location and \f$\phi\f$ the position angle with respect to
 * the model centre, counted counterclockwise from North.
 *
 * The function \f$S_{\rm p}(\theta, \phi | E, t)\f$ is given by
 *
 * \f[
 * S_{\rm p}(\theta, \phi | E, t) = {\tt m\_norm} \times
 * \exp \left( -\left(\frac{\theta}{r_{\rm eff}}\right)^1/\eta \right)
 * \f]
 *
 * where the effective ellipse radius \f$r_{\rm eff}\f$ towards a given
 * position angle is given by
 *
 * \f[
 * r_{\rm eff} = \frac{ab}
 *                   {\sqrt{\left( a \sin (\phi - \phi_0) \right)^2} +
 *                    \sqrt{\left( b \cos (\phi - \phi_0) \right)^2}}
 * \f]
 * and  
 * \f$a\f$ is the semi-major axis of the ellipse,
 * \f$b\f$ is the semi-minor axis,
 * \f$\phi_0\f$ is the position angle of the ellipse, counted
 * counterclockwise from North, and
 * \f$\et\f$ is the reciprocal of the radial profile index
 * \f${\tt m\_norm}\f$ is a normalization constant given by
 *
 * \f[
 *    {\tt m\_norm} = \frac{1}{2 \pi a b \eta \Gamma(\eta)}
 * \f]
 *
 * The method will not compute analytical parameter gradients, even if the
 * @p gradients argument is set to true. Parameter gradients
 * need to be computed numerically.
 *
 * @warning
 * The above formula for an elliptical generalised Gaussian are accurate for small
 * angles, with semimajor and semiminor axes below a few degrees. This
 * should be fine for most use cases, but you should be aware that the
 * ellipse will get distorted for larger angles.
 * Note that the model normalisation is only correct in the small angle
 * approximation and for  \f$\eta\f$ order unity or smaller.
 ***************************************************************************/
double GModelSpatialEllipticalGeneralGauss::eval(const double&  theta,
                                                 const double&  posangle,
                                                 const GEnergy& energy,
                                                 const GTime&   time,
                                                 const bool&    gradients) const
{
    // Initialise value
    double value = 0.0;

    // Continue only if we're inside circle enclosing the ellipse
    if (theta <= theta_max()) {

        // Update precomputation cache
        update();

        // Perform computations to compute exponent
        #if defined(G_SMALL_ANGLE_APPROXIMATION)
        double rel_posangle  = posangle - this->posangle() * gammalib::deg2rad;
        double cosinus       = std::cos(rel_posangle);
        double sinus         = std::sin(rel_posangle);
        double arg1          = m_minor_rad * cosinus;
        double arg2          = m_major_rad * sinus;
        double r_ellipse     = m_minor_rad * m_major_rad /
                               std::sqrt(arg1*arg1 + arg2*arg2); //!< small angle
        double r_relative    = theta/r_ellipse;
        double exponent      = std::pow(r_relative,m_inv_ridx);
        #else
        // Perform computations
        double sinphi = std::sin(posangle);
        double cosphi = std::cos(posangle);

        // Compute help term
        double term1 = m_term1 * cosphi * cosphi;
        double term2 = m_term2 * sinphi * sinphi;
        double term3 = m_term3 * sinphi * cosphi;

        // Compute exponent
        double exponent = std::pow(theta * std::sqrt(term1 + term2 + term3), m_inv_ridx) ;        
        #endif

        // Set value
        value = m_norm * std::exp(-exponent);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GModelSpatialEllipticalGeneralGauss::eval";
            std::cout << "(theta=" << theta << "): NaN/Inf encountered";
            std::cout << "(posangle=" << posangle << "): NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", MinorAxis^2=" << m_minor2;
            std::cout << ", MaxjorAxis^2=" << m_major2;
            std::cout << ", m_norm=" << m_norm;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: position was inside enclosing circle

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
 * Draws an arbitrary sky direction from the elliptical generalised
 * Gaussian model.
 ***************************************************************************/
GSkyDir GModelSpatialEllipticalGeneralGauss::mc(const GEnergy& energy,
                                                const GTime&   time,
                                                GRan&          ran) const
{
    // Update precomputation cache
    update();

    // Initialise simulation
    #if defined(G_DEBUG_MC)
    int    n_samples = 0;
    #endif
    double u_max     = m_norm * std::sin(theta_max());
    double value     = 0.0;
    double u         = 1.0;
    double theta     = 0.0;
    double phi       = 0.0;

    // Throw an exception if the maximum value is zero
    if (u_max <= 0.0) {
        std::string msg = "Non positive maximum ellpitical general Gauss model value.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Simulate offset from photon arrival direction using a rejection method
    do {
        theta = ran.uniform() * theta_max();
        phi   = ran.uniform() * gammalib::twopi;
        value = eval(theta, phi, energy, time) * std::sin(theta);
        u     = ran.uniform() * u_max;
        #if defined(G_DEBUG_MC)
        n_samples++;
        #endif
    } while (u > value);
    #if defined(G_DEBUG_MC)
    std::cout << "#=" << n_samples << " ";
    #endif

    // Rotate pointing direction by offset and azimuth angle
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi * gammalib::rad2deg, theta * gammalib::rad2deg);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Checks whether model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (deg).
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the elliptical gauss
 * model.
 *
 * @todo Implement correct evaluation of effective ellipse radius.
 ***************************************************************************/
bool GModelSpatialEllipticalGeneralGauss::contains(const GSkyDir& dir,
                                                   const double&  margin) const
{
    // Compute distance to centre (radian)
    double distance = dir.dist(this->dir());

    // Return flag
    return (distance <= theta_max() + margin*gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Returns maximum model radius (radians).
 *
 * Returns the maximum of 3.0 semimajor() and 3.0 semiminor() as
 * approximate edge of the Gaussian. This limit is of course arbitrary, but
 * allows to limit the integration region for response computation. The value
 * of 3.0 has been determined by experiment (#1561).
 ***************************************************************************/
double GModelSpatialEllipticalGeneralGauss::theta_max(void) const
{
    // Set maximum model radius
    double theta_max = (semimajor() > semiminor())
                       ? semimajor() * gammalib::deg2rad * c_theta_max
                       : semiminor() * gammalib::deg2rad * c_theta_max;

    // Return value
    return theta_max;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the elliptical gauss model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="EllipticalGeneralGaussian">
 *       <parameter name="RA"          scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *       <parameter name="R_Index"     scale="1.0" value="0.5"     min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="EllipticalGeneralGaussian">
 *       <parameter name="GLON"        scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"        scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *       <parameter name="R_Index"     scale="1.0" value="0.5"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement a test of the ellipse boundary. The axes
 *       and axes minimum should be >0.
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 6);

    // Read gauss location
    GModelSpatialElliptical::read(xml);

    // Get parameters
    const GXmlElement* minor = gammalib::xml_get_par(G_READ, xml, m_semiminor.name());
    const GXmlElement* major = gammalib::xml_get_par(G_READ, xml, m_semimajor.name());
    const GXmlElement* ridx  = gammalib::xml_get_par(G_READ, xml, m_ridx.name());

    // Read parameters
    m_semiminor.read(*minor);
    m_semimajor.read(*major);
    m_ridx.read(*ridx);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Write the elliptical gauss model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="EllipticalGeneralGaussian">
 *       <parameter name="RA"          scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *       <parameter name="R_Index"     scale="1.0" value="0.5"     min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Write disk location
    GModelSpatialElliptical::write(xml);

    // Get or create parameters
    GXmlElement* minor = gammalib::xml_need_par(G_WRITE, xml, m_semiminor.name());
    GXmlElement* major = gammalib::xml_need_par(G_WRITE, xml, m_semimajor.name());
    GXmlElement* ridx  = gammalib::xml_need_par(G_WRITE, xml, m_ridx.name());

    // Write parameters
    m_semiminor.write(*minor);
    m_semimajor.write(*major);
    m_ridx.write(*ridx);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialEllipticalGeneralGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialEllipticalGeneralGauss ===");

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
void GModelSpatialEllipticalGeneralGauss::init_members(void)
{
    // Initialise model type
    m_type = "EllipticalGeneralGaussian";

    // Initialise Reciprocal index
    m_ridx.clear();
    m_ridx.name("R_Index");
    m_ridx.value(0.5);    // Gaussian
    m_ridx.min(1.e-15);   // need to be > 0
    m_ridx.free();
    m_ridx.scale(1.0);
    m_ridx.gradient(0.0);
    m_ridx.has_grad(false);

    // Set parameter pointer(s)
    m_pars.push_back(&m_ridx);

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius is not meaningful
    m_last_minor    = 0.0;
    m_last_major    = 0.0;
    m_minor_rad     = 0.0;
    m_major_rad     = 0.0;
    m_norm          = 0.0;
    m_last_posangle = 9999.0; // Signals that has not been initialised
    m_sin2pos       = 0.0;
    m_cospos2       = 0.0;
    m_sinpos2       = 0.0;
    m_minor2        = 0.0;
    m_major2        = 0.0;
    m_term1         = 0.0;
    m_term2         = 0.0;
    m_term3         = 0.0;
    m_last_ridx     = 0.0;
    m_inv_ridx      = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Elliptical Gaussian model.
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::copy_members(const GModelSpatialEllipticalGeneralGauss& model)
{
    // Copy members
    m_type = model.m_type;   // Needed to conserve model type
    m_ridx = model.m_ridx;

    // Copy precomputation cache
    m_last_minor    = model.m_last_minor;
    m_last_major    = model.m_last_major;
    m_minor_rad     = model.m_minor_rad;
    m_major_rad     = model.m_major_rad;
    m_norm          = model.m_norm;
    m_last_posangle = model.m_last_posangle;
    m_sin2pos       = model.m_sin2pos;
    m_cospos2       = model.m_cospos2;
    m_sinpos2       = model.m_sinpos2;
    m_minor2        = model.m_minor2;
    m_major2        = model.m_major2;
    m_term1         = model.m_term1;
    m_term2         = model.m_term2;
    m_term3         = model.m_term3;
    m_last_ridx     = model.m_last_ridx;
    m_inv_ridx      = model.m_inv_ridx;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the normalization
 * \f[
 *    {\tt m\_norm} = \frac{1}{2 \pi \times a \times b}
 * \f]
 * where 
 * \f$a\f$ is the semi-major axis and
 * \f$b\f$ is the semi-minor axis of the ellipse.
 *
 * @warning
 * The normalization of the elliptical Gaussian is only valid in the small
 * angle approximation.
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::update() const
{
    // Initialise flag if something has changed
    bool changed = false;

    // Update if one axis has changed or radial profile index has changed
    if (m_last_minor != semiminor() || m_last_major != semimajor() || m_last_ridx != ridx()) {

        // Signal parameter changes
        changed = true;

        // Store last values
        m_last_minor = semiminor();
        m_last_major = semimajor();
        m_last_ridx = ridx();

        // Compute axes in radians
        m_minor_rad = m_last_minor * gammalib::deg2rad;
        m_major_rad = m_last_major * gammalib::deg2rad;

        // Perform precomputations. Note also
        // that the model normalisation is only correct in the small angle
        // approximation and for  ridx order unity or smaller.
        if (m_minor_rad > 0.0 and m_major_rad > 0.0 and ridx() > 0.0) {
            m_inv_ridx = 1.0 / ridx() ;
            m_norm     = 1.0 / (gammalib::twopi * m_minor_rad * m_major_rad * ridx() *
                                std::exp(gammalib::gammln(2.0 * ridx())) );
        }
        else {
            m_inv_ridx = 1.0;
            m_norm     = 0.0;
        }

    } // endif: update required

    // Update chache if position angle changed
    #if !defined(G_SMALL_ANGLE_APPROXIMATION)
    if (m_last_posangle != posangle()) {

        // Signal parameter changes
        changed = true;

        // Store last value
        m_last_posangle = posangle();

        // Compute angle in radians
        double posangle_rad = m_last_posangle * gammalib::deg2rad;

        // Compute sine and cosine
        double cospos = std::cos(posangle_rad);
        double sinpos = std::sin(posangle_rad);

        // Cache important values for further computations
        m_cospos2 = cospos * cospos;
        m_sinpos2 = sinpos * sinpos;
        m_sin2pos = std::sin(2.0 * posangle_rad);

    } // endif: position angle update required

    // Perform precomputations in case anything has changed
    if (changed) {

        // Compute help terms
        m_term1 = 0.5 * (m_cospos2 / m_minor2 + m_sinpos2 / m_major2);
        m_term2 = 0.5 * (m_cospos2 / m_major2 + m_sinpos2 / m_minor2);
        m_term3 = 0.5 * (m_sin2pos / m_major2 - m_sin2pos / m_minor2);

    } // endif: something has changed
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialEllipticalGeneralGauss::set_region(void) const
{
    // Set maximum model radius
    double max_radius = (semimajor() > semiminor()) ? semimajor() : semiminor();

    // Set sky region circle (maximum Gaussian sigma times a scaling
    // factor (actually 5))
    GSkyRegionCircle region(dir(), max_radius * c_theta_max);

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialEllipticalGeneralGauss*>(this)->m_region = region;

    // Return
    return;
}
