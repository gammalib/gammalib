/***************************************************************************
 *  GModelSpatialEllipticalGauss.cpp - Elliptical gauss source model class *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Michael Mayer                               *
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
 * @file GModelSpatialEllipticalGauss.cpp
 * @brief Elliptical gauss model class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialEllipticalGauss.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */
namespace {
    const double c_theta_max    = 3.0; //!< semiaxis multiplied for theta_max()
    const double c_max_exponent = 0.5 * c_theta_max * c_theta_max;
    const double c_fraction     = 1.0 - std::exp(-c_max_exponent);
}

/* __ Globals ____________________________________________________________ */
const GModelSpatialEllipticalGauss g_elliptical_gauss_seed;
const GModelSpatialRegistry        g_elliptical_gauss_registry(&g_elliptical_gauss_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialEllipticalGauss g_elliptical_gauss_legacy_seed(true, "EllipticalGauss");
const GModelSpatialRegistry        g_elliptical_gauss_legacy_registry(&g_elliptical_gauss_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_READ             "GModelSpatialEllipticalGauss::read(GXmlElement&)"
#define G_WRITE           "GModelSpatialEllipticalGauss::write(GXmlElement&)"

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
GModelSpatialEllipticalGauss::GModelSpatialEllipticalGauss(void) :
                              GModelSpatialElliptical()
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
 * Constructs empty elliptical Gaussian model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialEllipticalGauss::GModelSpatialEllipticalGauss(const bool&        dummy,
                                                           const std::string& type) :
                              GModelSpatialElliptical()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Elliptical Gaussian constructor
 *
 * @param[in] dir Centre of elliptical Gaussian.
 * @param[in] semimajor Semi-major axis (degrees).
 * @param[in] semiminor Semi-minor axis (degrees).
 * @param[in] posangle Position angle of semi-major axis (degrees).
 *
 * Construct elliptical Gaussian model from the ellipse centre direction
 * (@p dir), the @p semimajor and @p semiminor axes, and the position
 * angle (@p posangle).
 ***************************************************************************/
GModelSpatialEllipticalGauss::GModelSpatialEllipticalGauss(const GSkyDir& dir,
                                                           const double&  semimajor,
                                                           const double&  semiminor,
                                                           const double&  posangle) :
                              GModelSpatialElliptical()
{
    // Initialise members
    init_members();

    // Assign parameters
    this->dir(dir);
    this->semiminor(semiminor);
    this->semimajor(semimajor);
    this->posangle(posangle);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs elliptical Gaussian model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialEllipticalGauss::GModelSpatialEllipticalGauss(const GXmlElement& xml) :
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
 * @param[in] model Elliptical Gaussian model.
 ***************************************************************************/
GModelSpatialEllipticalGauss::GModelSpatialEllipticalGauss(const GModelSpatialEllipticalGauss& model) :
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
GModelSpatialEllipticalGauss::~GModelSpatialEllipticalGauss(void)
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
 * @param[in] model Elliptical gauss model.
 * @return Elliptical Gaussian model.
 ***************************************************************************/
GModelSpatialEllipticalGauss& GModelSpatialEllipticalGauss::operator=(const GModelSpatialEllipticalGauss& model)
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
 * @brief Clear elliptical Gaussian model
 ***************************************************************************/
void GModelSpatialEllipticalGauss::clear(void)
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
 * @brief Clone elliptical Gaussian model
 *
 * @return Pointer to deep copy of elliptical Gaussian model.
 ***************************************************************************/
GModelSpatialEllipticalGauss* GModelSpatialEllipticalGauss::clone(void) const
{
    // Clone elliptical Gaussian model
    return new GModelSpatialEllipticalGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from disk centre (radians).
 * @param[in] posangle Position angle (counterclockwise from North) (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the spatial component for an elliptical Gaussian source model.
 *
 * The elliptical Gaussian function is defined by
 * \f$S_{\rm p}(\theta, \phi | E, t)\f$, where
 * \f$\theta\f$ is the angular separation between elliptical Gaussian centre
 * and the actual location and \f$\phi\f$ the position angle with respect to
 * the model centre, counted counterclockwise from North.
 *
 * The function \f$S_{\rm p}(\theta, \phi | E, t)\f$ is given by
 *
 * \f[
 * S_{\rm p}(\theta, \phi | E, t) = {\tt m\_norm} \times
 * \exp \left( -\frac{\theta^2}{2 r_{\rm eff}^2} \right)
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
 * \f$b\f$ is the semi-minor axis, and
 * \f$\phi_0\f$ is the position angle of the ellipse, counted
 * counterclockwise from North.
 * \f${\tt m\_norm}\f$ is a normalization constant given by
 *
 * \f[
 *    {\tt m\_norm} = \frac{1}{2 \pi \times a \times b}
 * \f]
 *
 * @warning
 * The above formula for an elliptical Gaussian are accurate for small
 * angles, with semimajor and semiminor axes below a few degrees. This
 * should be fine for most use cases, but you should be aware that the
 * ellipse will get distorted for larger angles.
 *
 * @warning
 * For numerical reasons the elliptical Gaussian will be truncated for
 * \f$\theta\f$ angles that correspond to 3.0 times the effective ellipse
 * radius.
 ***************************************************************************/
double GModelSpatialEllipticalGauss::eval(const double&  theta,
                                          const double&  posangle,
                                          const GEnergy& energy,
                                          const GTime&   time) const
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
        double exponent      = 0.5*r_relative*r_relative;
        #else
        // Perform computations
        double sinphi = std::sin(posangle);
        double cosphi = std::cos(posangle);

        // Compute help term
        double term1 = m_term1 * cosphi * cosphi;
        double term2 = m_term2 * sinphi * sinphi;
        double term3 = m_term3 * sinphi * cosphi;

        // Compute exponent
        double exponent = theta * theta * (term1 + term2 + term3);        
        #endif

        // Set value if the exponent is within the truncation region
        if (exponent <= c_max_exponent) {
            value = m_norm * std::exp(-exponent);
        }

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GModelSpatialEllipticalGauss::eval";
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
 * @brief Evaluate function and gradients (in units of sr^-1)
 *
 * @param[in] theta Angular distance from gaussian centre (radians).
 * @param[in] posangle Position angle (counterclockwise from North) (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the function value. No gradient computation is implemented as
 * elliptical models will be convolved with the instrument response and thus
 * require the numerical computation of the derivatives.
 *
 * See the eval() method for more information.
 ***************************************************************************/
double GModelSpatialEllipticalGauss::eval_gradients(const double&  theta,
                                                    const double&  posangle,
                                                    const GEnergy& energy,
                                                    const GTime&   time) const
{
    // Return value
    return (eval(theta, posangle, energy, time));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction from the elliptical Gaussian model.
 *
 * @warning
 * For numerical reasons the elliptical Gaussian will be truncated for
 * \f$\theta\f$ angles that correspond to 3.0 times the effective ellipse
 * radius.
 ***************************************************************************/
GSkyDir GModelSpatialEllipticalGauss::mc(const GEnergy& energy,
                                         const GTime&   time,
                                         GRan&          ran) const
{
    // Update precomputation cache
    update();

    // Initialise photon
    GPhoton photon;
    photon.energy(energy);
    photon.time(time);

    // Draw gaussian offset from each axis
    double ran_major;
    double ran_minor;
    do {
        ran_major = ran.normal();
    } while (ran_major > c_theta_max);
    do {
        ran_minor = ran.normal();
    } while (ran_minor > c_theta_max);
    double theta1 = semimajor() * ran_major;
    double theta2 = semiminor() * ran_minor;

    // Compute total offset from model centre in small angle approximation
    double theta = std::sqrt(theta1 * theta1 + theta2 * theta2);

    // Compute rotation angle, taking into account given position angle
    double phi = gammalib::atan2d(theta2, theta1) + posangle();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi , theta);

    // Set photon sky direction
    photon.dir(sky_dir);

    // Return photon direction
    return (photon.dir());
}


/***********************************************************************//**
 * @brief Checks whether model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the elliptical gauss
 * model.
 *
 * @todo Implement correct evaluation of effective ellipse radius.
 ***************************************************************************/
bool GModelSpatialEllipticalGauss::contains(const GSkyDir& dir,
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
 * @return Returns maximum model radius.
 *
 * Returns the maximum of 3.0 semimajor() and 3.0 semiminor() as
 * approximate edge of the Gaussian. This limit is of course arbitrary, but
 * allows to limit the integration region for response computation. The value
 * of 3.0 has been determined by experiment (#1561).
 ***************************************************************************/
double GModelSpatialEllipticalGauss::theta_max(void) const
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
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Reads the elliptical gauss model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="EllipticalGaussian">
 *       <parameter name="RA"          scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="EllipticalGaussian">
 *       <parameter name="GLON"        scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"        scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement a test of the ellipse boundary. The axes
 *       and axes minimum should be >0.
 ***************************************************************************/
void GModelSpatialEllipticalGauss::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 5 parameters
    if (xml.elements() != 5 || npars != 5) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Elliptical gauss model requires exactly 5 parameters.");
    }

    // Read gauss location
    GModelSpatialElliptical::read(xml);

    // Extract model parameters
    int  npar[2] = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle semiminor radius
        if (par->attribute("name") == "MinorRadius") {
            
            // Read parameter
            m_semiminor.read(*par);
            
            // Increment parameter counter
            npar[0]++;
        }

        // Handle semimajor radius
        else if (par->attribute("name") == "MajorRadius") {

        	// Read parameter
        	m_semimajor.read(*par);

        	// Increment parameter counter
        	npar[1]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Elliptical disk model requires \"MinorRadius\" and"
              " \"MajorRadius\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the elliptical gauss model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="EllipticalGaussian">
 *       <parameter name="RA"          scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialEllipticalGauss::write(GXmlElement& xml) const
{
    // Write disk location
    GModelSpatialElliptical::write(xml);

    // If XML element has 3 nodes (which should be the location and PA nodes)
    // then append 2 parameter nodes
    if (xml.elements() == 3) {
        xml.append(GXmlElement("parameter name=\"MinorRadius\""));
        xml.append(GXmlElement("parameter name=\"MajorRadius\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 5 parameters
    if (xml.elements() != 5 || npars != 5) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Elliptical gauss model requires exactly 5 parameters.");
    }

    // Set or update model parameter attributes
    int npar[2] = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle semiminor radius
        if (par->attribute("name") == "MinorRadius") {

        	// Write parameter
            m_semiminor.write(*par);

            // Increment parameter counter
            npar[0]++;
        }

        // Handle semimajor radius
        else if (par->attribute("name") == "MajorRadius") {

        	// Write parameter
            m_semimajor.write(*par);

            // Increment parameter counter
            npar[1]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Elliptical gauss model requires \"MinorRadius\" and"
              " \"MajorRadius\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialEllipticalGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialEllipticalGauss ===");

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
void GModelSpatialEllipticalGauss::init_members(void)
{
    // Initialise model type
    m_type = "EllipticalGaussian";

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius is not meaningful
    m_last_minor        = 0.0;
    m_last_major        = 0.0;
    m_minor_rad         = 0.0;
    m_major_rad         = 0.0;
    m_norm              = 0.0;
    m_last_posangle     = 9999.0; // Signals that has not been initialised
    m_sin2pos           = 0.0;
    m_cospos2           = 0.0;
    m_sinpos2           = 0.0;
    m_minor2            = 0.0;
    m_major2            = 0.0;
    m_term1             = 0.0;
    m_term2             = 0.0;
    m_term3             = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Elliptical Gaussian model.
 ***************************************************************************/
void GModelSpatialEllipticalGauss::copy_members(const GModelSpatialEllipticalGauss& model)
{
    // Copy members
    m_type = model.m_type;

    // Copy precomputation cache
    m_last_minor        = model.m_last_minor;
    m_last_major        = model.m_last_major;
    m_minor_rad         = model.m_minor_rad;
    m_major_rad         = model.m_major_rad;
    m_norm              = model.m_norm;
    m_last_posangle     = model.m_last_posangle;
    m_sin2pos           = model.m_sin2pos;
    m_cospos2           = model.m_cospos2;
    m_sinpos2           = model.m_sinpos2;
    m_minor2            = model.m_minor2;
    m_major2            = model.m_major2;
    m_term1             = model.m_term1;
    m_term2             = model.m_term2;
    m_term3             = model.m_term3;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialEllipticalGauss::free_members(void)
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
void GModelSpatialEllipticalGauss::update() const
{
    // Initialise flag if something has changed
    bool changed = false;

    // Update if one axis has changed
    if (m_last_minor != semiminor() || m_last_major != semimajor()) {

        // Signal parameter changes
        changed = true;

        // Store last values
        m_last_minor = semiminor();
        m_last_major = semimajor();

        // Compute axes in radians
        m_minor_rad = m_last_minor * gammalib::deg2rad;
        m_major_rad = m_last_major * gammalib::deg2rad;

        // Perform precomputations. Note that I verified the normalization
        // by numerical integration of the resulting Gaussian. Note also
        // also that the normalization is only correct in the small angle
        // approximation.
        m_minor2     = m_minor_rad * m_minor_rad;
        m_major2     = m_major_rad * m_major_rad;
        double denom = gammalib::twopi * m_minor_rad * m_major_rad *
                       c_fraction;
        m_norm       = (denom > 0.0) ? 1.0 / denom : 0.0;

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
