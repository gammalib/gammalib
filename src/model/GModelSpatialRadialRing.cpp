/***************************************************************************
 *      GModelSpatialRadialRing.cpp - Radial ring source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Pierrick Martin                             *
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
 * @file GModelSpatialRadialRing.cpp
 * @brief Radial ring model class implementation
 * @author Pierrick Martin
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialRadialRing.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialRing g_radial_ring_seed;
const GModelSpatialRegistry   g_radial_ring_registry(&g_radial_ring_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialRadialRing g_radial_ring_legacy_seed(true, "RingFunction");
const GModelSpatialRegistry   g_radial_ring_legacy_registry(&g_radial_ring_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_READ                  "GModelSpatialRadialRing::read(GXmlElement&)"
#define G_WRITE                "GModelSpatialRadialRing::write(GXmlElement&)"

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
 * Constructs empty radial ring model.
 ***************************************************************************/
GModelSpatialRadialRing::GModelSpatialRadialRing(void) : GModelSpatialRadial()
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
 * Constructs empty radial ring model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialRadialRing::GModelSpatialRadialRing(const bool&        dummy,
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
 * @brief Ring constructor
 *
 * @param[in] dir Sky position of ring centre.
 * @param[in] radius Ring outer radius (degrees).
 * @param[in] width Ring width (degrees).
 *
 * Constructs radial ring model from the sky position of the ring centre
 * (@p dir) and the ring outer radius @p radius and width 
 * @p width in degrees.
 ***************************************************************************/
GModelSpatialRadialRing::GModelSpatialRadialRing(const GSkyDir& dir,
                                                 const double&  radius,
                                                 const double&  width) :
                         GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Assign center location
    this->dir(dir);
    
    // Assign radius and width after check
    if (width <= radius)  {
        this->radius(radius);
        this->width(width);
    } else {
        this->radius(radius);
        this->width(radius);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of radial ring model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialRing::GModelSpatialRadialRing(const GXmlElement& xml) :
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
 * @param[in] model Radial ring model.
 ***************************************************************************/
GModelSpatialRadialRing::GModelSpatialRadialRing(const GModelSpatialRadialRing& model) :
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
GModelSpatialRadialRing::~GModelSpatialRadialRing(void)
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
 * @param[in] model Radial ring model.
 * @return Radial ring model.
 ***************************************************************************/
GModelSpatialRadialRing& GModelSpatialRadialRing::operator=(const GModelSpatialRadialRing& model)
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
 * @brief Clear radial ring model
 ***************************************************************************/
void GModelSpatialRadialRing::clear(void)
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
 * @brief Clone radial ring model
 *
 * @return Pointer to deep copy of radial ring model.
 ***************************************************************************/
GModelSpatialRadialRing* GModelSpatialRadialRing::clone(void) const
{
    // Clone radial ring model
    return new GModelSpatialRadialRing(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from ring centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * Evaluates the spatial component for a ring source model using
 *
 * \f[
 *    S_{\rm p}(\vec{p} | E, t) = \left \{
 *    \begin{array}{l l}
 *       \displaystyle
 *       {\tt m\_norm}
 *       & \mbox{if $\theta \le $ outer radius and $\theta \ge $ inner radius} \\
 *       \\
 *      \displaystyle
 *      0 & \mbox{if $\theta > $ outer radius or $\theta < $ inner radius}
 *    \end{array}
 *    \right .
 * \f]
 *
 * where
 * - \f$\theta\f$ is the angular separation between ring centre and the
 *   sky direction of interest, and
 * - \f${\tt m\_norm} = \frac{1}{2 \pi (\cos r_i - \cos r_o)} \f$ is a normalization
 *   constant (see the update() method for details).
 *
 * The method will not compute analytical parameter gradients, even if the
 * @p gradients argument is set to true. Radial ring parameter gradients
 * need to be computed numerically.
 ***************************************************************************/
double GModelSpatialRadialRing::eval(const double&  theta,
                                     const GEnergy& energy,
                                     const GTime&   time,
                                     const bool&    gradients) const
{
    // Update precomputation cache
    update();

    // Set value
    double value = 0.0;
    value = (theta <= m_outer_radius_rad) ? 1.0 : 0.0;
    value *= (theta >= m_inner_radius_rad) ? m_norm : 0.0;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpatialRadialRing::eval";
        std::cout << "(theta=" << theta << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", m_inner_radius_rad=" << m_inner_radius_rad;
        std::cout << ", m_outer_radius_rad=" << m_outer_radius_rad;
        std::cout << ", m_norm=" << m_norm;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction from the 2D ring distribution as function
 * of the photon @p energy and arrival @p time.
 ***************************************************************************/
GSkyDir GModelSpatialRadialRing::mc(const GEnergy& energy,
                                    const GTime&   time,
                                    GRan&          ran) const
{
    // Update precomputation cache
    update();
    
    // Simulate offset from photon arrival direction
    double cosirad = std::cos(m_inner_radius_rad);
    double cosorad = std::cos(m_outer_radius_rad);
    double theta  = std::acos(cosirad - ran.uniform() * (cosirad - cosorad)) *
                    gammalib::rad2deg;
    double phi    = 360.0 * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Checks whether model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the radial ring model.
 ***************************************************************************/
bool GModelSpatialRadialRing::contains(const GSkyDir& dir,
                                       const double&  margin) const
{
    // Update precomputation cache
    update();
    
    // Compute distance to centre (radians)
    double distance = dir.dist(this->dir());
    
    // Set flag
    bool flag;
    flag = (distance <= m_outer_radius_rad + margin*gammalib::deg2rad)
           *(distance >= m_inner_radius_rad - margin*gammalib::deg2rad);
    
    // Return flag
    return flag;
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialRing::theta_max(void) const
{
    // Return value
    return m_outer_radius_rad;
}


/***********************************************************************//**
 * @brief Return minimum model radius (in radians)
 *
 * @return Minimum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialRing::theta_min(void) const
{
    // Return value
    return m_inner_radius_rad;
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
 * Reads the radial ring model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="RadialRing">
 *       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Width"  scale="1.0" value="0.15"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="RadialRing">
 *       <parameter name="GLON"   scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Width"  scale="1.0" value="0.15"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 * 
 * @todo Implement a test of the radius and radius boundary. The radius
 *       and radius minimum should be >0.
 ***************************************************************************/
void GModelSpatialRadialRing::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || npars != 4) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Ring model requires exactly 4 parameters.");
    }

    // Read ring center location
    GModelSpatialRadial::read(xml);

    // Get parameters
    const GXmlElement* radius = gammalib::xml_get_par(G_READ, xml, m_radius.name());
    const GXmlElement* width = gammalib::xml_get_par(G_READ, xml, m_width.name());

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
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the radial ring model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="RadialRing">
 *       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Width"  scale="1.0" value="0.15"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadialRing::write(GXmlElement& xml) const
{
    // Write ring center location
    GModelSpatialRadial::write(xml);

    // Get or create parameters
    GXmlElement* radius = gammalib::xml_need_par(G_WRITE, xml, m_radius.name());
    GXmlElement* width = gammalib::xml_need_par(G_WRITE, xml, m_width.name());

    // Write parameters
    m_radius.write(*radius);
    m_width.write(*width);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialRing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialRing ===");

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
void GModelSpatialRadialRing::init_members(void)
{
    // Initialise model type
    m_type = "RadialRing";

    // Initialise Radius
    m_radius.clear();
    m_radius.name("Radius");
    m_radius.unit("deg");
    m_radius.value(2.778e-4); // 1 arcsec
    m_radius.min(2.778e-4);   // 1 arcsec
    m_radius.free();
    m_radius.scale(1.0);
    m_radius.gradient(0.0);
    m_radius.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_radius);
    
    // Initialise Width
    m_width.clear();
    m_width.name("Width");
    m_width.unit("deg");
    m_width.value(2.778e-4); // 1 arcsec
    m_width.min(2.778e-4);   // 1 arcsec
    m_width.free();
    m_width.scale(1.0);
    m_width.gradient(0.0);
    m_width.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_width);
    
    // Initialise precomputation cache
    m_last_radius       = 0.0;
    m_last_width        = 0.0;
    m_inner_radius_rad  = 0.0;
    m_outer_radius_rad  = 0.0;
    m_width_rad         = 0.0;
    m_norm              = 0.0;

    // Initialise other members
    m_region.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial ring model.
 *
 * We do not have to push back the members on the parameter stack as this
 * should have been done by init_members() that was called before. Otherwise
 * we would have the radius twice on the stack.
 ***************************************************************************/
void GModelSpatialRadialRing::copy_members(const GModelSpatialRadialRing& model)
{
    // Copy members
    m_type   = model.m_type;
    m_radius = model.m_radius;
    m_width  = model.m_width;
    m_region = model.m_region;

    // Copy precomputation cache
    m_last_radius       = model.m_last_radius;
    m_last_width        = model.m_last_width;
    m_inner_radius_rad  = model.m_inner_radius_rad;
    m_outer_radius_rad  = model.m_outer_radius_rad;
    m_width_rad         = model.m_width_rad;
    m_norm              = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialRing::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the normalization
 * \f[
 *    {\tt m\_norm} = \frac{1}{2 \pi (\cos r_i - \cos r_o)}
 * \f]
 *
 * This is the correct normalization on the sphere for any radii.
 ***************************************************************************/
void GModelSpatialRadialRing::update() const
{
    // Update if radius or width has changed
    if ((m_last_radius != radius()) || (m_last_width != width())) {

        // Store last values
        m_last_radius = radius();
        m_last_width  = width();

        // Compute ring parameters in radians
        m_outer_radius_rad = radius() * gammalib::deg2rad;
        m_inner_radius_rad = (radius()-width()) * gammalib::deg2rad;
        m_width_rad = width() * gammalib::deg2rad;

        // Perform precomputations
        double denom = gammalib::twopi * 
                       (std::cos(m_inner_radius_rad) - 
                        std::cos(m_outer_radius_rad));
        m_norm       = (denom > 0.0) ? 1.0 / denom : 0.0;

    } // endif: update required
    
    // Also update region
    set_region();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialRadialRing::set_region(void) const
{
    // Set sky region centre to ring centre
    m_region.centre(m_ra.value(), m_dec.value());

    // Set sky region radius to ring radius
    m_region.radius(m_radius.value());

    // Return
    return;
}
