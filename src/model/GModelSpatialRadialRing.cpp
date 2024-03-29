/***************************************************************************
 *      GModelSpatialRadialRing.cpp - Radial ring source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020-2022 by Pierrick Martin                             *
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

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR     "GModelSpatialRadialRing::GModelSpatialRadialRing("\
                                  "GSkyDir&, double&, double&, std::string&)"
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
 * @brief Ring constructor
 *
 * @param[in] dir Sky position of ring centre.
 * @param[in] radius Ring inner radius (degrees).
 * @param[in] width Ring width (degrees).
 * @param[in] coordsys Coordinate system (either "CEL" or "GAL")
 *
 * @exception GException::invalid_argument
 *            Invalid @p coordsys argument specified.
 *
 * Constructs radial ring model from the sky position of the ring centre
 * (@p dir) and the ring inner radius (@p radius) and width (@p width) in
 * degrees. The @p coordsys parameter specifies whether the sky direction
 * should be interpreted in the celestial or Galactic coordinate system.
 ***************************************************************************/
GModelSpatialRadialRing::GModelSpatialRadialRing(const GSkyDir&     dir,
                                                 const double&      radius,
                                                 const double&      width,
                                                 const std::string& coordsys) :
                         GModelSpatialRadial()
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

    // Assign center location, radius and width
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
    if ((theta >= m_inner_radius_rad) && (theta <= m_outer_radius_rad)) {
        value = m_norm;
    }

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
    double theta   = std::acos(m_cos_inner_radius_rad - ran.uniform() *
                               (m_cos_inner_radius_rad - m_cos_outer_radius_rad)) *
                     gammalib::rad2deg;
    double phi     = 360.0 * ran.uniform();

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
    // Compute distance to centre (radians)
    double distance = dir.dist(this->dir());

    // Compute margin in radians
    double margin_rad   = margin * gammalib::deg2rad;
    double distance_min = theta_min() - margin_rad;
    double distance_max = theta_max() + margin_rad;
    if (distance_min < 0.0) {
        distance_min = 0.0;
    }

    // Return flag
    return ((distance >= distance_min) && (distance <= distance_max));
}


/***********************************************************************//**
 * @brief Return minimum model radius (in radians)
 *
 * @return Minimum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialRing::theta_min(void) const
{
    // Return value
    return (radius() * gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialRing::theta_max(void) const
{
    // Return value
    return ((radius()+width()) * gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
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
 * @todo Implement a test of the radius and width. Both parameters should
 * be >0.
 ***************************************************************************/
void GModelSpatialRadialRing::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 4);

    // Read ring center location
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
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Write ring center location
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
    m_radius.has_grad(false); // Radial components never have gradients

    // Initialise Width
    m_width.clear();
    m_width.name("Width");
    m_width.unit("deg");
    m_width.value(2.778e-4);  // 1 arcsec
    m_width.min(2.778e-4);    // 1 arcsec
    m_width.free();
    m_width.scale(1.0);
    m_width.gradient(0.0);
    m_width.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_radius);
    m_pars.push_back(&m_width);

    // Initialise precomputation cache
    m_last_radius          = 0.0;
    m_last_width           = 0.0;
    m_inner_radius_rad     = 0.0;
    m_outer_radius_rad     = 0.0;
    m_cos_inner_radius_rad = 1.0;
    m_cos_outer_radius_rad = 1.0;
    m_norm                 = 0.0;

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
    m_type   = model.m_type;   // Needed to conserve model type
    m_radius = model.m_radius;
    m_width  = model.m_width;

    // Copy precomputation cache
    m_last_radius          = model.m_last_radius;
    m_last_width           = model.m_last_width;
    m_inner_radius_rad     = model.m_inner_radius_rad;
    m_outer_radius_rad     = model.m_outer_radius_rad;
    m_cos_inner_radius_rad = model.m_cos_inner_radius_rad;
    m_cos_outer_radius_rad = model.m_cos_outer_radius_rad;
    m_norm                 = model.m_norm;

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
        m_inner_radius_rad     = theta_min();
        m_outer_radius_rad     = theta_max();
        m_cos_inner_radius_rad = std::cos(m_inner_radius_rad);
        m_cos_outer_radius_rad = std::cos(m_outer_radius_rad);

        // Perform precomputations
        double denom = gammalib::twopi * (m_cos_inner_radius_rad -
                                          m_cos_outer_radius_rad);
        m_norm       = (denom > 0.0) ? 1.0 / denom : 0.0;

    } // endif: update required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialRadialRing::set_region(void) const
{
    // Set sky region circle (maximum Gaussian sigma times a scaling
    // factor (actually 3))
    GSkyRegionCircle region(dir(), radius()+width());

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialRadialRing*>(this)->m_region = region;

    // Return
    return;
}
