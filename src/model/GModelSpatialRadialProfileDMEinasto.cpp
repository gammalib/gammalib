/***************************************************************************
 * GModelSpatialRadialProfileDMEinasto.cpp - DMEinasto radial profile class *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialRadialProfileDMEinasto.cpp
 * @brief Radial DM Einasto profile model class implementation
 * @author Nathan Kelley-Hoskins
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GXmlElement.hpp"
#include "GModelSpatialRadialProfileDMEinasto.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */
const double mass_factor = 10.0 ; // determines volume integration radius

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialProfileDMEinasto g_radial_disk_seed;
const GModelSpatialRegistry           g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ   "GModelSpatialRadialProfileDMEinasto::read(GXmlElement&)"
#define G_WRITE  "GModelSpatialRadialProfileDMEinasto::write(GXmlElement&)"

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
 * Constructs empty radial DMEinasto profile
 ***************************************************************************/
GModelSpatialRadialProfileDMEinasto::GModelSpatialRadialProfileDMEinasto(void) :
                                 GModelSpatialRadialProfile()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs radial DMEinasto profile model by extracting information from
 * an XML element. See the read() method for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialProfileDMEinasto::GModelSpatialRadialProfileDMEinasto(const GXmlElement& xml) :
                                 GModelSpatialRadialProfile()
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
 * @param[in] model Radial DMEinasto profile model.
 *
 * Copies radial DMEinasto profile model from another radial profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMEinasto::GModelSpatialRadialProfileDMEinasto(const GModelSpatialRadialProfileDMEinasto& model) :
                                 GModelSpatialRadialProfile(model)
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
 * Destructs radial DMEinasto profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMEinasto::~GModelSpatialRadialProfileDMEinasto(void)
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
 * @param[in] model Radial DMEinasto profile model.
 * @return Radial DMEinasto profile model.
 *
 * Assigns radial DMEinasto profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMEinasto& GModelSpatialRadialProfileDMEinasto::operator=(const GModelSpatialRadialProfileDMEinasto& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialRadialProfile::operator=(model);

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
 * @brief Clear radial DMEinasto profile model
 *
 * Clears radial DMEinasto profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::clear(void)
{
    // Free class members
    free_members();
    this->GModelSpatialRadialProfile::free_members();
    this->GModelSpatialRadial::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialRadial::init_members();
    this->GModelSpatialRadialProfile::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone radial DMEinasto profile model
 *
 * @return Pointer to deep copy of radial DMEinasto profile model.
 *
 * Returns a deep copy of the radial DMEinasto profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMEinasto* GModelSpatialRadialProfileDMEinasto::clone(void) const
{
    // Clone radial disk model
    return new GModelSpatialRadialProfileDMEinasto(*this);
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::theta_max(void) const
{
    
    // update precomputation cache
    update();
    
    double theta = 0.0 ;
    
    // if earth is within the significant radius, then theta_max must
    // contain the entire profile (180deg)
    if ( m_halo_distance.value() < m_mass_radius ) {
      theta = gammalib::pi ;
    
    // if the halo is far enough away (further than the mass radius)
    // then we just need to deal with the angles within the sphere of the
    // significant radius.
    } else {
      theta = gammalib::acosd( m_mass_radius / m_halo_distance.value() ) ;
    }
    
    // Return value
    return theta ;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the DMEinasto radial profile model information from an XML element.
 * The XML element shall have either the format 
 *
 *     <spatialModel type="DMEinastoProfile">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Scale Radius"  scale="1.0" value="21.5" min="1.0e-6" free="1"/>
 *       <parameter name="Halo Distance" scale="1.0" value="7.94" min="1.0e-6" free="1"/>
 *       <parameter name="Alpha"         scale="1.0" value="0.17" min="0.01"   max="10" free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="DMEinastoProfile">
 *       <parameter name="GLON"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"  scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Scale Radius"  scale="1.0" value="21.5" min="1.0e-6" free="1"/>
 *       <parameter name="Halo Distance" scale="1.0" value="7.94" min="1.0e-6" free="1"/>
 *       <parameter name="Alpha"         scale="1.0" value="0.17" min="0.01"   max="10" free="1"/>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::read(const GXmlElement& xml)
{
    // Read DMEinasto location
    GModelSpatialRadial::read(xml);

    const GXmlElement* par1 = gammalib::xml_get_par(G_READ, xml, "Scale Radius");
    m_scale_radius.read(*par1);
    
    const GXmlElement* par2 = gammalib::xml_get_par(G_READ, xml, "Halo Distance");
    m_halo_distance.read(*par2);

    const GXmlElement* par3 = gammalib::xml_get_par(G_READ, xml, "Alpha");
    m_alpha.read(*par3);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Writes the DMEinasto radial profile model information into an XML element.
 * The XML element will have the format 
 *
 *     <spatialModel type="DMEinastoProfile">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::write(GXmlElement& xml) const
{
    // Write DMEinasto location
    GModelSpatialRadial::write(xml);

    // Write Scale Radius parameter
    GXmlElement* par1 = gammalib::xml_need_par(G_WRITE, xml, "Scale Radius");
    m_scale_radius.write(*par1);

    // Write Halo Distance parameter
    GXmlElement* par2 = gammalib::xml_need_par(G_WRITE, xml, "Halo Distance");
    m_halo_distance.write(*par2);

    // Write Alpha parameter
    GXmlElement* par3 = gammalib::xml_need_par(G_WRITE, xml, "Alpha");
    m_alpha.write(*par3);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialProfileDMEinasto::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialProfileDMEinasto ===");

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
void GModelSpatialRadialProfileDMEinasto::init_members(void)
{

    // Initialise scale radius
    m_scale_radius.clear();
    m_scale_radius.name("Scale Radius");
    m_scale_radius.unit("kpc");
    m_scale_radius.value(21.5); // default to GC scale radius
    m_scale_radius.min(1.0e-6);   // arbitrarily chosen :/
    m_scale_radius.free();
    m_scale_radius.scale(1.0);
    m_scale_radius.gradient(0.0);
    m_scale_radius.has_grad(false);  // Radial components never have gradients

    // Initialise halo distance
    m_halo_distance.clear();
    m_halo_distance.name("Halo Distance");
    m_halo_distance.unit("kpc");
    m_halo_distance.value(7.94); // default to GC halo distance
    m_halo_distance.min(1.0e-6); // arbitrarily chosen
    m_halo_distance.free();
    m_halo_distance.scale(1.0);
    m_halo_distance.gradient(0.0);
    m_halo_distance.has_grad(false);  // Radial components never have gradients

    // Initialise alpha
    m_alpha.clear();
    m_alpha.name("Alpha");
    m_alpha.unit("unitless");
    m_alpha.value(0.17); // default to GC alpha
    m_alpha.min(0.01);   // arbitrarily chosen
    m_alpha.max(10.0);   // arbitrarily chosen
    m_alpha.free();
    m_alpha.scale(1.0);
    m_alpha.gradient(0.0);
    m_alpha.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_scale_radius );
    m_pars.push_back(&m_halo_distance);
    m_pars.push_back(&m_alpha        );
    
    // Initialize precomputation cache. Note that zero values flag
    // uninitialised, as a zero radius is not meaningful
    m_last_scale_radius = 0.0 ;
    m_mass_radius = 0.0 ;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial DMEinasto model.
 *
 * Copies class members from another radial profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::copy_members(const GModelSpatialRadialProfileDMEinasto& model)
{
    // Copy members. We do not have to push back the members on the parameter
    // stack as this should have been done by init_members() that was called
    // before. Otherwise we would have sigma twice on the stack.
    m_scale_radius  = model.m_scale_radius ;
    m_halo_distance = model.m_halo_distance ;
    m_alpha         = model.m_alpha   ;

    // copy cache values
    m_last_scale_radius = model.m_last_scale_radius ;
    m_mass_radius       = model.m_mass_radius       ;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Radial profile
 *
 * @param[in] theta Angular distance from DMEinasto centre (radians).
 * @return Profile value.
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::profile_value(const double& theta) const
{
    
    // update precomputation cache
    update();
    
    // initialize integral value
    double value = 0.0 ;
    
    // Set up integration limits
    double los_min = m_halo_distance.value() - m_mass_radius ;
    double los_max = m_halo_distance.value() + m_mass_radius ;
    
    // handle case where observer is within halo mass radius
    if ( los_min < 0.0 ) {
      los_min = 0.0 ;
    }
    
    // Set up integral
    halo_kernel_los integrand( m_scale_radius.value(),
                               m_halo_distance.value(),
                               m_alpha.value(),
                               theta ) ;
    GIntegral integral(&integrand) ;
    integral.max_iter( 30 ) ;
    
    // also play with eps() or max_iter(), bounds()
    
    // Set up integration boundaries
    // As there is usually an infinity at the halo center, this splits
    // the integral at the m_halo_distance.
    std::vector<double> bounds ;
    bounds.push_back( los_min ) ;
    bounds.push_back( los_max ) ;
    bounds.push_back( m_halo_distance.value() );
    
    // Compute value
    //value = integral.romberg( los_min, los_max ) ;
    value = integral.romberg( bounds ) ;

    // Return value
    return value;
}

/***********************************************************************//**
 * @brief Kernel for halo density profile squared
 *
 * @param[in] distance from observer to point in space (meters)
 *
 * Computes the value of an einasto halo density profile squared, 
 * at distance l from observer, at angle \f[\theta\f] from the halo center:
 * 
 * \f[
 *    f(\theta, l) = E^{ -\frac{2}{\alpha} \left( g^{\alpha} - 1 \right) }
 * \f]
 *
 * where
 *
 * \f[
 *    g = \frac{ \sqrt{l^2+d^2-2ldCos(\theta)} }{r_s}
 * \f]
 *
 * This profile is detailed in:
 *   Springel et al, 2008
 *   "The Aquarius Project: the subhaloes of galactic haloes"
 *   Mon. Not. R. Astron. Soc. 391, 1685–1711
 *   http://mnras.oxfordjournals.org/content/391/4/1685
 *
 * which cites:
 *   J. Einasto, 1965
 *   "Kinematics and dynamics of stellar systems"
 *   Trudy Inst. Astrofiz. Alma-Ata 5, 87
 *
 * @return unit
 *
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::halo_kernel_los::eval( const double &los )
{
  
  double g = 0.0 ;
  g  = los * los ;
  g += m_halo_distance * m_halo_distance ;
  g -= 2.0 * los * m_halo_distance * gammalib::cosd(m_theta) ;
  g  = sqrt(g) ;
  g /= m_scale_radius ;
  
  double f = 0.0 ;
  f  = pow( g , m_alpha ) ;
  f -= 1.0 ;
  f *= -2.0 / m_alpha ;
  f  = std::exp( f ) ;
  
  return f;

}

/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the m_mass_radius calculation, determining the radius around
 * the halo that contains 99.99% of the mass.  For an einasto halo profile,
 * this is just 10.0 * scale_radius .
 *
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::update() const
{
  
  // Update if scale radius has changed
  if ( m_last_scale_radius != scale_radius() ) {
    
    // Store last values
    m_last_scale_radius = scale_radius() ;
    
    // perform precomputations
    m_mass_radius = 10.0 * scale_radius() ;

  }
}
