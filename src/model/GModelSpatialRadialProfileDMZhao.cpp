/***************************************************************************
 * GModelSpatialRadialProfileDMZhao.cpp - DMZhao radial profile class *
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
 * @file GModelSpatialRadialProfileDMZhao.cpp
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
#include "GModelSpatialRadialProfileDMZhao.hpp"
#include "GModelSpatialRegistry.hpp"
#include <iomanip>

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialProfileDMZhao g_radial_disk_seed;
const GModelSpatialRegistry           g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ   "GModelSpatialRadialProfileDMZhao::read(GXmlElement&)"
#define G_WRITE  "GModelSpatialRadialProfileDMZhao::write(GXmlElement&)"

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
 * Constructs empty radial DMZhao profile
 ***************************************************************************/
GModelSpatialRadialProfileDMZhao::GModelSpatialRadialProfileDMZhao(void) :
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
 * Constructs radial DMZhao profile model by extracting information from
 * an XML element. See the read() method for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialProfileDMZhao::GModelSpatialRadialProfileDMZhao(const GXmlElement& xml) :
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
 * @param[in] model Radial DMZhao profile model.
 *
 * Copies radial DMZhao profile model from another radial profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMZhao::GModelSpatialRadialProfileDMZhao(const GModelSpatialRadialProfileDMZhao& model) :
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
 * Destructs radial DMZhao profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMZhao::~GModelSpatialRadialProfileDMZhao(void)
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
 * @param[in] model Radial DMZhao profile model.
 * @return Radial DMZhao profile model.
 *
 * Assigns radial DMZhao profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMZhao& GModelSpatialRadialProfileDMZhao::operator=(const GModelSpatialRadialProfileDMZhao& model)
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
 * @brief Clear radial DMZhao profile model
 *
 * Clears radial DMZhao profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::clear(void)
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
 * @brief Clone radial DMZhao profile model
 *
 * @return Pointer to deep copy of radial DMZhao profile model.
 *
 * Returns a deep copy of the radial DMZhao profile model.
 ***************************************************************************/
GModelSpatialRadialProfileDMZhao* GModelSpatialRadialProfileDMZhao::clone(void) const
{
    // Clone radial disk model
    return new GModelSpatialRadialProfileDMZhao(*this);
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::theta_max(void) const
{
    
    // update precomputation cache
    update();
    
    double theta = 0.0 ;
    
    // if earth is within the significant radius, then theta_max must
    // contain the entire profile (180deg)
    if ( m_halo_distance.value() < m_mass_radius ) {
      theta = gammalib::pi ;
    
    // if the halo is far enough away (further than the significant radius)
    // then we just need to deal with the angles within the sphere of the
    // significant radius.
    } else {
      theta = std::atan( m_mass_radius / m_halo_distance.value() ) ;
    }
    
    // Return value
    return theta ;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the DMZhao radial profile model information from an XML element.
 * The XML element shall have either the format 
 *
 *     <spatialModel type="DMZhaoProfile">
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
 *     <spatialModel type="DMZhaoProfile">
 *       <parameter name="GLON"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"  scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *       <parameter name="Scale Radius"  scale="1.0" value="21.5" min="1.0e-6" free="1"/>
 *       <parameter name="Halo Distance" scale="1.0" value="7.94" min="1.0e-6" free="1"/>
 *       <parameter name="Alpha"         scale="1.0" value="0.17" min="0.01"   max="10" free="1"/>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::read(const GXmlElement& xml)
{
    // Read halo sky coordinate
    GModelSpatialRadial::read(xml);

    // Read Scale Radius parameter
    const GXmlElement* par1 = gammalib::xml_get_par(G_READ, xml, "Scale Radius");
    m_scale_radius.read(*par1);
    
    // Read Halo Distance parameter
    const GXmlElement* par2 = gammalib::xml_get_par(G_READ, xml, "Halo Distance");
    m_halo_distance.read(*par2);

    // Read Alpha parameter
    const GXmlElement* par3 = gammalib::xml_get_par(G_READ, xml, "Alpha");
    m_alpha.read(*par3);

    // Read Beta parameter
    const GXmlElement* par4 = gammalib::xml_get_par(G_READ, xml, "Beta");
    m_beta.read(*par4);

    // Read Gamma parameter
    const GXmlElement* par5 = gammalib::xml_get_par(G_READ, xml, "Gamma");
    m_gamma.read(*par5);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Writes the DMZhao radial profile model information into an XML element.
 * The XML element will have the format 
 *
 *     <spatialModel type="DMZhaoProfile">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::write(GXmlElement& xml) const
{
    // Write DMZhao location
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

    // Write Beta parameter
    GXmlElement* par4 = gammalib::xml_need_par(G_WRITE, xml, "Beta");
    m_beta.write(*par4);

    // Write Gamma parameter
    GXmlElement* par5 = gammalib::xml_need_par(G_WRITE, xml, "Gamma");
    m_gamma.write(*par5);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialProfileDMZhao::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialProfileDMZhao ===");

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
void GModelSpatialRadialProfileDMZhao::init_members(void)
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
    m_alpha.value(1.0);  // default to NFW alpha
    m_alpha.min(0.01);   // arbitrarily chosen
    m_alpha.max(10.0);   // arbitrarily chosen
    m_alpha.free();
    m_alpha.scale(1.0);
    m_alpha.gradient(0.0);
    m_alpha.has_grad(false);  // Radial components never have gradients

    // Initialise 
    m_beta.clear();
    m_beta.name("Beta");
    m_beta.unit("unitless");
    m_beta.value(3.0);  // default to NFW beta
    m_beta.min(0.01);   // arbitrarily chosen
    m_beta.max(10.0);   // arbitrarily chosen
    m_beta.free();
    m_beta.scale(1.0);
    m_beta.gradient(0.0);
    m_beta.has_grad(false);  // Radial components never have gradients

    // Initialise 
    m_gamma.clear();
    m_gamma.name("Gamma");
    m_gamma.unit("unitless");
    m_gamma.value(1.0);  // default to NFW gamma
    m_gamma.min(0.01);   // arbitrarily chosen
    m_gamma.max(10.0);   // arbitrarily chosen
    m_gamma.free();
    m_gamma.scale(1.0);
    m_gamma.gradient(0.0);
    m_gamma.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_scale_radius );
    m_pars.push_back(&m_halo_distance);
    m_pars.push_back(&m_alpha        );
    m_pars.push_back(&m_beta         );
    m_pars.push_back(&m_gamma        );

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial DMZhao model.
 *
 * Copies class members from another radial profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::copy_members(const GModelSpatialRadialProfileDMZhao& model)
{
    // Copy members. We do not have to push back the members on the parameter
    // stack as this should have been done by init_members() that was called
    // before. Otherwise we would have sigma twice on the stack.
    m_scale_radius  = model.m_scale_radius  ;
    m_halo_distance = model.m_halo_distance ;
    m_alpha         = model.m_alpha         ;
    m_beta          = model.m_beta          ;
    m_gamma         = model.m_gamma         ;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Radial profile
 *
 * @param[in] theta Angular distance from DMZhao centre (radians).
 * @return Profile value.
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::profile_value(const double& theta) const
{
    
    // update precomputation cache
    update();
    
    // initialize integral value
    double value   = 0.0 ;
    
    // Set up integration limits
    double los_min = m_halo_distance.value() - m_mass_radius ;
    double los_max = m_halo_distance.value() + m_mass_radius ;
    
    // Set up integral
    halo_kernel_los integrand( m_scale_radius.value(),
                               m_halo_distance.value(),
                               m_alpha.value(),
                               m_beta.value(),
                               m_gamma.value(),
                               theta ) ;
    GIntegral integral(&integrand) ;
    
    // Set up integration boundaries
    // As there is usually an infinity at the halo center, this splits
    // the integral at the m_halo_distance.
    std::vector<double> bounds ;
    bounds.push_back( los_min ) ;
    bounds.push_back( los_max ) ;
    bounds.push_back( m_halo_distance.value() );

    // Compute value
    value = integral.romberg( bounds ) ;
    
    //std::cout << "  theta=" << theta << "  profile_value=" << value << std::endl;

    // Return value
    return value;
}

/***********************************************************************//**
 * @brief Kernel for zhao halo density profile squared
 *
 * @param[in] distance from observer to point in space (meters)
 *
 * Computes the value of a zhao halo density profile squared, 
 * at distance l from observer, at angle \f[\theta\f] from the halo center,
 * with a halo posessing a scale radius of \f[r_s\f] :
 * 
 * \f[
 *    f(\theta, l) = \frac{1}{ g^{\gamma} {\left( g^{\alpha} + 1 \right)}^{ \frac{\beta-\gamma}{\alpha} } }
 * \f]
 *
 * where
 *
 * \f[
 *    g = \frac{ \sqrt{l^2+d^2-2ldCos(\theta)} }{r_s}
 * \f]
 *
 * \f[ \beta \f] is the slope of the density at radii much bigger than \f[r_s\f]
 * \f[ \gamma \f] is the slope of the density profile at radii much smaller than \f[r_s\f]
 * \f[ \alpha \f] governs how large the transition region is between these two slopes.
 *   larger \f[ \alpha \f] means a smaller transition region.
 *
 * This profile is detailed in:
 *   Zhao, 1996
 *   "Analytical models for galactic nuclei"
 *   Monthly Notices of the Royal Astronomical Society, 278, 488-49
 *   http://mnras.oxfordjournals.org/content/278/2/488.short
 *
 * @return unit
 *
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::halo_kernel_los::eval( const double &los )
{
  
  double g = 0.0 ;
  g  = los * los ;
  g += m_halo_distance * m_halo_distance ;
  g -= 2.0 * los * m_halo_distance * std::cos(m_theta) ;
  g  = sqrt(g) ;
  g /= m_scale_radius ;
  
  double f = 0.0 ;
  f  = pow( g , m_alpha ) ;
  f += 1 ;
  f  = pow( f, (m_beta-m_gamma)/m_alpha ) ;
  f *= pow( g, m_gamma ) ;
  f = 1 / f ;

  // squared, for annihilating dm
  // would just be f if it was decaying dm
  f = f * f ;
  
  //std::cout << std::setprecision(10) << "kern_los::eval  los=" << los << "  hd=" << m_halo_distance << "  theta=" << m_theta << "  alpha=" << m_alpha << "  beta=" << m_beta << "  gamma=" << m_gamma << "  g=" << g << "  f=" << f << std::endl ;
  
  return f;

}

/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the m_mass_radius calculation, determining the radius around
 * the halo that contains 99.99% of the mass.  For a zhao halo profile,
 * this is just 10.0 * scale_radius .
 *
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::update() const
{
  
  // Update if scale radius has changed
  if ( m_last_scale_radius != scale_radius() ) {
    
    // Store last values
    m_last_scale_radius = scale_radius() ;
    
    // perform precomputations
    m_mass_radius = 10.0 * scale_radius() ;

  }
}

double GModelSpatialRadialProfileDMZhao::prof_val( const double& theta )
{
  return this->profile_value(theta) ;
}
