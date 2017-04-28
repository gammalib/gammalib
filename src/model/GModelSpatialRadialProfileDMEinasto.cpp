/***************************************************************************
 * GModelSpatialRadialProfileDMEinasto.cpp - Einasto radial profile class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Nathan Kelley-Hoskins                            *
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
//#include <iomanip>

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialProfileDMEinasto g_radial_disk_seed;
const GModelSpatialRegistry               g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ      "GModelSpatialRadialProfileDMEinasto::read(GXmlElement&)"
#define G_WRITE    "GModelSpatialRadialProfileDMEinasto::write(GXmlElement&)"

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
 * @brief Return minimum model radius (in radians)
 *
 * @return Minimum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::theta_min(void) const
{
    
    // update precomputation cache
    update();
    
    // Return value
    return m_theta_min.value();
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::theta_max(void) const
{
    
    // Update precomputation cache
    update();

    // Initialise maximum theta angle
    double theta = 0.0;
    
    // If Earth is within the significant radius, then theta_max must
    // contain the entire profile (180deg) ...
    if (m_halo_distance.value() < m_mass_radius) {
        theta = gammalib::pi;
    }
    
    // ... otherwise if the halo is far enough away (further than the mass
    // radius) then we just need to deal with the angles within the sphere
    // of the significant radius.
    else {
        theta = std::atan(m_mass_radius / m_halo_distance.value());
    }

    // Always chose the lesser of mass_radius theta and theta_max
    double theta_max = m_theta_max.value() * gammalib::deg2rad;
    if (theta > theta_max) {
        theta = theta_max;
    }
    
    // Return maximum theta value
    return theta;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the DMEinasto radial profile model information from an XML element.
 * The XML element shall have the format
 *
 *     <spatialModel type="DMEinastoProfile">
 *       <parameter name="RA"           scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="DEC"          scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleRadius"  scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleDensity" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="HaloDistance" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Alpha"        scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMin"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMax"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="CoreRadius"   scale=.. value=.. min=.. max=.. free=../>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::read(const GXmlElement& xml)
{
    // Read DMEinasto location
    GModelSpatialRadial::read(xml);

    // Read ScaleRadius parameter
    const GXmlElement* par1 = gammalib::xml_get_par(G_READ, xml, "ScaleRadius");
    m_scale_radius.read(*par1);
    
    // Read ScaleDensity parameter
    const GXmlElement* par2 = gammalib::xml_get_par(G_READ, xml, "ScaleDensity");
    m_scale_density.read(*par2);
    
    // Read HaloDistance parameter
    const GXmlElement* par3 = gammalib::xml_get_par(G_READ, xml, "HaloDistance");
    m_halo_distance.read(*par3);

    // Read Alpha parameter
    const GXmlElement* par4 = gammalib::xml_get_par(G_READ, xml, "Alpha");
    m_alpha.read(*par4);

    // Read ThetaMin parameter
    const GXmlElement* par5 = gammalib::xml_get_par(G_READ, xml, "ThetaMin");
    m_theta_min.read(*par5);
    
    // Read ThetaMax parameter
    const GXmlElement* par6 = gammalib::xml_get_par(G_READ, xml, "ThetaMax");
    m_theta_max.read(*par6);
    
    // Read CoreRadius parameter
    const GXmlElement* par7 = gammalib::xml_get_par(G_READ, xml, "CoreRadius");
    m_core_radius.read(*par7);

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
 *       <parameter name="RA"           scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="DEC"          scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleRadius"  scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleDensity" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="HaloDistance" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Alpha"        scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMin"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMax"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="CoreRadius"   scale=.. value=.. min=.. max=.. free=../>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::write(GXmlElement& xml) const
{
    // Write DMEinasto location
    GModelSpatialRadial::write(xml);

    // Write ScaleRadius parameter
    GXmlElement* par1 = gammalib::xml_need_par(G_WRITE, xml, "ScaleRadius");
    m_scale_radius.write(*par1);
    
    // Write ScaleDensity parameter
    GXmlElement* par2 = gammalib::xml_need_par(G_WRITE, xml, "ScaleDensity");
    m_scale_density.write(*par2);

    // Write HaloDistance parameter
    GXmlElement* par3 = gammalib::xml_need_par(G_WRITE, xml, "HaloDistance");
    m_halo_distance.write(*par3);

    // Write Alpha parameter
    GXmlElement* par4 = gammalib::xml_need_par(G_WRITE, xml, "Alpha");
    m_alpha.write(*par4);

    // Write ThetaMin parameter
    GXmlElement* par5 = gammalib::xml_need_par(G_WRITE, xml, "ThetaMin");
    m_theta_min.write(*par5);
    
    // Write ThetaMax parameter
    GXmlElement* par6 = gammalib::xml_need_par(G_WRITE, xml, "ThetaMax");
    m_theta_max.write(*par6);
    
    // Write CoreRadius parameter
    GXmlElement* par7 = gammalib::xml_need_par(G_WRITE, xml, "CoreRadius");
    m_core_radius.write(*par7);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness.
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
    m_scale_radius.name("ScaleRadius");
    m_scale_radius.unit("kpc");
    m_scale_radius.value(21.5);       // default to GC scale radius
    m_scale_radius.min(1.0e-6);       // arbitrarily chosen :/
    m_scale_radius.free();
    m_scale_radius.scale(1.0);
    m_scale_radius.gradient(0.0);
    m_scale_radius.has_grad(false);   // Radial components never have gradients
    
    // Initialise scale density
    m_scale_density.clear();
    m_scale_density.name("ScaleDensity");
    m_scale_density.unit("GeV/cm^3");
    m_scale_density.value(0.2);       // GeV/cm3, default to GC scale density
    m_scale_density.min(1.0e-6);
    m_scale_density.fix();
    m_scale_density.scale(1.0);
    m_scale_density.gradient(0.0);
    m_scale_density.has_grad(false);  // Radial components never have gradients

    // Initialise halo distance
    m_halo_distance.clear();
    m_halo_distance.name("HaloDistance");
    m_halo_distance.unit("kpc");
    m_halo_distance.value(7.94);     // default to GC halo distance
    m_halo_distance.min(1.0e-6);     // arbitrarily chosen
    m_halo_distance.free();
    m_halo_distance.scale(1.0);
    m_halo_distance.gradient(0.0);
    m_halo_distance.has_grad(false); // Radial components never have gradients

    // Initialise alpha
    m_alpha.clear();
    m_alpha.name("Alpha");
    m_alpha.unit("unitless");
    m_alpha.value(0.17);             // default to GC alpha
    m_alpha.min(0.01);               // arbitrarily chosen
    m_alpha.max(10.0);               // arbitrarily chosen
    m_alpha.free();
    m_alpha.scale(1.0);
    m_alpha.gradient(0.0);
    m_alpha.has_grad(false);         // Radial components never have gradients

    // Initialise theta max 
    m_theta_min.clear();
    m_theta_min.name("ThetaMin");
    m_theta_min.unit("degrees");
    m_theta_min.value(180.0);        // can only go from halo center to opposite halo center
    m_theta_min.min(1.0e-6);         // arbitrarily chosen
    m_theta_min.fix();               // should always be fixed!
    m_theta_min.scale(1.0);
    m_theta_min.gradient(0.0);
    m_theta_min.has_grad(false);     // Radial components never have gradients

    // Initialise theta max 
    m_theta_max.clear();
    m_theta_max.name("ThetaMax");
    m_theta_max.unit("degrees");
    m_theta_max.value(1.0e-6);       // can only go from halo center to opposite halo center
    m_theta_max.min(1.0e-10);        // arbitrarily chosen
    m_theta_max.fix();               // should always be fixed!
    m_theta_max.scale(1.0);
    m_theta_max.gradient(0.0);
    m_theta_max.has_grad(false);     // Radial components never have gradients
    
    // Initialise core radius
    m_core_radius.clear();
    m_core_radius.name("CoreRadius");
    m_core_radius.unit("kpc");
    m_core_radius.value(0.5);        // example: galactic center core
    m_core_radius.min(0.0);
    m_core_radius.fix();
    m_core_radius.scale(1.0);
    m_core_radius.gradient(0.0);
    m_core_radius.has_grad(false);

    // Set parameter pointer(s)
    m_pars.push_back(&m_scale_radius);
    m_pars.push_back(&m_scale_density);
    m_pars.push_back(&m_halo_distance);
    m_pars.push_back(&m_alpha);
    m_pars.push_back(&m_theta_min);
    m_pars.push_back(&m_theta_max);
    m_pars.push_back(&m_core_radius);
    
    // Initialize precomputation cache. Note that zero values flag
    // uninitialised, as a zero radius is not meaningful
    m_last_scale_radius     = 0.0;
    m_last_scale_density    = 0.0;
    m_mass_radius           = 0.0;
    m_scale_density_squared = 0.0;

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
    // before.
    m_scale_radius  = model.m_scale_radius;
    m_scale_density = model.m_scale_density;
    m_halo_distance = model.m_halo_distance;
    m_alpha         = model.m_alpha;
    m_theta_min     = model.m_theta_min;
    m_theta_max     = model.m_theta_max;
    m_core_radius   = model.m_core_radius;

    // copy cache values
    m_last_scale_radius     = model.m_last_scale_radius;
    m_last_scale_density    = model.m_last_scale_density;
    m_mass_radius           = model.m_mass_radius;
    m_scale_density_squared = model.m_scale_density_squared;
    
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
    // Update precomputation cache
    update();
    
    // Initialize integral value
    double value = 0.0;
    
    // Set up integration limits
    double los_min = m_halo_distance.value() - m_mass_radius;
    double los_max = m_halo_distance.value() + m_mass_radius;
    
    // handle case where observer is within halo mass radius
    if (los_min < 0.0) {
        los_min = 0.0;
    }
    
    // Set up integral
    halo_kernel_los integrand(m_scale_radius.value(),
                              m_halo_distance.value(),
                              m_alpha.value(),
                              theta,
                              m_core_radius.value());
    GIntegral integral(&integrand);
    integral.max_iter(30);
    
    // Set up integration boundaries. As there is usually an infinity at the
    // halo center, this splits the integral at the m_halo_distance.
    std::vector<double> bounds;
    bounds.push_back(los_min);
    bounds.push_back(los_max);
    bounds.push_back(m_halo_distance.value());
    
    // Compute value
    value = integral.romberg(bounds);

    // Multiply in the density^2
    value *= m_scale_density_squared;

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
    // Calculate the scale distance g, the ( distance from integration point
    // to the halo center ) divided by ( the halo scale radius )
    
    // First calculate the distance of the integration point from the halo
    // center via the law of cosines
    double g = los * los;
    g       += m_halo_distance * m_halo_distance;
    g       -= 2.0 * los * m_halo_distance * std::cos(m_theta);
    g        = std::sqrt(g);
    
    // If we have a core radius specified, all halo values inside this core
    // radius should be the same as at the core radius itself.
    if (g < m_core_radius) {
        g = m_core_radius;
    }
    
    // Finish scaling the integration point by the halo's scale radius
    g /= m_scale_radius;
  
    // Calculate the halo density f at the scale distance g
    double f  = std::pow(g, m_alpha);
    f        -= 1.0;
    f        *= -2.0 / m_alpha;
    f         = std::exp(f);
  
    // Squared, for annihilating DM, would just be f if it was decaying DM
    f = f * f;

    // Return function value
    return f;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the m_mass_radius calculation, determining the radius around
 * the halo that contains 99.99% of the mass. For an Einasto halo profile,
 * this is just 10.0 * scale_radius .
 ***************************************************************************/
void GModelSpatialRadialProfileDMEinasto::update() const
{
    // Update if scale radius has changed
    if (m_last_scale_radius  != scale_radius() ||
        m_last_scale_density != scale_density()) {
    
        // Store last values
        m_last_scale_radius = scale_radius();
        m_last_scale_density = scale_density();
    
        // Perform precomputations
        m_mass_radius           = 10.0 * scale_radius();
        m_scale_density_squared = scale_density() * scale_density();

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Calculate halo mass density
 *
 * @param[in] radius from halo center (kpc)
 * @return Halo mass density.
 *
 * Calculates the halo's mass density at a given radial distance from the
 * halo center.
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::mass_density( const double& radius ) const
{
    // Set-up kernel
    halo_kernel_los halo_shape(m_scale_radius.value(),
                               m_halo_distance.value(),
                               m_alpha.value(),
                               0.0,
                               m_core_radius.value());
  
  
    // eval produces a unitless density^2, so we must take its square root
    double density = std::sqrt(halo_shape.eval(m_halo_distance.value() + radius));
  
    // Multiply in the missing scale density
    density *= m_scale_density.value();

    // Return halo mass density
    return density;
}


/***********************************************************************//**
 * @brief Calculate J-factor
 *
 * @param[in] angle from halo center (radians)
 * @return J-factor.
 *
 * Calculates the halo's J-factor at an angle from the halo center.
 ***************************************************************************/
double GModelSpatialRadialProfileDMEinasto::jfactor(const double& angle) const
{
    // Integration settings
    double minradian = 0.0;
    int    npoints   = 200;

    // Initialize other variables
    double jfactor   = 0.0 ;
    double dr        = (angle - minradian) / npoints;
    double r         = 0.0;

    // Loop over different radii in the profile
    for (int i = 0; i < npoints; ++i) {
        r      = minradian + (i * dr);
        jfactor += profile_value(r) * r * dr;
    }
  
    // J-factor = 2 * pi * Int[ profile(r) * r * dr , {r,minradian,angle} ]
    jfactor *= gammalib::twopi;
  
    // Return J-factor
    return jfactor;
}
