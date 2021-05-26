/***************************************************************************
 *    GModelSpatialRadialProfileDMZhao.cpp - Zhao radial profile class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2021 by Nathan Kelley-Hoskins                       *
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

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialProfileDMZhao g_radial_disk_seed;
const GModelSpatialRegistry            g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ         "GModelSpatialRadialProfileDMZhao::read(GXmlElement&)"
#define G_WRITE       "GModelSpatialRadialProfileDMZhao::write(GXmlElement&)"

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
 * @brief Return minimum model radius (in radians)
 *
 * @return Minimum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::theta_min(void) const
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
double GModelSpatialRadialProfileDMZhao::theta_max(void) const
{
    // Update precomputation cache
    update();

    // Initialize maximum theta angle
    double theta = 0.0;

    // If Earth is within the significant radius, then theta_max must
    // contain the entire profile (180deg) ...
    if (m_halo_distance.value() < m_mass_radius) {
        theta = gammalib::pi;
    }

    // ... otherwise, if the halo is far enough away (further than the
    // significant radius) then we just need to deal with the angles within
    // the sphere of the significant radius.
    else {
        theta = std::atan(m_mass_radius / m_halo_distance.value());
    }

    // Always chose the lesser of ( mass_radius theta, theta_max )
    double theta_max = m_theta_max.value() * gammalib::deg2rad;
    if (theta_max < theta) {
        theta = theta_max;
    }

    // Return value
    return theta;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the DMZhao radial profile model information from an XML element.
 * The XML element shall have the format
 *
 *     <spatialModel type="DMZhaoProfile">
 *       <parameter name="RA"           scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="DEC"          scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleRadius"  scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleDensity" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="HaloDistance" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Alpha"        scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Beta"         scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Gamma"        scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMin"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMax"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="CoreRadius"   scale=.. value=.. min=.. max=.. free=../>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 11);

    // Read halo sky coordinate
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

    // Read Beta parameter
    const GXmlElement* par5 = gammalib::xml_get_par(G_READ, xml, "Beta");
    m_beta.read(*par5);

    // Read Gamma parameter
    const GXmlElement* par6 = gammalib::xml_get_par(G_READ, xml, "Gamma");
    m_gamma.read(*par6);

    // Read ThetaMin parameter
    const GXmlElement* par7 = gammalib::xml_get_par(G_READ, xml, "ThetaMin");
    m_theta_min.read(*par7);

    // Read ThetaMax parameter
    const GXmlElement* par8 = gammalib::xml_get_par(G_READ, xml, "ThetaMax");
    m_theta_max.read(*par8);
    
    // Read CoreRadius parameter
    const GXmlElement* par9 = gammalib::xml_get_par(G_READ, xml, "CoreRadius");
    m_core_radius.read(*par9);

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
 *       <parameter name="RA"           scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="DEC"          scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleRadius"  scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ScaleDensity" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="HaloDistance" scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Alpha"        scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Beta"         scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="Gamma"        scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMin"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="ThetaMax"     scale=.. value=.. min=.. max=.. free=../>
 *       <parameter name="CoreRadius"   scale=.. value=.. min=.. max=.. free=../>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Write DMZhao location
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

    // Write Beta parameter
    GXmlElement* par5 = gammalib::xml_need_par(G_WRITE, xml, "Beta");
    m_beta.write(*par5);

    // Write Gamma parameter
    GXmlElement* par6 = gammalib::xml_need_par(G_WRITE, xml, "Gamma");
    m_gamma.write(*par6);

    // Write ThetaMax parameter
    GXmlElement* par7 = gammalib::xml_need_par(G_WRITE, xml, "ThetaMin");
    m_theta_min.write(*par7);

    // Write ThetaMax parameter
    GXmlElement* par8 = gammalib::xml_need_par(G_WRITE, xml, "ThetaMax");
    m_theta_max.write(*par8);
    
    // Write CoreRadius parameter
    GXmlElement* par9 = gammalib::xml_need_par(G_WRITE, xml, "CoreRadius");
    m_core_radius.write(*par9);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness.
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
    // Initialise model type
    m_type = "DMZhaoProfile";

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
    m_scale_density.value(0.2);       // GeV/cm^3 default to GC scale density
    m_scale_density.free();
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
    m_alpha.value(1.0);              // default to NFW alpha
    m_alpha.min(0.01);               // arbitrarily chosen
    m_alpha.max(10.0);               // arbitrarily chosen
    m_alpha.free();
    m_alpha.scale(1.0);
    m_alpha.gradient(0.0);
    m_alpha.has_grad(false);         // Radial components never have gradients

    // Initialise 
    m_beta.clear();
    m_beta.name("Beta");
    m_beta.unit("unitless");
    m_beta.value(3.0);               // default to NFW beta
    m_beta.min(0.01);                // arbitrarily chosen
    m_beta.max(10.0);                // arbitrarily chosen
    m_beta.free();
    m_beta.scale(1.0);
    m_beta.gradient(0.0);
    m_beta.has_grad(false);          // Radial components never have gradients

    // Initialise 
    m_gamma.clear();
    m_gamma.name("Gamma");
    m_gamma.unit("unitless");
    m_gamma.value(1.0);              // default to NFW gamma
    m_gamma.min(0.01);               // arbitrarily chosen
    m_gamma.max(10.0);               // arbitrarily chosen
    m_gamma.free();
    m_gamma.scale(1.0);
    m_gamma.gradient(0.0);
    m_gamma.has_grad(false);         // Radial components never have gradients

    // Initialise theta max 
    m_theta_min.clear();
    m_theta_min.name("ThetaMin");
    m_theta_min.unit("degrees");
    m_theta_min.value(1.0e-6);      // can only go from halo center to opposite halo center
    m_theta_min.min(1.0e-10);       // arbitrarily chosen, some halos don't converge at theta=0.0
    m_theta_min.fix();              // should always be fixed!
    m_theta_min.scale(1.0);
    m_theta_min.gradient(0.0);
    m_theta_min.has_grad(false);    // Radial components never have gradients

    // Initialise theta max 
    m_theta_max.clear();
    m_theta_max.name("ThetaMax");
    m_theta_max.unit("degrees");
    m_theta_max.value(180.0);       // can only go from halo center to opposite halo center
    m_theta_max.min(1.0e-6);        // arbitrarily chosen
    m_theta_max.fix();              // should always be fixed!
    m_theta_max.scale(1.0);
    m_theta_max.gradient(0.0);
    m_theta_max.has_grad(false);   // Radial components never have gradients
    
    // Initialise core radius
    m_core_radius.clear();
    m_core_radius.name("CoreRadius");
    m_core_radius.unit("kpc");
    m_core_radius.value(10.0);
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
    m_pars.push_back(&m_beta);
    m_pars.push_back(&m_gamma);
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
 * @param[in] model Radial DMZhao model.
 *
 * Copies class members from another radial profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::copy_members(const GModelSpatialRadialProfileDMZhao& model)
{
    // Copy members. We do not have to push back the members on the parameter
    // stack as this should have been done by init_members() that was called
    // before.
    m_type          = model.m_type;   // Needed to conserve model type
    m_scale_radius  = model.m_scale_radius;
    m_scale_density = model.m_scale_density;
    m_halo_distance = model.m_halo_distance;
    m_alpha         = model.m_alpha;
    m_beta          = model.m_beta;
    m_gamma         = model.m_gamma;
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
    // Update precomputation cache
    update();

    // Set up integration limits
    double los_min = m_halo_distance.value() - m_mass_radius;
    double los_max = m_halo_distance.value() + m_mass_radius;

    // Set up integral
    halo_kernel_los integrand(m_scale_radius.value(),
                              m_halo_distance.value(),
                              m_alpha.value(),
                              m_beta.value(),
                              m_gamma.value(),
                              theta,
                              m_core_radius.value());
    GIntegral integral(&integrand);
    integral.max_iter(25);

    // Set up integration boundaries
    // As there is usually an infinity at the halo center, this splits
    // the integral at the m_halo_distance.
    std::vector<double> bounds;
    bounds.push_back(los_min);
    bounds.push_back(los_max);
    bounds.push_back(m_halo_distance.value());

    // Compute value
    double value = integral.romberg(bounds);

    // apply scale density squared
    // TODO: must be multiplied by the particle physics factor
    value *= m_scale_density_squared;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for zhao halo density profile squared
 *
 * @param[in] los Distance from observer to point in space (meters)
 * @return Halo density.
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
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::halo_kernel_los::eval(const double &los)
{
    // Calculate the scale distance g, the ( distance from integration point
    // to the halo center ) divided by ( the halo scale radius )
    
    // First calculate the distance of the integration point from the halo
    // center via the law of cosines
    double g  = los * los;
    g        += m_halo_distance * m_halo_distance;
    g        -= 2.0 * los * m_halo_distance * std::cos(m_theta);
    g         = std::sqrt(g);
    
    // If we have a core radius specified, all halo values inside this core
    // radius should be the same as at the core radius itself.
    if (g < m_core_radius) {
        g = m_core_radius;
    }
    
    // Finish scaling the integration point by the halo's scale radius
    g /= m_scale_radius ;
  
    // apply the zhao profile from the scaled distance to the halo center
    double f  = std::pow(g, m_alpha);
    f        += 1.0;
    f         = std::pow(f, (m_beta-m_gamma)/m_alpha);
    f        *= std::pow(g, m_gamma);
    f         = 1.0 / f;

    // Squared, for annihilating DM, would just be f if it was decaying DM
    f = f * f;

    // Return function value
    return f;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the m_mass_radius calculation, determining the radius around
 * the halo that contains 99.99% of the mass. For a Zhao halo profile,
 * this is just 10.0 * scale_radius .
 ***************************************************************************/
void GModelSpatialRadialProfileDMZhao::update() const
{
    // Update if scale radius has changed
    if (m_last_scale_radius  != scale_radius() ||
        m_last_scale_density != scale_density() ) {

        // Store last values
        m_last_scale_radius = scale_radius();
        m_last_scale_density = scale_density();

        // perform precomputations
        m_mass_radius = 10.0 * scale_radius();
        m_scale_density_squared = scale_density() * scale_density() ;

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Calculate Halo Mass Density
 *
 * @param[in] radius Distance from halo center (kpc).
 * @return Halo mass density.
 *
 * Calculates the halo's mass density at a given radial distance from the halo
 * center.
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::mass_density(const double& radius) const
{
    // Initialize halo kernel with stored values. We use theta=0 to effectivly
    // sample the halo radially
    halo_kernel_los halo_shape(m_scale_radius.value(),
                               m_halo_distance.value(),
                               m_alpha.value(),
                               m_beta.value(),
                               m_gamma.value(),
                               0.0,
                               m_core_radius.value());
  
    // eval produces a unitless density^2, so we must take its square root
    double density = std::sqrt(halo_shape.eval(m_halo_distance.value() + radius));
    
    // Multiply in the missing scale density
    density *= m_scale_density.value();
    
    return density ;
}


/***********************************************************************//**
 * @brief Calculate J-factor
 *
 * @param[in] angle from halo center (radians)
 * @return J-factor
 *
 * Calculates the halo's J-factor at an angle from the halo center.
 ***************************************************************************/
double GModelSpatialRadialProfileDMZhao::jfactor(const double& angle) const
{
    // Integration settings
    const double minradian = 0.0;
    const int    npoints   = 200;
  
    // Initialize other variables
    double jfactor = 0.0;
    double dr      = (angle - minradian) / npoints;
  
    // Loop over different radii in the profile
    for (int i = 0; i < npoints; ++i) {
        double r = minradian + (i * dr);
        jfactor += profile_value(r) * r * dr;
    }
  
    // J-factor = 2 * pi * Int[ profile(r) * r * dr , {r,minradian,angle} ]
    jfactor *= gammalib::twopi;

    // Return J-factor
    return jfactor;
}
