/***************************************************************************
 *       GCTAModelRadialProfile.cpp - Radial Profile CTA model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialProfile.cpp
 * @brief Radial Profile model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GIntegral.hpp"
#include "GCTAObservation.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAModelRadialProfile.hpp"
#include "GCTAModelRadialRegistry.hpp"
#include "GCTAModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRadialProfile   g_cta_radial_profile_seed;
const GCTAModelRadialRegistry  g_cta_radial_profile_registry(&g_cta_radial_profile_seed);
const GCTAModelSpatialRegistry g_cta_radial_profile_spatial_registry(&g_cta_radial_profile_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                   "GCTAModelRadialProfile::read(GXmlElement&)"
#define G_WRITE                 "GCTAModelRadialProfile::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_DEBUG_MC                                     //!< Debug MC method

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelRadialProfile::GCTAModelRadialProfile(void) : GCTAModelRadial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] width Width.
 * @param[in] core Core size.
 * @param[in] tail Tail size.
 ***************************************************************************/
GCTAModelRadialProfile::GCTAModelRadialProfile(const double& width,
                                               const double& core,
                                               const double& tail) : GCTAModelRadial()
{
    // Initialise members
    init_members();

    // Assign values
    this->width(width);
    this->core(core);
    this->tail(tail);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a radial Profile model by extracting information from
 * an XML element. See GCTAModelRadialProfile::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GCTAModelRadialProfile::GCTAModelRadialProfile(const GXmlElement& xml)
                                                          : GCTAModelRadial()
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
 * @param[in] model Radial Profile model.
 ***************************************************************************/
GCTAModelRadialProfile::GCTAModelRadialProfile(const GCTAModelRadialProfile& model)
                                                     : GCTAModelRadial(model)
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
GCTAModelRadialProfile::~GCTAModelRadialProfile(void)
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
 * @param[in] model Radial Profile model.
 ***************************************************************************/
GCTAModelRadialProfile& GCTAModelRadialProfile::operator=(const GCTAModelRadialProfile& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GCTAModelRadial::operator=(model);

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
 * @brief Clear instance
***************************************************************************/
void GCTAModelRadialProfile::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAModelRadial::free_members();

    // Initialise members
    this->GCTAModelRadial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAModelRadialProfile* GCTAModelRadialProfile::clone(void) const
{
    return new GCTAModelRadialProfile(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] offset Offset angle [degrees].
 * @param[in] gradients Compute gradients?
 * @return Function value
 *
 * Evaluates the Profile model for a given offset. The Profile model is
 * defined as
 * \f[f(\theta) = (1 + (\theta/c_0)^{c_1})^{-c_2/c_1}\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees),
 * \f$c_0\f$ is the width of the profile (width),
 * \f$c_1\f$ is the width of the central plateau (core), and
 * \f$c_2\f$ is the size of the tail (tail).
 *
 * Note that no analytical partial derivatives are implemented for this
 * function.
 ***************************************************************************/
double GCTAModelRadialProfile::eval(const double& offset,
                                    const bool&   gradients) const
{
    // Compute value
    double arg   = 1.0 + std::pow(offset / width(), core());
    double value = std::pow(arg, -tail()/core());

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAModelRadialProfile::eval";
        std::cout << "(offset=" << offset << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", width=" << width();
        std::cout << ", core=" << core();
        std::cout << ", tail=" << tail();
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in,out] ran Random number generator.
 * @return CTA instrument direction.
 *
 * Draws an arbitrary CTA instrument position from
 * \f[f(\theta) = \sin(\theta) (1 + (\theta/c_0)^{c_1})^{-c_2/c_1}\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees),
 * \f$c_0\f$ is the width of the profile (width),
 * \f$c_1\f$ is the width of the central plateau (core), and
 * \f$c_2\f$ is the size of the tail (tail).
 * using the rejection method.
 *
 * The maximum offset angle that is thrown by this method is fixed to 10
 * degrees.
 *
 * @todo Method can be optimised by using a random deviate of sin instead
 *       of the uniform random deviate which leads to many unnecessary
 *       rejections.
 ***************************************************************************/
GCTAInstDir GCTAModelRadialProfile::mc(GRan& ran) const
{
    // Set constants
    const double offset_max = 10.0;
    const double u_max      = std::sin(offset_max * gammalib::deg2rad);

    // Debug option: initialise number if samples thrown for one value
    #if defined(G_DEBUG_MC)
    int n_samples = 0;
    #endif
    
    // Simulate offset from photon arrival direction until we're not
    // rejected anymore
    double value  = 0.0;
    double u      = 1.0;
    double offset = 0.0;
    do {
        // Throw offset angle
        offset = ran.uniform() * offset_max;
        
        // Compute function value at this offset angle
        value  = std::sin(offset * gammalib::deg2rad) * eval(offset);
        
        // Throw value for rejection method
        u = ran.uniform() * u_max;

        // Debug option: update number of samples
        #if defined(G_DEBUG_MC)
        n_samples++;
        #endif
        
    } while (u > value);

    // Debug option: print number if samples thrown for one value
    #if defined(G_DEBUG_MC)
    std::cout << "#=" << n_samples << " ";
    #endif

    // Simulate azimuth angle
    double phi = 360.0 * ran.uniform();

    // Convert from degrees to radians
    offset *= gammalib::deg2rad;
    phi    *= gammalib::deg2rad;

    // Compute DETX and DETY coordinates
    double detx(0.0);
    double dety(0.0);
	if (offset > 0.0 ) {
		detx = offset * std::cos(phi);
		dety = offset * std::sin(phi);
	}

    // Set instrument direction
    GCTAInstDir dir(detx, dety);

    // Return instrument direction
    return dir;
}


/***********************************************************************//**
 * @brief Return maximum function value for Monte Carlo simulations
 *
 * @param[in] obs CTA Observation (not used).
 * @return Maximum function value for Monte Carlo simulations.
 ***************************************************************************/
double GCTAModelRadialProfile::mc_max_value(const GCTAObservation& obs) const
{
    // Set constants
    const double offset_max = 10.0;
    const double u_max      = std::sin(offset_max * gammalib::deg2rad);

    // Return maximum value
    return u_max;
}


/***********************************************************************//**
 * @brief Returns integral over radial model (in steradians)
 *
 * Computes
 * \f[\Omega = 2 \pi \int_0^{\pi} \sin \theta f(\theta) d\theta\f]
 * where
 * \f[f(\theta) = (1 + (\theta/c_0)^{c_1})^{-c_2/c_1}\f]
 * \f$\theta\f$ is the offset angle (in degrees),
 * \f$c_0\f$ is the width of the profile (width),
 * \f$c_1\f$ is the width of the central plateau (core), and
 * \f$c_2\f$ is the size of the tail (tail).
 *
 * The integration is performed numerically, and the upper integration bound
 * \f$\pi\f$
 * is fixed to 10 degrees.
 ***************************************************************************/
double GCTAModelRadialProfile::omega(void) const
{
    // Set constants
    const double offset_max_rad = 10.0 * gammalib::deg2rad;

    // Allocate integrand
    GCTAModelRadialProfile::integrand integrand(this);

    // Allocate intergal
    GIntegral integral(&integrand);

    // Perform numerical integration
    double omega = integral.romberg(0.0, offset_max_rad) * gammalib::twopi;

    // Return integral
    return omega;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the Profile radial model information from an XML element in the
 * following format
 *
 *     <radialModel type="...">
 *       <parameter name="Width"  scale="1.0" value="1.5" min="0.1" max="10.0" free="1"/>
 *       <parameter name="Core"   scale="1.0" value="3.0" min="1.0" max="10.0" free="1"/>
 *       <parameter name="Tail"   scale="1.0" value="5.0" min="1.0" max="10.0" free="1"/>
 *     </radialModel>
 ***************************************************************************/
void GCTAModelRadialProfile::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameters
    gammalib::xml_check_parnum(G_READ, xml, 3);

    // Get parameters
    const GXmlElement* width = gammalib::xml_get_par(G_READ, xml, m_width.name());
    const GXmlElement* core  = gammalib::xml_get_par(G_READ, xml, m_core.name());
    const GXmlElement* tail  = gammalib::xml_get_par(G_READ, xml, m_tail.name());

    // Read parameters
    m_width.read(*width);
    m_core.read(*core);
    m_tail.read(*tail);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write the radial Profile model information into an XML element in the
 * following format
 *
 *     <radialModel type="...">
 *       <parameter name="Width"  scale="1.0" value="1.5" min="0.1" max="10.0" free="1"/>
 *       <parameter name="Core"   scale="1.0" value="3.0" min="1.0" max="10.0" free="1"/>
 *       <parameter name="Tail"   scale="1.0" value="5.0" min="1.0" max="10.0" free="1"/>
 *     </radialModel>
 ***************************************************************************/
void GCTAModelRadialProfile::write(GXmlElement& xml) const
{
    // Check model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get or create parameters
    GXmlElement* width = gammalib::xml_need_par(G_WRITE, xml, m_width.name());
    GXmlElement* core  = gammalib::xml_need_par(G_WRITE, xml, m_core.name());
    GXmlElement* tail  = gammalib::xml_need_par(G_WRITE, xml, m_tail.name());

    // Write parameters
    m_width.write(*width);
    m_core.write(*core);
    m_tail.write(*tail);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point source information.
 ***************************************************************************/
std::string GCTAModelRadialProfile::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelRadialProfile ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
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
void GCTAModelRadialProfile::init_members(void)
{
    // Initialise width parameter
    m_width.clear();
    m_width.name("Width");
    m_width.unit("deg");
    m_width.scale(1.0);
    m_width.value(1.5);          // default: 1.5 deg
    m_width.min(0.1);            // min:     0.1 deg
    m_width.free();
    m_width.gradient(0.0);
    m_width.has_grad(false);

    // Initialise core parameter
    m_core.clear();
    m_core.name("Core");
    m_core.scale(1.0);
    m_core.value(3.0);           // default: 3.0
    m_core.min(1.0);             // min:     1.0 (could even be larger)
    m_core.free();
    m_core.gradient(0.0);
    m_core.has_grad(false);

    // Initialise tail parameter
    m_tail.clear();
    m_tail.name("Tail");
    m_tail.scale(1.0);
    m_tail.value(5.0);           // default: 5.0
    m_tail.min(1.0);             // min:     1.0
    m_tail.fix();
    m_tail.gradient(0.0);
    m_tail.has_grad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_width);
    m_pars.push_back(&m_core);
    m_pars.push_back(&m_tail);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial Profileian model.
 ***************************************************************************/
void GCTAModelRadialProfile::copy_members(const GCTAModelRadialProfile& model)
{
    // Copy members
    m_width = model.m_width;
    m_core  = model.m_core;
    m_tail  = model.m_tail;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_width);
    m_pars.push_back(&m_core);
    m_pars.push_back(&m_tail);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRadialProfile::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
