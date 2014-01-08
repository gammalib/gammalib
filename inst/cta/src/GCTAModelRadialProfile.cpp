/***************************************************************************
 *       GCTAModelRadialProfile.cpp - Radial Profile CTA model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
#include "GIntegral.hpp"
#include "GCTAModelRadialProfile.hpp"
#include "GCTAModelRadialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRadialProfile  g_cta_radial_profile_seed;
const GCTAModelRadialRegistry g_cta_radial_profile_registry(&g_cta_radial_profile_seed);

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
 * Note that this method implements a function which is unity for
 * \f$\theta=0\f$.
 ***************************************************************************/
double GCTAModelRadialProfile::eval(const double& offset) const
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
 * @brief Evaluate function and gradients
 *
 * @param[in] offset Offset angle [degrees].
 *
 * This method simply calls GCTAModelRadialProfile::eval() as no analytical
 * gradients will be computed. See GCTAModelRadialProfile::eval() for details
 * about the implemented method.
 ***************************************************************************/
double GCTAModelRadialProfile::eval_gradients(const double& offset) const
{
    // Return value
    return (eval(offset));
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] dir Pointing direction.
 * @param[in] ran Random number generator.
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
GCTAInstDir GCTAModelRadialProfile::mc(const GCTAInstDir& dir, GRan& ran) const
{
    // Set constants
    const double offset_max = 10.0;
    const double u_max      = sin(offset_max * gammalib::deg2rad);

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
        value  = sin(offset * gammalib::deg2rad) * eval(offset);
        
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

    // Rotate pointing direction by offset and azimuth angle
    GCTAInstDir mc_dir = dir;
    mc_dir.dir().rotate_deg(phi, offset);

    // Return MC direction
    return mc_dir;
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
    double omega = integral.romb(0.0, offset_max_rad) * gammalib::twopi;

    // Return integral
    return omega;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parvalue
 *            Non-positive parameter value found.
 * @exception GException::model_invalid_parlimit
 *            Missing or non-positive minimum parameter boundary.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the Profile radial model information from an XML element. The XML
 * element is required to have 3 parameter named "Width", "Core", and "Tail". 
 ***************************************************************************/
void GCTAModelRadialProfile::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameter
    if (xml.elements() != 3 || xml.elements("parameter") != 3) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Radial Profile model requires exactly 3 parameters.");
    }

    // Extract model parameters
    int  npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle Width
        if (par->attribute("name") == "Width") {
            m_width.read(*par);
            if (m_width.value() <= 0.0) {
                throw GException::model_invalid_parvalue(G_READ, xml,
                      "\"Width\" parameter is required to be positive.");
            }
            if (!m_width.has_min() || m_width.min() <= 0.0) {
                throw GException::model_invalid_parlimit(G_READ, xml,
                      "\"Width\" parameter requires positive minimum boundary.");
            }
            npar[0]++;
        }

        // Handle Core
        if (par->attribute("name") == "Core") {
            m_core.read(*par);
            if (m_core.value() <= 0.0) {
                throw GException::model_invalid_parvalue(G_READ, xml,
                      "\"Core\" parameter is required to be positive.");
            }
            if (!m_core.has_min() || m_core.min() <= 0.0) {
                throw GException::model_invalid_parlimit(G_READ, xml,
                      "\"Core\" parameter requires positive minimum boundary.");
            }
            npar[1]++;
        }

        // Handle Tail
        if (par->attribute("name") == "Tail") {
            m_tail.read(*par);
            if (m_tail.value() <= 0.0) {
                throw GException::model_invalid_parvalue(G_READ, xml,
                      "\"Core\" parameter is required to be positive.");
            }
            if (!m_tail.has_min() || m_tail.min() <= 0.0) {
                throw GException::model_invalid_parlimit(G_READ, xml,
                      "\"Core\" parameter requires positive minimum boundary.");
            }
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Width\", \"Core\" and \"Tail\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type 'ProfileFunction'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the radial Profile model information into an XML element. The XML
 * element will have 3 parameter leafs named "Width", "Core" and "Tail".
 ***************************************************************************/
void GCTAModelRadialProfile::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Radial Profile model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append 3 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Width\""));
        xml.append(GXmlElement("parameter name=\"Core\""));
        xml.append(GXmlElement("parameter name=\"Tail\""));
    }

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Radial Profile model requires exactly 3 parameters.");
    }

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Width") {
            m_width.write(*par);
            npar[0]++;
        }

        // Handle index
        else if (par->attribute("name") == "Core") {
            m_core.write(*par);
            npar[1]++;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Tail") {
            m_tail.write(*par);
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Width\", \"Core\" and \"Tail\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
