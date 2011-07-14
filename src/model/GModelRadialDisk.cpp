/***************************************************************************
 *      GModelRadialDisk.cpp  -  Radial disk source model class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Christoph Deil                                   *
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
 * @file GModelRadialDisk.cpp
 * @brief Radial disk model class implementation
 * @author C. Deil
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelRadialDisk.hpp"
#include "GModelRadialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelRadialDisk     g_radial_disk_seed;
const GModelRadialRegistry g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                         "GModelRadialDisk::read(GXmlElement&)"
#define G_WRITE                       "GModelRadialDisk::write(GXmlElement&)"

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
 ***************************************************************************/
GModelRadialDisk::GModelRadialDisk(void) : GModelRadial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Disk constructor
 *
 * @param[in] dir Sky position of disk centre.
 * @param[in] radius Disk radius (degrees).
 ***************************************************************************/
GModelRadialDisk::GModelRadialDisk(const GSkyDir& dir,
                                   const double&  radius) : GModelRadial()
{
    // Initialise members
    init_members();

    // Assign parameters
    this->dir(dir);
    this->radius(radius);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of radial disk model by extracting information from an
 * XML element. See GModelRadialDisk::read() for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelRadialDisk::GModelRadialDisk(const GXmlElement& xml) : GModelRadial()
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
 * @param[in] model Radial disk model.
 ***************************************************************************/
GModelRadialDisk::GModelRadialDisk(const GModelRadialDisk& model)
                                                        : GModelRadial(model)
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
GModelRadialDisk::~GModelRadialDisk(void)
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
 * @param[in] model Radial disk model.
 ***************************************************************************/
GModelRadialDisk& GModelRadialDisk::operator= (const GModelRadialDisk& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelRadial::operator=(model);

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
void GModelRadialDisk::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelRadial::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelRadial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelRadialDisk* GModelRadialDisk::clone(void) const
{
    return new GModelRadialDisk(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from shell centre (radians).
 *
 * Evaluates the spatial part for a disk source model. The disk source
 * model is a radial function \f$f(\theta)\f$, where \f$\theta\f$ is the
 * angular separation between shell centre and the actual location:
 * \f[
 * f(\theta) = \left \{
 *  \begin{array}{l l}
 *     \displaystyle
 *     2 \pi ( 1 - \cos(\theta))
 *     & \mbox{if $\theta \le $ radius} \\
 *     \\
 *    \displaystyle
 *    0 & \mbox{if $\theta > $ radius}
 *  \end{array}
 *  \right .
 * \f]
 ***************************************************************************/
double GModelRadialDisk::eval(const double& theta) const
{
    // Update precomputation cache
    update();

    // Set value
    double value = (theta <= m_radius_rad) ? m_norm : 0.0;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnan(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelRadialDisk::eval";
        std::cout << "(theta=" << theta << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", m_radius_rad=" << m_radius_rad;
        std::cout << ", m_norm=" << m_norm;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients (in units of sr^-1)
 *
 * @param[in] theta Angular distance from shell centre (radians).
 *
 * Evaluates the function value. No gradient computation is implemented as
 * radial models will be convolved with the instrument response and thus
 * require the numerical computation of the derivatives.
 ***************************************************************************/
double GModelRadialDisk::eval_gradients(const double& theta) const
{
    // Return value
    return (eval(theta));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary sky position from the 2D disk distribution.
 ***************************************************************************/
GSkyDir GModelRadialDisk::mc(GRan& ran) const
{
    // Simulate offset from photon arrival direction
    double cosrad = std::cos(radius()*deg2rad);
    double theta  = std::acos(1.0 - ran.uniform() * (1.0 - cosrad)) * rad2deg;
    double phi    = 360.0 * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 ***************************************************************************/
double GModelRadialDisk::theta_max(void) const
{
    // Return value
    return (radius()*deg2rad);
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
 * Read the point source information from an XML element. The XML element
 * is required to have 3 parameters.
 * The position is named either "RA" and "DEC" or "GLON" and "GLAT", the
 * disk radius is named "Radius".
 *
 * @todo Implement a test of the radius and radius boundary. The radius
 *       and radius minimum should be >0.
 ***************************************************************************/
void GModelRadialDisk::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Disk model requires exactly 3 parameters.");

    // Read disk location
    GModelRadial::read(xml);

    // Extract model parameters
    int  npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle Radius
        if (par->attribute("name") == "Radius") {
            
            // Read parameter
            m_radius.read(*par);
            
            //TODO: Check parameter
            
            // Increment parameter counter
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Radius\" parameter.");

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
 * Write the Disk source information into an XML element. The XML element
 * will have 3 parameter leafs named "RA", "DEC" and "Radius"
 ***************************************************************************/
void GModelRadialDisk::write(GXmlElement& xml) const
{
    // Write disk location
    GModelRadial::write(xml);

    // If XML element has 2 nodes (which should be the location nodes)
    // then append 1 parameter node
    if (xml.elements() == 2) {
        xml.append(new GXmlElement("parameter name=\"Radius\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Disk model requires exactly 3 parameters.");

    // Set or update model parameter attributes
    int npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle Radius
        if (par->attribute("name") == "Radius") {
            m_radius.write(*par);
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Radius\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 ***************************************************************************/
std::string GModelRadialDisk::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelRadialDisk ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

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
void GModelRadialDisk::init_members(void)
{
    // Initialise Radius
    m_radius.clear();
    m_radius.name("Radius");
    m_radius.unit("deg");
    m_radius.value(2.778e-4); // 1 arcsec
    m_radius.min(2.778e-4);   // 1 arcsec
    m_radius.free();
    m_radius.scale(1.0);
    m_radius.gradient(0.0);
    m_radius.hasgrad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_radius);

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius and width shell is not meaningful
    m_last_radius = 0.0;
    m_radius_rad  = 0.0;
    m_norm        = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial disk model.
 *
 * We do not have to push back the members on the parameter stack as this
 * should have been done by init_members() that was called before. Otherwise
 * we would have the radius twice on the stack.
 ***************************************************************************/
void GModelRadialDisk::copy_members(const GModelRadialDisk& model)
{
    // Copy members
    m_radius = model.m_radius;

    // Copy precomputation cache
    m_last_radius = model.m_last_radius;
    m_radius_rad  = model.m_radius_rad;
    m_norm        = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelRadialDisk::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the normalization
 * \f[{\tt m\_norm} = \frac{1}{2 \pi (1 - \cos r)}\f]
 ***************************************************************************/
void GModelRadialDisk::update() const
{
    // Update if radius has changed
    if (m_last_radius != radius()) {

        // Store last values
        m_last_radius = radius();

        // Compute disk radius in radians
        m_radius_rad = radius() * deg2rad;

        // Perform precomputations
        double denom = twopi * (1 - std::cos(m_radius_rad));
        m_norm       = (denom > 0.0) ? 1.0 / denom : 0.0;

    } // endif: update required

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
