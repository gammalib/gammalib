/***************************************************************************
 *      GModelSpatialRadialDisk.cpp - Radial disk source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Christoph Deil                              *
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
 * @file GModelSpatialRadialDisk.cpp
 * @brief Radial disk model class implementation
 * @author Christoph Deil
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialDisk g_radial_disk_seed;
const GModelSpatialRegistry   g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                  "GModelSpatialRadialDisk::read(GXmlElement&)"
#define G_WRITE                "GModelSpatialRadialDisk::write(GXmlElement&)"

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
GModelSpatialRadialDisk::GModelSpatialRadialDisk(void) : GModelSpatialRadial()
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
 *
 * Constructs radial disk model from the sky position of the disk centre
 * (@p dir) and the disk @p radius in degrees.
 ***************************************************************************/
GModelSpatialRadialDisk::GModelSpatialRadialDisk(const GSkyDir& dir,
                                                 const double&  radius) :
                         GModelSpatialRadial()
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
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialDisk::GModelSpatialRadialDisk(const GXmlElement& xml) :
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
 * @param[in] model Radial disk model.
 ***************************************************************************/
GModelSpatialRadialDisk::GModelSpatialRadialDisk(const GModelSpatialRadialDisk& model) :
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
GModelSpatialRadialDisk::~GModelSpatialRadialDisk(void)
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
 * @return Radial disk model.
 ***************************************************************************/
GModelSpatialRadialDisk& GModelSpatialRadialDisk::operator=(const GModelSpatialRadialDisk& model)
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
 * @brief Clear radial disk model
 ***************************************************************************/
void GModelSpatialRadialDisk::clear(void)
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
 * @brief Clone radial disk model
 *
 * @return Pointer to deep copy of radial disk model.
 ***************************************************************************/
GModelSpatialRadialDisk* GModelSpatialRadialDisk::clone(void) const
{
    // Clone radial disk model
    return new GModelSpatialRadialDisk(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from disk centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the spatial component for a disk source model using
 *
 * \f[
 *    S_{\rm p}(\vec{p} | E, t) = \left \{
 *    \begin{array}{l l}
 *       \displaystyle
 *       {\tt m\_norm}
 *       & \mbox{if $\theta \le $ radius} \\
 *       \\
 *      \displaystyle
 *      0 & \mbox{if $\theta > $ radius}
 *    \end{array}
 *    \right .
 * \f]
 *
 * where
 * - \f$\theta\f$ is the angular separation between disk centre and the
 *   sky direction of interest, and
 * - \f${\tt m\_norm} = \frac{1}{2 \pi (1 - \cos r)} \f$ is a normalization
 *   constant (see the update() method for details).
 ***************************************************************************/
double GModelSpatialRadialDisk::eval(const double&  theta,
                                     const GEnergy& energy,
                                     const GTime&   time) const
{
    // Update precomputation cache
    update();

    // Set value
    double value = (theta <= m_radius_rad) ? m_norm : 0.0;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpatialRadialDisk::eval";
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
 * @param[in] theta Angular distance from disk centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the function value. No gradient computation is implemented as
 * radial models will be convolved with the instrument response and thus
 * require the numerical computation of the derivatives.
 *
 * See the eval() method for more information.
 ***************************************************************************/
double GModelSpatialRadialDisk::eval_gradients(const double&  theta,
                                               const GEnergy& energy,
                                               const GTime&   time) const
{
    // Return value
    return (eval(theta, energy, time));
}


/***********************************************************************//**
 * @brief Return MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction from the 2D disk distribution as function
 * of the photon @p energy and arrival @p time.
 ***************************************************************************/
GSkyDir GModelSpatialRadialDisk::mc(const GEnergy& energy,
                                    const GTime&   time,
                                    GRan&          ran) const
{
    // Simulate offset from photon arrival direction
    double cosrad = std::cos(radius() * gammalib::deg2rad);
    double theta  = std::acos(1.0 - ran.uniform() * (1.0 - cosrad)) *
                    gammalib::rad2deg;
    double phi    = 360.0 * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialDisk::theta_max(void) const
{
    // Return value
    return (radius() * gammalib::deg2rad);
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
 * Reads the radial disk model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="GLON"   scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 * 
 * @todo Implement a test of the radius and radius boundary. The radius
 *       and radius minimum should be >0.
 ***************************************************************************/
void GModelSpatialRadialDisk::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Disk model requires exactly 3 parameters.");
    }

    // Read disk location
    GModelSpatialRadial::read(xml);

    // Extract model parameters
    int  npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

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
    if (npar[0] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Radius\" parameter.");
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
 * Writes the radial disk model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Radius" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadialDisk::write(GXmlElement& xml) const
{
    // Write disk location
    GModelSpatialRadial::write(xml);

    // If XML element has 2 nodes (which should be the location nodes)
    // then append 1 parameter node
    if (xml.elements() == 2) {
        xml.append(GXmlElement("parameter name=\"Radius\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Disk model requires exactly 3 parameters.");
    }

    // Set or update model parameter attributes
    int npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle Radius
        if (par->attribute("name") == "Radius") {
            m_radius.write(*par);
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Radius\" parameter.");
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
std::string GModelSpatialRadialDisk::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialDisk ===");

        // Append parameters
        result.append("\n"+parformat("Number of parameters")+str(size()));
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
void GModelSpatialRadialDisk::init_members(void)
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
    // uninitialised as a zero radius is not meaningful
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
void GModelSpatialRadialDisk::copy_members(const GModelSpatialRadialDisk& model)
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
void GModelSpatialRadialDisk::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the normalization
 * \f[
 *    {\tt m\_norm} = \frac{1}{2 \pi (1 - \cos r)}
 * \f]
 *
 * Note that this is the correct normalization on the sphere for any
 * disk radius r. For small r it is very similar to the cartesian
 * approximation you might have expected:
 * \f[{\tt m\_norm} = \frac{1}{\pi r ^ 2}\f]
 ***************************************************************************/
void GModelSpatialRadialDisk::update() const
{
    // Update if radius has changed
    if (m_last_radius != radius()) {

        // Store last values
        m_last_radius = radius();

        // Compute disk radius in radians
        m_radius_rad = radius() * gammalib::deg2rad;

        // Perform precomputations
        double denom = gammalib::twopi * (1 - std::cos(m_radius_rad));
        m_norm       = (denom > 0.0) ? 1.0 / denom : 0.0;

    } // endif: update required

    // Return
    return;
}
