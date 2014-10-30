/***************************************************************************
 *   GModelSpatialEllipticalDisk.cpp - Elliptical disk source model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Michael Mayer                               *
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
 * @file GModelSpatialEllipticalDisk.cpp
 * @brief Elliptical disk model class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialEllipticalDisk.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialEllipticalDisk g_elliptical_disk_seed;
const GModelSpatialRegistry       g_elliptical_disk_registry(&g_elliptical_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ              "GModelSpatialEllipticalDisk::read(GXmlElement&)"
#define G_WRITE            "GModelSpatialEllipticalDisk::write(GXmlElement&)"

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
GModelSpatialEllipticalDisk::GModelSpatialEllipticalDisk(void) :
                             GModelSpatialElliptical()
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
 * @param[in] semimajor Semi-major axis (degrees).
 * @param[in] semiminor Semi-minor axis (degrees).
 * @param[in] posangle Position angle of semi-major axis (degrees).
 *
 * Construct elliptical disk model from sky position of the ellipse centre
 * (@p dir), the @p semimajor and @p semiminor axes, and the position
 * angle (@p posangle).
 ***************************************************************************/
GModelSpatialEllipticalDisk::GModelSpatialEllipticalDisk(const GSkyDir& dir,
                                                         const double&  semimajor,
                                                         const double&  semiminor,
                                                         const double&  posangle) :
                             GModelSpatialElliptical()
{
    // Initialise members
    init_members();

    // Assign parameters
    this->dir(dir);
    this->semiminor(semiminor);
    this->semimajor(semimajor);
    this->posangle(posangle);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs elliptical disk model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialEllipticalDisk::GModelSpatialEllipticalDisk(const GXmlElement& xml) :
                             GModelSpatialElliptical()
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
 * @param[in] model Elliptical disk model.
 ***************************************************************************/
GModelSpatialEllipticalDisk::GModelSpatialEllipticalDisk(const GModelSpatialEllipticalDisk& model) :
                             GModelSpatialElliptical(model)
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
GModelSpatialEllipticalDisk::~GModelSpatialEllipticalDisk(void)
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
 * @param[in] model Elliptical disk model.
 * @return Elliptical disk model.
 ***************************************************************************/
GModelSpatialEllipticalDisk& GModelSpatialEllipticalDisk::operator=(const GModelSpatialEllipticalDisk& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialElliptical::operator=(model);

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
 * @brief Clear elliptical disk model
 ***************************************************************************/
void GModelSpatialEllipticalDisk::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatialElliptical::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialElliptical::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone elliptical disk model
 *
 * @return Pointer to deep copy of elliptical disk model.
 ***************************************************************************/
GModelSpatialEllipticalDisk* GModelSpatialEllipticalDisk::clone(void) const
{
    // Clone elliptical disk model
    return new GModelSpatialEllipticalDisk(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from disk centre (radians).
 * @param[in] posangle Position angle (counterclockwise from North) (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the spatial component for an elliptical disk source model. The
 * disk source model is an elliptical function
 * \f$S_{\rm p}(\theta, \phi | E, t)\f$, where
 * \f$\theta\f$ is the angular separation between elliptical disk centre and
 * the actual location and \f$\phi\f$ the position angle with respect to the
 * ellipse centre, counted counterclockwise from North.
 *
 * The function \f$f(\theta, \phi)\f$ is given by
 *
 * \f[
 * S_{\rm p}(\theta, \phi | E, t) = \left \{
 *  \begin{array}{l l}
 *     {\tt m\_norm}
 *     & \mbox{if} \, \, \theta \le \theta_0 \\
 *     \\
 *     0 & \mbox{else}
 *  \end{array}
 *  \right .
 * \f]
 *
 * where \f$\theta_0\f$ is the effective radius of the ellipse on the sphere
 * given by
 *
 * \f[\theta_0\ =
 *    \frac{ab}{\sqrt{b^2 \cos^2(\phi-\phi_0) + a^2 \sin^2(\phi-\phi_0)}}\f]
 *
 * and
 * \f$a\f$ is the semi-major axis of the ellipse,
 * \f$b\f$ is the semi-minor axis, and
 * \f$\phi_0\f$ is the position angle of the ellipse, counted
 * counterclockwise from North.
 *
 * The normalisation constant \f${\tt m\_norm}\f$ which is the inverse of the
 * solid angle subtended by an ellipse is given by
 *
 * @todo Quote formula for ellipse solid angle
 *
 * (see the update() method).
 ***************************************************************************/
double GModelSpatialEllipticalDisk::eval(const double&  theta,
                                         const double&  posangle,
                                         const GEnergy& energy,
                                         const GTime&   time) const
{
    // Initialise value
    double value = 0.0;

    // Continue only if we're inside circle enclosing the ellipse
    if (theta <= theta_max()) {

        // Update precomputation cache
        update();

        // Perform computations
        double diff_angle = posangle - m_posangle.value() * gammalib::deg2rad;
        double cosinus    = std::cos(diff_angle);
        double sinus      = std::sin(diff_angle);
        double arg1       = m_semiminor_rad * cosinus;
        double arg2       = m_semimajor_rad * sinus;
        double r_ell      = m_semiminor_rad * m_semimajor_rad /
                            std::sqrt(arg1*arg1 + arg2*arg2);

        // Set value
        value = (theta <= r_ell) ? m_norm : 0.0;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GModelSpatialEllipticalDisk::eval";
            std::cout << "(theta=" << theta << "): NaN/Inf encountered";
            std::cout << "(posangle=" << posangle << "): NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", R_ellipse=" << r_ell;
            std::cout << ", diff_angle=" << diff_angle;
            std::cout << ", m_norm=" << m_norm;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: position was inside enclosing circle

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients (in units of sr^-1)
 *
 * @param[in] theta Angular distance from disk centre (radians).
 * @param[in] posangle Position angle (counterclockwise from North) (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the function value. No gradient computation is implemented as
 * Elliptical models will be convolved with the instrument response and thus
 * require the numerical computation of the derivatives.
 *
 * See the eval() method for more information.
 ***************************************************************************/
double GModelSpatialEllipticalDisk::eval_gradients(const double&  theta,
                                                   const double&  posangle,
                                                   const GEnergy& energy,
                                                   const GTime&   time) const
{
    // Return value
    return (eval(theta, posangle, energy, time));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky position from the 2D disk distribution.
 *
 * @todo Test function
 ***************************************************************************/
GSkyDir GModelSpatialEllipticalDisk::mc(const GEnergy& energy,
                                        const GTime&   time,
                                        GRan&          ran) const
{
    // Update precomputation cache
	update();

	// Initialise photon
    GPhoton photon;
    photon.energy(energy);
    photon.time(time);

	// Draw randomly from the radial disk
	// and reject the value if its outside the ellipse
	do {

        // Simulate offset from photon arrival direction
        double cosrad = std::cos(semimajor() * gammalib::deg2rad);
        double theta  = std::acos(1.0 - ran.uniform() * (1.0 - cosrad)) *
                        gammalib::rad2deg;
        double phi    = 360.0 * ran.uniform();

        // Rotate sky direction by offset
        GSkyDir sky_dir = dir();
        sky_dir.rotate_deg(phi, theta);

        // Set photon sky direction
        photon.dir(sky_dir);

	} while(GModelSpatialElliptical::eval(photon) <= 0.0);

	// Return photon direction
	return (photon.dir());

}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Returns maximum model radius.
 ***************************************************************************/
double GModelSpatialEllipticalDisk::theta_max(void) const
{
    // Set maximum model radius
    double theta_max = (semimajor() > semiminor())
                       ? semimajor() * gammalib::deg2rad
                       : semiminor() * gammalib::deg2rad;

    // Return value
    return theta_max;
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
 * Reads the elliptical disk model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="RA"          scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="GLON"        scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"        scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement a test of the ellipse boundary. The axes
 *       and axes minimum should be >0.
 ***************************************************************************/
void GModelSpatialEllipticalDisk::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 5 parameters
    if (xml.elements() != 5 || npars != 5) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Elliptical disk model requires exactly 5 parameters.");
    }

    // Read disk location
    GModelSpatialElliptical::read(xml);

    // Extract model parameters
    int  npar[2] = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle semiminor radius
        if (par->attribute("name") == "MinorRadius") {
            
            // Read parameter
            m_semiminor.read(*par);
            
            // Increment parameter counter
            npar[0]++;
        }

        // Handle semimajor radius
        else if (par->attribute("name") == "MajorRadius") {

        	// Read parameter
        	m_semimajor.read(*par);

        	// Increment parameter counter
        	npar[1]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Elliptical disk model requires \"MinorRadius\" and"
              " \"MajorRadius\" parameters.");
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
 * Write the elliptical disk model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="RA"          scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
 *       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
 *       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialEllipticalDisk::write(GXmlElement& xml) const
{
    // Write disk location
    GModelSpatialElliptical::write(xml);

    // If XML element has 3 nodes (which should be the location and PA nodes)
    // then append 2 parameter nodes
    if (xml.elements() == 3) {
        xml.append(GXmlElement("parameter name=\"MinorRadius\""));
        xml.append(GXmlElement("parameter name=\"MajorRadius\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 5 parameters
    if (xml.elements() != 5 || npars != 5) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Elliptical Disk model requires exactly 5 parameters.");
    }

    // Set or update model parameter attributes
    int npar[2] = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle semiminor radius
        if (par->attribute("name") == "MinorRadius") {

        	// Write parameter
            m_semiminor.write(*par);

            // Increment parameter counter
            npar[0]++;
        }

        // Handle semimajor radius
        else if (par->attribute("name") == "MajorRadius") {

        	// Write parameter
            m_semimajor.write(*par);

            // Increment parameter counter
            npar[1]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Elliptical disk model requires \"MinorRadius\" and"
              " \"MajorRadius\" parameters.");
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
std::string GModelSpatialEllipticalDisk::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialEllipticalDisk ===");

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
void GModelSpatialEllipticalDisk::init_members(void)
{
    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius is not meaningful
    m_last_semiminor = 0.0;
    m_last_semimajor = 0.0;
    m_semiminor_rad  = 0.0;
    m_semimajor_rad  = 0.0;
    m_norm           = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Elliptical disk model.
 *
 * We do not have to push back the members on the parameter stack as this
 * should have been done by init_members() that was called before. Otherwise
 * we would have the radius twice on the stack.
 ***************************************************************************/
void GModelSpatialEllipticalDisk::copy_members(const GModelSpatialEllipticalDisk& model)
{
    // Copy precomputation cache
    m_last_semiminor = model.m_last_semiminor;
    m_last_semimajor = model.m_last_semimajor;
    m_semiminor_rad  = model.m_semiminor_rad;
    m_semimajor_rad  = model.m_semimajor_rad;
    m_norm           = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialEllipticalDisk::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the normalization
 * \f[
 *    {\tt m\_norm} = \frac{1}{2 \pi (1 - \cos a) (1 - \cos b)}
 * \f]
 *
 * @todo check this formula
 ***************************************************************************/
void GModelSpatialEllipticalDisk::update() const
{
    // Update if one axis has changed
    if (m_last_semiminor != semiminor() || m_last_semimajor != semimajor()) {

        // Store last values
        m_last_semiminor = semiminor();
        m_last_semimajor = semimajor();

        // Compute axes in radians
        m_semiminor_rad = semiminor() * gammalib::deg2rad;
        m_semimajor_rad = semimajor() * gammalib::deg2rad;

        // Perform precomputations
        double denom = gammalib::twopi *
                       std::sqrt((1.0 - std::cos(m_semiminor_rad)) *
                                 (1.0 - std::cos(m_semimajor_rad)));
        m_norm       = (denom > 0.0) ? 1.0 / denom : 0.0;

    } // endif: update required

    // Return
    return;
}
