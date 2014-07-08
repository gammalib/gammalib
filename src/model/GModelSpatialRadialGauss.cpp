/***************************************************************************
 *    GModelSpatialRadialGauss.cpp - Radial Gaussian source model class    *
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
 * @file GModelSpatialRadialGauss.cpp
 * @brief Radial Gaussian model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialRadialGauss.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialGauss g_radial_gauss_seed;
const GModelSpatialRegistry    g_radial_gauss_registry(&g_radial_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GModelSpatialRadialGauss::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialRadialGauss::write(GXmlElement&)"

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
GModelSpatialRadialGauss::GModelSpatialRadialGauss(void) : GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] dir Sky position of Gaussian.
 * @param[in] sigma Width of Gaussian (in degrees).
 *
 * Constructs a Gaussian spatial model using a sky direction (@p dir) and
 * a Gaussian width parameter @p sigma in degrees.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const GSkyDir& dir,
                                                   const double&  sigma) :
                          GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Assign direction and sigma
    this->dir(dir);
    this->sigma(sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a Gaussian spatial model by extracting information from an XML
 * element. See the method read() for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const GXmlElement& xml) :
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
 * @param[in] model Radial Gaussian model.
 ***************************************************************************/
GModelSpatialRadialGauss::GModelSpatialRadialGauss(const GModelSpatialRadialGauss& model) : 
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
GModelSpatialRadialGauss::~GModelSpatialRadialGauss(void)
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
 * @param[in] model Radial Gaussian model.
 * @return Radial Gaussian model.
 ***************************************************************************/
GModelSpatialRadialGauss& GModelSpatialRadialGauss::operator=(const GModelSpatialRadialGauss& model)
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
 * @brief Clear radial Gauss model
 ***************************************************************************/
void GModelSpatialRadialGauss::clear(void)
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
 * @brief Clone radial Gauss model
 ***************************************************************************/
GModelSpatialRadialGauss* GModelSpatialRadialGauss::clone(void) const
{
    // Clone radial Gauss model
    return new GModelSpatialRadialGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] theta Angular distance from Gaussian centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * Evaluates the spatial component for a Gaussian source model using
 *
 * \f[
 *    S_{\rm p}(\vec{p} | E, t) =
 *    \frac{1}{2 \pi \sigma^2} \exp 
 *    \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)
 * \f]
 *
 * where
 * - \f$\theta\f$ is the angular separation from the source direction, and
 * - \f$\sigma\f$ is the Gaussian width.
 *
 * @todo The Gaussian function is only correct in the small angle
 *       approximation.
 ***************************************************************************/
double GModelSpatialRadialGauss::eval(const double&  theta,
                                      const GEnergy& energy,
                                      const GTime&   time) const
{
    // Compute value
    double sigma_rad = sigma() * gammalib::deg2rad;
    double sigma2    = sigma_rad * sigma_rad;
    double theta2    = theta   * theta;
    double value     = std::exp(-0.5 * theta2 / sigma2) /
                       (gammalib::twopi * sigma2) / 0.8;
    if ( theta > theta_max() ) {
    	value = 0.0;
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpatialRadialGauss::eval";
        std::cout << "(theta=" << theta << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", sigma_rad=" << sigma_rad;
        std::cout << ", sigma2=" << sigma2;
        std::cout << ", theta2=" << theta2;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] theta Angular distance from Gaussian centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @return Model value.
 *
 * This method simply calls the eval() method as no analytical gradients
 * will be computed. See the eval() method for more details.
 ***************************************************************************/
double GModelSpatialRadialGauss::eval_gradients(const double&  theta,
                                                const GEnergy& energy,
                                                const GTime&   time) const
{
    // Return value
    return (eval(theta, energy, time));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction from the 2D Gaussian distribution as
 * function of the photon @p energy and arrival @p time.
 *
 * @todo This method is only valid in the small angle approximation.
 ***************************************************************************/
GSkyDir GModelSpatialRadialGauss::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Simulate offset from photon arrival direction
    double theta = sigma() * ran.chisq2();
    double phi   = 360.0   * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * Returns \f$5 \sigma\f$ as approximate edge of the Gaussian. This limit
 * is of course arbitrary, but allows to limit the integration region for
 * response computation.
 ***************************************************************************/
double GModelSpatialRadialGauss::theta_max(void) const
{
    // Return value
    return (sigma() * gammalib::deg2rad * 5.0);
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
 * Reads the radial Gauss model information from an XML element. The XML
 * element shall have either the format 
 *
 *     <spatialModel type="GaussFunction">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="GaussFunction">
 *       <parameter name="GLON"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"  scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * @todo Implement a test of the sigma and sigma boundary. The sigma
 *       and sigma minimum should be >0.
 ***************************************************************************/
void GModelSpatialRadialGauss::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Gaussian source model requires exactly 3 parameters.");
    }

    // Read Gaussian location
    GModelSpatialRadial::read(xml);

    // Extract model parameters
    int  npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle Gaussian width
        if (par->attribute("name") == "Sigma") {
            
            // Read parameter
            m_sigma.read(*par);
            
            //TODO: Check parameter
            
            // Increment parameter counter
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Sigma\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type 'GaussFunction'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the radial disk model information into an XML element. The XML
 * element will have the format 
 *
 *     <spatialModel type="DiskFunction">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadialGauss::write(GXmlElement& xml) const
{
    // Write Gaussian location
    GModelSpatialRadial::write(xml);

    // If XML element has 2 nodes (which should be the location nodes)
    // then append 1 parameter node
    if (xml.elements() == 2) {
        xml.append(GXmlElement("parameter name=\"Sigma\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Point source model requires exactly 3 parameters.");
    }

    // Set or update model parameter attributes
    int npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle Sigma
        if (par->attribute("name") == "Sigma") {
            m_sigma.write(*par);
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Sigma\" parameter.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Gaussian source information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialGauss ===");

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
 *
 * Note that this model implements no gradients, as spatial models will
 * always require numerical gradient computations. The minimum Gaussian
 * width is set to 1 arcsec.
 ***************************************************************************/
void GModelSpatialRadialGauss::init_members(void)
{

    // Initialise Gaussian sigma
    m_sigma.clear();
    m_sigma.name("Sigma");
    m_sigma.unit("deg");
    m_sigma.value(2.778e-4); // 1 arcsec
    m_sigma.min(2.778e-4);   // 1 arcsec
    m_sigma.free();
    m_sigma.scale(1.0);
    m_sigma.gradient(0.0);
    m_sigma.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Gaussian spatial model.
 *
 * We do not have to push back the members on the parameter stack as this
 * should have been done by init_members() that was called before. Otherwise
 * we would have sigma twice on the stack.
 ***************************************************************************/
void GModelSpatialRadialGauss::copy_members(const GModelSpatialRadialGauss& model)
{
    // Copy members
    m_sigma = model.m_sigma;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialGauss::free_members(void)
{
    // Return
    return;
}
