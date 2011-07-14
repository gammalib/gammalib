/***************************************************************************
 *      GModelRadialGauss.cpp  -  Radial Gaussian source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GModelRadialGauss.cpp
 * @brief Radial Gaussian model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelRadialGauss.hpp"
#include "GModelRadialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelRadialGauss    g_radial_gauss_seed;
const GModelRadialRegistry g_radial_gauss_registry(&g_radial_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                        "GModelRadialGauss::read(GXmlElement&)"
#define G_WRITE                      "GModelRadialGauss::write(GXmlElement&)"

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
GModelRadialGauss::GModelRadialGauss(void) : GModelRadial()
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
 * Creates instance of a Gaussian spatial model using a sky direction and
 * a Gaussian width parameter \f$\sigma\f$ (in degrees).
 ***************************************************************************/
GModelRadialGauss::GModelRadialGauss(const GSkyDir& dir,
                                     const double&  sigma) : GModelRadial()
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
 * Creates instance of a Gaussian spatial model by extracting information
 * from an XML element. See GModelRadialGauss::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelRadialGauss::GModelRadialGauss(const GXmlElement& xml) : GModelRadial()
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
GModelRadialGauss::GModelRadialGauss(const GModelRadialGauss& model)
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
GModelRadialGauss::~GModelRadialGauss(void)
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
 ***************************************************************************/
GModelRadialGauss& GModelRadialGauss::operator=(const GModelRadialGauss& model)
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
void GModelRadialGauss::clear(void)
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
GModelRadialGauss* GModelRadialGauss::clone(void) const
{
    return new GModelRadialGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] theta Angular distance from Gaussian centre (radians).
 *
 * Evaluates the spatial part for a Gaussian source model. The Gaussian
 * source model is defined as
 * \f[f(\theta)=\frac{1}{2 \pi \sigma^2} \exp 
 *    \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)\f]
 * where
 * \f$\theta\f$ is the angular separation from the source direction, and
 * \f$\sigma\f$ is the Gaussian width.
 *
 * @todo The Gaussian function is only correct in the small angle
 *       approximation.
 ***************************************************************************/
double GModelRadialGauss::eval(const double& theta) const
{
    // Compute value
    double sigma_rad = sigma() * deg2rad;
    double sigma2    = sigma_rad * sigma_rad;
    double theta2    = theta   * theta;
    double value     = exp(-0.5 * theta2 / sigma2) / (twopi * sigma2);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelRadialGauss::eval";
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
 *
 * This method simply calls GModelRadialShell::eval() as no analytical
 * gradients will be computed. See GModelRadialShell::eval() for details
 * about the implemented method.
 ***************************************************************************/
double GModelRadialGauss::eval_gradients(const double& theta) const
{
    // Return value
    return (eval(theta));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary sky position from the 2D Gaussian distribution.
 *
 * @todo This method is only valid in the small angle approximation.
 ***************************************************************************/
GSkyDir GModelRadialGauss::mc(GRan& ran) const
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
double GModelRadialGauss::theta_max(void) const
{
    // Return value
    return (sigma()*deg2rad*5.0);
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
 * Read the Gaussian source information from an XML element. The XML element
 * is required to have 3 parameters. 
 * The position is named either "RA" and "DEC" or "GLON" and "GLAT", the
 * Gaussian width is named "Sigma".
 *
 * @todo Implement a test of the sigma and sigma boundary. The sigma
 *       and sigma minimum should be >0.
 ***************************************************************************/
void GModelRadialGauss::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Gaussian source model requires exactly 3 parameters.");

    // Read Gaussian location
    GModelRadial::read(xml);

    // Extract model parameters
    int  npar[1]   = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

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
    if (npar[0] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Sigma\" parameters.");

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
 * Write the Gaussian source information into an XML element. The XML element
 * will have 3 parameter leafs named "RA", "DEC" and "Sigma". The location
 * leafs are handled by the GModelRadial base class.
 ***************************************************************************/
void GModelRadialGauss::write(GXmlElement& xml) const
{
    // Write shell location
    GModelRadial::write(xml);

    // If XML element has 2 nodes (which should be the location nodes)
    // then append 1 parameter node
    if (xml.elements() == 2) {
        xml.append(new GXmlElement("parameter name=\"Sigma\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || npars != 3)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Point source model requires exactly 3 parameters.");

    // Set or update model parameter attributes
    int npar[1] = {0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle Sigma
        if (par->attribute("name") == "Sigma") {
            m_sigma.write(*par);
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Sigma\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Gaussian source information
 ***************************************************************************/
std::string GModelRadialGauss::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelRadialGauss ===\n");
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
 *
 * Note that this model implements no gradients, as spatial models will
 * always require numerical gradient computations. The minimum Gaussian
 * width is set to 1 arcsec.
 ***************************************************************************/
void GModelRadialGauss::init_members(void)
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
    m_sigma.hasgrad(false);  // Radial components never have gradients

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
void GModelRadialGauss::copy_members(const GModelRadialGauss& model)
{
    // Copy members
    m_sigma = model.m_sigma;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelRadialGauss::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
