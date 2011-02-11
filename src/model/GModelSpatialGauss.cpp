/***************************************************************************
 *     GModelSpatialGauss.cpp  -  Spatial Gaussian source model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialGauss.cpp
 * @brief Gaussian spatial model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialGauss.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialGauss    g_spatial_gauss_seed;
const GModelSpatialRegistry g_spatial_gauss_registry(&g_spatial_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                       "GModelSpatialGauss::read(GXmlElement&)"
#define G_WRITE                     "GModelSpatialGauss::write(GXmlElement&)"

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
GModelSpatialGauss::GModelSpatialGauss(void) : GModelSpatial()
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
GModelSpatialGauss::GModelSpatialGauss(const GSkyDir& dir,
                                       const double& sigma) : GModelSpatial()
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
 * from an XML element. See GModelSpatialGauss::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelSpatialGauss::GModelSpatialGauss(const GXmlElement& xml) : GModelSpatial()
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
 * @param[in] model Gaussian spatial model.
 ***************************************************************************/
GModelSpatialGauss::GModelSpatialGauss(const GModelSpatialGauss& model) :
                                       GModelSpatial(model)
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
GModelSpatialGauss::~GModelSpatialGauss(void)
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
 * @param[in] model Gaussian spatial model.
 ***************************************************************************/
GModelSpatialGauss& GModelSpatialGauss::operator=(const GModelSpatialGauss& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

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
void GModelSpatialGauss::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpatialGauss* GModelSpatialGauss::clone(void) const
{
    return new GModelSpatialGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for a Gaussian source model. The Gaussian
 * source model is defined as
 * \f[f(\theta)=\frac{1}{2 \pi \sigma^2} \exp 
 *    \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)\f]
 * where
 * \f$\theta\f$ is the angular separation from the source direction, and
 * \f$\sigma\f$ is the Gaussian width.
 ***************************************************************************/
double GModelSpatialGauss::eval(const GSkyDir& srcDir) const
{
    // Compute distance from source
    double theta = srcDir.dist(dir());
    
    // Compute value
    double sigma_rad = sigma() * deg2rad;
    double sigma2    = sigma_rad * sigma_rad;
    double theta2    = theta   * theta;
    double value     = exp(-0.5 * theta2 / sigma2) / (twopi * sigma2);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for a Gaussian source model and the gradient
 * with respect to the source extent.  The Gaussian source model is defined
 * as
 * \f[f(\theta)=\frac{1}{2 \pi \sigma^2} \exp 
 *    \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)\f]
 * where
 * \f$\theta\f$ is the angular separation from the source direction, and
 * \f$\sigma\f$ is the Gaussian width.
 *
 * The partial derivative of the Gaussian width is given by
 * \f[\frac{df}{d\sigma_v} = 
 *    f(\theta) (\frac{\theta^2}{\sigma^3} - \frac{2}{\sigma}) \sigma_s\f]
 * where 
 * \f$\sigma_v\f$ is the value part, 
 * \f$\sigma_s\f$ is the scaling part, and 
 * \f$\sigma = \sigma_v \sigma_s\f$. 
 *
 * This method only provides the parameter gradient for the Gaussian width
 * \f$\sigma\f$. No gradients are provided for the position of the Gaussian.
 ***************************************************************************/
double GModelSpatialGauss::eval_gradients(const GSkyDir& srcDir) const
{
    // Compute distance from source (in radians)
    double theta = srcDir.dist(dir());
    
    // Compute value (using radians)
    double sigma_rad = sigma() * deg2rad;
    double sigma2    = sigma_rad * sigma_rad;
    double theta2    = theta * theta;
    double arg       = theta2 / sigma2;
    double value     = exp(-0.5 * arg) / (twopi * sigma2);

    // Compute partial derivatives of the sigma parameter. Note that this
    // parameter is given in degrees
    double g_sigma = value * (arg - 2.0) / sigma() * m_sigma.scale();

    // Set gradients
    ((GModelSpatialGauss*)this)->m_sigma.gradient(g_sigma);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary sky position from the 2D Gaussian distribution.
 ***************************************************************************/
GSkyDir GModelSpatialGauss::mc(GRan& ran) const
{
    // Simulate offset from photon arrival direction
    double theta = sigma() * ran.chisq2();
    double phi   = 360.0 * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate(phi, theta);
    
    // Return sky direction
    return sky_dir;
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
 * Gaussian width is named "Sigma".
 ***************************************************************************/
void GModelSpatialGauss::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Point source model requires exactly 3 parameters.");

    // Extract model parameters
    bool has_glon = false;
    bool has_glat = false;
    int  npar[]   = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle RA/GLON
        if (par->attribute("name") == "RA") {
            m_ra.read(*par);
            npar[0]++;
        }
        else if (par->attribute("name") == "GLON") {
            m_ra.read(*par);
            npar[0]++;
            has_glon = true;
        }

        // Handle DEC/GLAT
        else if (par->attribute("name") == "DEC") {
            m_dec.read(*par);
            npar[1]++;
        }
        else if (par->attribute("name") == "GLAT") {
            m_dec.read(*par);
            npar[1]++;
            has_glat = true;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Sigma") {
            m_sigma.read(*par);
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Check if we have to convert GLON/GLAT into RA/DEC
    if (has_glon && has_glat) {
        GSkyDir dir;
        dir.lb_deg(ra(), dec()),
        m_ra.real_value(dir.ra_deg());
        m_dec.real_value(dir.dec_deg());
    }
    else if (has_glon || has_glat) {
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require either RA/DEC or GLON/GLAT.");
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
 * Write the Gaussian source information into an XML element. The XML element
 * has to be of type 'GaussFunction' and will have 3 parameter leafs
 * named "RA", "DEC" and "Sigma".
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialGauss::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", type());

    // Verify model type
    if (xml.attribute("type") != type())
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \""+type()+"\".");

    // If XML element has 0 nodes then append 3 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"RA\""));
        xml.append(new GXmlElement("parameter name=\"DEC\""));
        xml.append(new GXmlElement("parameter name=\"Sigma\""));
    }

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Point source model requires exactly 3 parameters.");

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "RA") {
            npar[0]++;
            m_ra.write(*par);
        }

        // Handle index
        else if (par->attribute("name") == "DEC") {
            npar[1]++;
            m_dec.write(*par);
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Sigma") {
            m_sigma.write(*par);
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"RA\", \"DEC\" and \"Sigma\" parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Gaussian source information
 ***************************************************************************/
std::string GModelSpatialGauss::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialGauss ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return position of point source
 ***************************************************************************/
GSkyDir GModelSpatialGauss::dir(void) const
{
    // Allocate sky direction
    GSkyDir srcDir;

    // Set sky direction
    srcDir.radec_deg(ra(), dec());

    // Return direction
    return srcDir;
}


/***********************************************************************//**
 * @brief Set position of point source
 ***************************************************************************/
void GModelSpatialGauss::dir(const GSkyDir& dir)
{
    // Assign Right Ascension and Declination
    m_ra.real_value(dir.ra_deg());
    m_dec.real_value(dir.dec_deg());

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialGauss::init_members(void)
{
    // Initialise Right Ascension
    m_ra.clear();
    m_ra.name("RA");
    m_ra.unit("deg");
    m_ra.fix();
    m_ra.scale(1.0);
    m_ra.gradient(0.0);
    m_ra.hasgrad(false);

    // Initialise Declination
    m_dec.clear();
    m_dec.name("DEC");
    m_dec.unit("deg");
    m_dec.fix();
    m_dec.scale(1.0);
    m_dec.gradient(0.0);
    m_dec.hasgrad(false);

    // Initialise Gaussian sigma
    m_sigma.clear();
    m_sigma.name("Sigma");
    m_sigma.unit("deg");
    m_sigma.free();
    m_sigma.scale(1.0);
    m_sigma.gradient(0.0);
    m_sigma.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Gaussian spatial model.
 ***************************************************************************/
void GModelSpatialGauss::copy_members(const GModelSpatialGauss& model)
{
    // Copy members
    m_ra    = model.m_ra;
    m_dec   = model.m_dec;
    m_sigma = model.m_sigma;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialGauss::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
