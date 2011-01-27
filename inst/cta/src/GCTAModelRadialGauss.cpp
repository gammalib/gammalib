/***************************************************************************
 *      GCTAModelRadialGauss.cpp  -  Radial Gaussian CTA model class       *
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
 * @file GCTAModelRadialGauss.cpp
 * @brief GCTAModelRadialGauss class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GCTAModelRadialGauss.hpp"
#include "GCTAModelRadialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRadialGauss    g_cta_radial_gauss_seed;
const GCTAModelRadialRegistry g_cta_radial_gauss_registry(&g_cta_radial_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                     "GCTAModelRadialGauss::operator() (int)"
#define G_READ                     "GCTAModelRadialGauss::read(GXmlElement&)"
#define G_WRITE                   "GCTAModelRadialGauss::write(GXmlElement&)"

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
GCTAModelRadialGauss::GCTAModelRadialGauss(void) : GCTAModelRadial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] sigma Gaussian width (degrees).
 ***************************************************************************/
GCTAModelRadialGauss::GCTAModelRadialGauss(const double& sigma) : GCTAModelRadial()
{
    // Initialise members
    init_members();

    // Assign sigma
    this->sigma(sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GCTAModelRadialGauss::GCTAModelRadialGauss(const GXmlElement& xml) : GCTAModelRadial()
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
 * @param[in] model Model.
 ***************************************************************************/
GCTAModelRadialGauss::GCTAModelRadialGauss(const GCTAModelRadialGauss& model) :
                                           GCTAModelRadial(model)
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
GCTAModelRadialGauss::~GCTAModelRadialGauss(void)
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
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar& GCTAModelRadialGauss::operator() (int index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS, index, 0, size()-1);
    #endif

    // Return pointer
    return *(m_par[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
const GModelPar& GCTAModelRadialGauss::operator() (int index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS, index, 0, size()-1);
    #endif

    // Return pointer
    return *(m_par[index]);
}


/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Model.
 ***************************************************************************/
GCTAModelRadialGauss& GCTAModelRadialGauss::operator= (const GCTAModelRadialGauss& model)
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
void GCTAModelRadialGauss::clear(void)
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
GCTAModelRadialGauss* GCTAModelRadialGauss::clone(void) const
{
    return new GCTAModelRadialGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] offset Offset angle (degrees).
 *
 * Evaluates the radial part for a Gaussian model. The Gaussian model is
 * defined as
 * \f[f(\theta)=\exp \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$\sigma\f$ is the Gaussian width (in degrees).
 *
 * Note that this method implements a Gaussian function which is unity for
 * \f$\theta=0\f$.
 ***************************************************************************/
double GCTAModelRadialGauss::eval(const double& offset)
{
    // Compute value
    double arg   = offset / sigma();
    double value = exp(-0.5 * arg * arg);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] offset Offset angle (degrees).
 *
 * Evaluates the radial part for a Gaussian model. The Gaussian model is
 * defined as
 * \f[f(\theta)=\exp \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right)\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$\sigma\f$ is the Gaussian width (in degrees).
 *
 * The partial derivative of the Gaussian width is given by
 * \f[\frac{df}{d\sigma_v} = f(\theta) \frac{\theta^2}{\sigma^3} \sigma_s\f]
 * where 
 * \f$\sigma_v\f$ is the value part, 
 * \f$\sigma_s\f$ is the scaling part, and 
 * \f$\sigma = \sigma_v \sigma_s\f$. 
 *
 * Note that this method implements a Gaussian function which is unity for
 * \f$\theta=0\f$.
 ***************************************************************************/
double GCTAModelRadialGauss::eval_gradients(const double& offset)
{
    // Compute value
    double arg   = offset / sigma();
    double arg2  = arg * arg;
    double value = exp(-0.5 * arg2);

    // Compute partial derivatives of the sigma parameter.
    double g_sigma = value * arg2 / sigma() * m_sigma.scale();

    // Set gradients
    m_sigma.gradient(g_sigma);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] dir Pointing direction.
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary CTA instrument position from a 2D Gaussian
 * distribution.
 ***************************************************************************/
GCTAInstDir GCTAModelRadialGauss::mc(const GCTAInstDir& dir, GRan& ran) const
{
    // Simulate offset from photon arrival direction
    double offset = sigma() * ran.chisq2();
    double phi    = 360.0 * ran.uniform();

    // Rotate pointing direction by offset and azimuth angle
    GCTAInstDir mc_dir = dir;
    mc_dir.rotate(phi, offset);
    
    // Return MC direction
    return mc_dir;
}


/***********************************************************************//**
 * @brief Returns integral over radial model (in steradians)
 *
 * Computes
 * \f[\Omega = 2 \pi \int f(\theta) sin(\theta) d\theta\f]
 * where
 * \f$f(\theta)\f$ is the radial Gaussian function, and
 * \f$\theta\f$ is the offset angle.
 ***************************************************************************/
double GCTAModelRadialGauss::omega(void) const
{
    // Analytically compute integral over the Gaussian
    double sigma_rad = sigma() * deg2rad;
    double omega     = twopi * sigma_rad * sigma_rad;
    
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
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the Gaussian radial model information from an XML element. The XML
 * element is required to have 1 parameter named 'Sigma'. 
 ***************************************************************************/
void GCTAModelRadialGauss::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Radial Gaussian model requires exactly 1 parameter.");

    // Extract model parameters
    int  npar[]   = {0};
    for (int i = 0; i < 1; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle Sigma
        if (par->attribute("name") == "Sigma") {
            m_sigma.read(*par);
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Radial Gaussian model requires \"Sigma\" parameter.");

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
 * Write the Gaussian radial model information into an XML element. The XML
 * element will have 1 parameter leaf named 'Sigma'.
 ***************************************************************************/
void GCTAModelRadialGauss::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", type());

    // Verify model type
    if (xml.attribute("type") != type())
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Radial Gaussian model is not of type \""+type()+"\".");

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"Sigma\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Radial Gaussian model requires exactly 1 parameter.");

    // Set or update model parameter attributes
    int npar[] = {0};
    for (int i = 0; i < 1; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle sigma
        if (par->attribute("name") == "Sigma") {
            m_sigma.write(*par);
            npar[0]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Radial Gaussian model requires \"Sigma\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 ***************************************************************************/
std::string GCTAModelRadialGauss::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAModelRadialGauss ===");
    result.append("\n"+parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_par[i]->print());

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
void GCTAModelRadialGauss::init_members(void)
{
    // Initialise parameters
    m_npars  = 1;
    m_par[0] = &m_sigma;

    // Initialise Gaussian sigma
    m_sigma = GModelPar();
    m_sigma.name("Sigma");
    m_sigma.unit("deg");
    m_sigma.free();
    m_sigma.scale(1.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GCTAModelRadialGauss::copy_members(const GCTAModelRadialGauss& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_sigma = model.m_sigma;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRadialGauss::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
