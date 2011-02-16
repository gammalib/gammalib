/***************************************************************************
 *     GModelSpatialShell.cpp  -  Spatial shell source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Christoph Deil                                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialShell.cpp
 * @brief GModelSpatialShell class implementation.
 * @author C. Deil
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialShell.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialShell    g_spatial_shell_seed;
const GModelSpatialRegistry g_spatial_shell_registry(&g_spatial_shell_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                       "GModelSpatialShell::read(GXmlElement&)"
#define G_WRITE                     "GModelSpatialShell::write(GXmlElement&)"

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
GModelSpatialShell::GModelSpatialShell(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Shell constructor
 *
 * @param[in] dir Sky position of shell centre.
 * @param[in] radius Inner shell radius (degrees).
 * @param[in] width Shell width (degrees).
 * @param [in] small_angle Use small angle approximation
 * (faster, but incorrect for very large shells (> a few degrees))
 ***************************************************************************/
GModelSpatialShell::GModelSpatialShell(const GSkyDir& dir,
                                       const double&  radius,
                                       const double&  width,
                                       const bool& small_angle) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Assign parameters
    this->dir(dir);
    this->radius(radius);
    this->width(width);
    this->small_angle(small_angle);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of shell spatial model by extracting information from an
 * XML element. See GModelSpatialShell::read() for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpatialShell::GModelSpatialShell(const GXmlElement& xml) : GModelSpatial()
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
 * @param[in] model Shell source spatial model.
 ***************************************************************************/
GModelSpatialShell::GModelSpatialShell(const GModelSpatialShell& model) :
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
GModelSpatialShell::~GModelSpatialShell(void)
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
 * @param[in] model Shell source spatial model.
 ***************************************************************************/
GModelSpatialShell& GModelSpatialShell::operator= (const GModelSpatialShell& model)
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
void GModelSpatialShell::clear(void)
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
GModelSpatialShell* GModelSpatialShell::clone(void) const
{
    return new GModelSpatialShell(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for a shell source model. The shell source
 * model is a radial function \f$f(\theta)\f$, where \f$\theta\f$ is the
 * angular separation between shell centre and the actual location.
 *
 * In the small angle approximation,
 * \f[
 * f(\theta) = {\tt m\_norm} \left \{
 *  \begin{array}{l l}
 *     \displaystyle
 *     \sqrt{ \theta_{\rm out}^2 - \theta^2 } -
       \sqrt{ \theta_{\rm in}^2  - \theta^2 }
 *     & \mbox{if $\theta \le \theta_{\rm in}$} \\
 *     \\
 *    \displaystyle
 *     \sqrt{ \theta_{\rm out}^2 - \theta^2 }
 *     & \mbox{if $\theta_{\rm in} < \theta \le \theta_{\rm out}$} \\
 *     \\
 *    \displaystyle
 *    0 & \mbox{if $\theta > \theta_{\rm out}$}
 *  \end{array}
 *  \right .
 * \f]
 * where \f${\tt m\_norm}\f$ is a normalization constant (see the update()
 * method).
 *
 * In the general form that is also valid for large angles
 * \f[
 * f(\theta) = {\tt m\_norm} \left \{
 *  \begin{array}{l l}
 *     \displaystyle
 *     \sqrt{ \sin^2 \theta_{\rm out} - \sin^2 \theta } -
       \sqrt{ \sin^2 \theta_{\rm in}  - \sin^2 \theta }
 *     & \mbox{if $\theta \le \theta_{\rm in}$} \\
 *     \\
 *    \displaystyle
 *     \sqrt{ \sin^2 \theta_{\rm out} - \sin^2 \theta }
 *     & \mbox{if $\theta_{\rm in} < \theta \le \theta_{\rm out}$} \\
 *     \\
 *    \displaystyle
 *    0 & \mbox{if $\theta > \theta_{\rm out}$}
 *  \end{array}
 *  \right .
 * \f]
 * 
 * Here, \f$\theta_{\rm in}\f$ and \f$\theta_{\rm out}\f$ are the shell inner
 * and outer radius.
 ***************************************************************************/
double GModelSpatialShell::eval(const GSkyDir& srcDir) const
{
    // Update precomputation cache
    update();

    // Compute distance from shell centre in radians
    double theta = dir().dist(srcDir);

    // Set x appropriately for the small angle approximation or not
    double x;
    if (m_small_angle)
    	x = theta * theta;
    else {
        x  = std::sin(theta);
    	x *= x;
    }

    // Compute value
  	double value = 0;
  	if (x < m_x_out) {
  		value = sqrt(m_x_out - x);
        if (x < m_x_in)
            value -= sqrt(m_x_in - x);
    }

    // Return normalised value
  	return (m_norm * value);
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * This method simply calls GModelSpatialShell::eval() as no analytical
 * gradients will be computed. See GModelSpatialShell::eval() for details
 * about the implemented method.
 ***************************************************************************/
double GModelSpatialShell::eval_gradients(const GSkyDir& srcDir) const
{
    // Return value
    return (eval(srcDir));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary sky position from the 2D shell distribution.
 *
 * @todo To be implemented.
 ***************************************************************************/
GSkyDir GModelSpatialShell::mc(GRan& ran) const
{
    // Simulate offset from photon arrival direction
	  // @todo How to draw theta for this model?
    double theta = (radius()+width()) * sqrt(ran.uniform());
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
 * is required to have 4 parameters.
 * The position is named either "RA" and "DEC" or "GLON" and "GLAT", the
 * shell inner and outer radius are named "Radius" and "Width".
 ***************************************************************************/
void GModelSpatialShell::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Shell model requires exactly 4 parameters.");

    // Extract model parameters
    bool has_glon = false;
    bool has_glat = false;
    int  npar[]   = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

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

        // Handle Radius / Width
        else if (par->attribute("name") == "Radius") {
            m_radius.read(*par);
            npar[2]++;
        }
        else if (par->attribute("name") == "Width") {
            m_width.read(*par);
            npar[3]++;
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

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"RA/GLON\", \"DEC/GLAT\", \"Radius\" and \"Width\" parameters.");

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
 * Write the Shell source information into an XML element. The XML element
 * will have 4 parameter leafs named "RA", "DEC", "Radius" and "Width".
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialShell::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", type());

    // Verify model type
    if (xml.attribute("type") != type())
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \""+type()+"\".");

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"RA\""));
        xml.append(new GXmlElement("parameter name=\"DEC\""));
        xml.append(new GXmlElement("parameter name=\"Radius\""));
        xml.append(new GXmlElement("parameter name=\"Width\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Point source model requires exactly 4 parameters.");

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle RA
        if (par->attribute("name") == "RA") {
            npar[0]++;
            m_ra.write(*par);
        }

        // Handle DEC
        else if (par->attribute("name") == "DEC") {
            npar[1]++;
            m_dec.write(*par);
        }

        // Handle Radius
        else if (par->attribute("name") == "Radius") {
            m_radius.write(*par);
            npar[2]++;
        }

        // Handle Width
        else if (par->attribute("name") == "Width") {
            m_width.write(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"RA\", \"DEC\", \"Radius\" and \"Width\" parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 ***************************************************************************/
std::string GModelSpatialShell::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialShell ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return position of point source
 ***************************************************************************/
GSkyDir GModelSpatialShell::dir(void) const
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
void GModelSpatialShell::dir(const GSkyDir& dir)
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
void GModelSpatialShell::init_members(void)
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

    // Initialise Radius
    m_radius.clear();
    m_radius.name("Radius");
    m_radius.unit("deg");
    m_radius.value(0.0);
    m_radius.free();
    m_radius.scale(1.0);
    m_radius.gradient(0.0);
    m_radius.hasgrad(false);

    // Initialise Width
    m_width.clear();
    m_width.name("Width");
    m_width.unit("deg");
    m_width.value(1.0);
    m_width.free();
    m_width.scale(1.0);
    m_width.gradient(0.0);
    m_width.hasgrad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_radius);
    m_pars.push_back(&m_width);

    // Initialise other members
    m_small_angle = true;

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius and width shell is not meaningful
    m_last_radius = 0.0;
    m_last_width  = 0.0;
    m_theta_in    = 0.0;
    m_x_in        = 0.0;
    m_theta_out   = 0.0;
    m_x_out       = 0.0;
    m_norm        = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Shell source spatial model.
 ***************************************************************************/
void GModelSpatialShell::copy_members(const GModelSpatialShell& model)
{
    // Copy members
    m_ra          = model.m_ra;
    m_dec         = model.m_dec;
    m_radius      = model.m_radius;
    m_width       = model.m_width;
    m_small_angle = model.m_small_angle;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_radius);
    m_pars.push_back(&m_width);

    // Copy precomputation cache
    m_last_radius = model.m_last_radius;
    m_last_width  = model.m_last_width;
    m_theta_in    = model.m_theta_in;
    m_x_in        = model.m_x_in;
    m_theta_out   = model.m_theta_out;
    m_x_out       = model.m_x_out;
    m_norm        = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialShell::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Performs precomputations that are valid for a given pair of radius and
 * width values. The following members are set by this method:
 * m_last_radius, m_last_width, m_theta_in, m_theta_out, m_x_in, m_x_out
 * and m_norm.
 *
 * m_theta_in contains the inner shell radius in radians, while m_theta_out
 * contains the outer shell radius in radians.
 *
 * In the small angle approximation, 
 * \f${\tt m\_x\_in} = {\tt m\_theta\_in}^2\f$,
 * \f${\tt m\_x\_out} = {\tt m\_theta\_out}^2\f$, and
 * \f[{\tt m\_norm} = \frac{3}{2 \pi}
 *    \frac{1}\left( {\tt m\_theta\_out}^3 - {\tt m\_theta\_in}^3 \right)\f].
 *
 * In the general case that is also valid for large angles,
 * \f${\tt m\_x\_in} = \sin^2 {\tt m\_theta\_in}\f$,
 * \f${\tt m\_x\_out} = \sin^2 {\tt m\_theta\_out}\f$, and
 * \f[{\tt m\_norm} = \frac{1}{2 \pi}
 *    \frac{\sqrt{1-\cos 2 {\tt m\_theta\_out}} -
 *          \sqrt{1-\cos 2 {\tt m\_theta\_in}}}{2 \sqrt{2}} +
 *    \frac{1+\cos 2 {\tt m\_theta\_out}}{4} \ln \left(
 *          \frac{\sqrt{2} \cos {\tt m\_theta\_out}}
 *               {\sqrt{2} + \sqrt{1 - \cos 2 {\tt m\_theta\_out}}} -
 *    \frac{1+\cos 2 {\tt m\_theta\_in}}{4} \ln \left(
 *          \frac{\sqrt{2} \cos {\tt m\_theta\_in}}
 *               {\sqrt{2} + \sqrt{1 - \cos 2 {\tt m\_theta\_in}}}\f]
 ***************************************************************************/
void GModelSpatialShell::update() const
{
    // Set constants
    const double c1 = twopi / 3.0;
    const double c2 = 1.0 / (2.0 * sqrt_two);

    // Update if radius or width have changed
    if (m_last_radius != radius() || m_last_width != width()) {

        // Store last values
        m_last_radius = radius();
        m_last_width  = width();

        // Perform precomputations
        m_theta_in   = radius()  * deg2rad;
        m_theta_out  = (radius() + width()) * deg2rad;
        if (m_small_angle) {
            m_x_in  = m_theta_in  * m_theta_in;
            m_x_out = m_theta_out * m_theta_out;
            m_norm  = 1.0 / (c1 * (m_x_out*m_theta_out - m_x_in*m_theta_in));
        } else {
            double sin_theta_in  = std::sin(m_theta_in);
            double sin_theta_out = std::sin(m_theta_out);
            double term1         = (f1(m_theta_out) - f1(m_theta_in)) * c2;
            double term2         = f2(m_theta_out);
            double term3         = f2(m_theta_in);
            m_norm  = 1.0 / (twopi * (term1 + term2 - term3));
            m_x_in  = sin_theta_in*sin_theta_in;
            m_x_out = sin_theta_out * sin_theta_out;
        }

    } // endif: update required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return function 1 value needed for precomputation
 *
 * Computes \f$f1(x) = \sqrt{1 - \cos 2 x}\f$.
 ***************************************************************************/
double GModelSpatialShell::f1(double x)
{
    // Compute value
    double f1 = std::sqrt(1.0 - std::cos(2.0 * x));

    // Return value
	return f1;
}


/***********************************************************************//**
 * @brief Return function 2 value needed for precomputation
 *
 * Compute
 * \f[f2(x) = \frac{1+\cos 2x}{4} 
 *    \ln \left( \frac{\sqrt{2} \cos x}{\sqrt{2} + \sqrt{ 1 - \cos 2 x}}
 *        \right)\f].
 ***************************************************************************/
double GModelSpatialShell::f2(double x)
{
    // Compute value
	double t1 = (1.0 + std::cos(2.0*x)) / 4.0;
	double t2 = sqrt_two * std::cos(x);
	double t3 = sqrt_two + f1(x);
	double f2 = t1 * std::log(t2 / t3);

    // Return value
	return f2;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
