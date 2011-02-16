/***************************************************************************
 *     GModelSpatialDisk.cpp  -  Spatial disk source model class           *
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
 * @file GModelSpatialDisk.cpp
 * @brief GModelSpatialDisk class implementation.
 * @author C. Deil
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialDisk.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDisk     g_spatial_disk_seed;
const GModelSpatialRegistry g_spatial_disk_registry(&g_spatial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                       "GModelSpatialDisk::read(GXmlElement&)"
#define G_WRITE                     "GModelSpatialDisk::write(GXmlElement&)"

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
GModelSpatialDisk::GModelSpatialDisk(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Disk constructor
 *
 * @param[in] dir Sky position of shell centre.
 * @param[in] radius Disk radius (degrees).
 ***************************************************************************/
GModelSpatialDisk::GModelSpatialDisk(const GSkyDir& dir,
                                     const double&  radius) : GModelSpatial()
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
 * Creates instance of disk spatial model by extracting information from an
 * XML element. See GModelSpatialDisk::read() for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpatialDisk::GModelSpatialDisk(const GXmlElement& xml) : GModelSpatial()
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
 * @param[in] model Disk source spatial model.
 ***************************************************************************/
GModelSpatialDisk::GModelSpatialDisk(const GModelSpatialDisk& model) :
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
GModelSpatialDisk::~GModelSpatialDisk(void)
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
 * @param[in] model Disk source spatial model.
 ***************************************************************************/
GModelSpatialDisk& GModelSpatialDisk::operator= (const GModelSpatialDisk& model)
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
void GModelSpatialDisk::clear(void)
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
GModelSpatialDisk* GModelSpatialDisk::clone(void) const
{
    return new GModelSpatialDisk(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
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
double GModelSpatialDisk::eval(const GSkyDir& srcDir) const
{
    // Update precomputation cache
    update();

    // Compute angular distance from disk centre
    double theta = dir().dist_deg(srcDir);

    // Set value
    double value = (theta <= radius()) ? m_norm : 0.0;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcDir True photon arrival direction.
 ***************************************************************************/
double GModelSpatialDisk::eval_gradients(const GSkyDir& srcDir) const
{
    // Update precomputation cache
    update();

    // Compute angular distance from disk centre
    double theta = dir().dist_deg(srcDir);

    // Set value
    double value = (theta <= radius()) ? m_norm : 0.0;

    // Compute partial derivative of radius
    double g_radius = (m_radius.isfree()) ? m_gradient : 0.0;

    // Set gradient (circumvent const correctness)
    ((GModelSpatialDisk*)this)->m_radius.gradient(g_radius);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary sky position from the 2D disk distribution.
 ***************************************************************************/
GSkyDir GModelSpatialDisk::mc(GRan& ran) const
{
    // Simulate offset from photon arrival direction
    double theta = std::acos(1.0 - ran.uniform() *
                             (1.0 - std::cos(radius()*deg2rad))) * rad2deg;
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
 * disk radius is named "Radius".
 ***************************************************************************/
void GModelSpatialDisk::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Disk model requires exactly 3 parameters.");

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

        // Handle Radius
        else if (par->attribute("name") == "Radius") {
            m_radius.read(*par);
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

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"RA/GLON\", \"DEC/GLAT\" and \"Radius\" parameters.");

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
 * Write the Disk source information into an XML element. The XML element
 * will have 3 parameter leafs named "RA", "DEC" and "Radius"
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialDisk::write(GXmlElement& xml) const
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
        xml.append(new GXmlElement("parameter name=\"Radius\""));
    }

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Disk model requires exactly 3 parameters.");

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

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

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"RA\", \"DEC\", \"Radius\" parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 ***************************************************************************/
std::string GModelSpatialDisk::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialDisk ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return centre of disk source
 ***************************************************************************/
GSkyDir GModelSpatialDisk::dir(void) const
{
    // Allocate sky direction
    GSkyDir srcDir;

    // Set sky direction
    srcDir.radec_deg(ra(), dec());

    // Return direction
    return srcDir;
}


/***********************************************************************//**
 * @brief Set centre of disk source
 ***************************************************************************/
void GModelSpatialDisk::dir(const GSkyDir& dir)
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
void GModelSpatialDisk::init_members(void)
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
    m_radius.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_radius);

    // Initialise precomputation cache. Note that zero values flag
    // uninitialised as a zero radius and width shell is not meaningful
    m_last_radius = 0.0;
    m_norm        = 0.0;
    m_gradient    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Disk source spatial model.
 ***************************************************************************/
void GModelSpatialDisk::copy_members(const GModelSpatialDisk& model)
{
    // Copy members
    m_ra     = model.m_ra;
    m_dec    = model.m_dec;
    m_radius = model.m_radius;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_radius);

    // Copy precomputation cache
    m_last_radius = model.m_last_radius;
    m_norm        = model.m_norm;
    m_gradient    = model.m_gradient;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialDisk::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 *
 * Computes the normalization
 * \f[{\tt m\_norm} = \frac{1}{2 \pi (1 - \cos r)}\f]
 * and the derivative with respect to the radius value
 * \f[\frac{\delta{\tt m\_norm}}{\delta r_v} = -2 \pi {\tt m\_norm}^2 
 *    \sin r \times r_s\f]
 ***************************************************************************/
void GModelSpatialDisk::update() const
{
    // Update if radius has changed
    if (m_last_radius != radius()) {

        // Store last values
        m_last_radius = radius();

        // Perform precomputations
        double r   = radius() * deg2rad;
        m_norm     = 1. / (twopi * (1 - std::cos(r)));
        m_gradient = -twopi * m_norm * m_norm *std::sin(r) * m_radius.scale();
        
    } // endif: update required

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
