/***************************************************************************
 *    GModelSpatialRadial.cpp - Abstract radial spatial model base class   *
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
 * @file GModelSpatialRadial.cpp
 * @brief Abstract radial spatial model base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelSpatialRadial.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                      "GModelSpatialRadial::read(GXmlElement&)"
#define G_WRITE                    "GModelSpatialRadial::write(GXmlElement&)"

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
GModelSpatialRadial::GModelSpatialRadial(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a radial spatial model by extracting information
 * from an XML element. See GModelSpatialRadial::read() for more information about
 * the expected structure of the XML element.
 ***************************************************************************/
GModelSpatialRadial::GModelSpatialRadial(const GXmlElement& xml) : 
                     GModelSpatial()
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
 * @param[in] model Radial spatial model.
 ***************************************************************************/
GModelSpatialRadial::GModelSpatialRadial(const GModelSpatialRadial& model) :
                     GModelSpatial()
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
GModelSpatialRadial::~GModelSpatialRadial(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Radial spatial model.
 ***************************************************************************/
GModelSpatialRadial& GModelSpatialRadial::operator=(const GModelSpatialRadial& model)
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return model value
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the radial spatial model for a given true photon arrival
 * direction.
 ***************************************************************************/
double GModelSpatialRadial::eval(const GSkyDir& srcDir) const
{
    // Compute distance from source (in radians)
    double theta = srcDir.dist(dir());

    // Evaluate model
    double value = eval(theta);

    // Return result
    return value;
}


/***********************************************************************//**
 * @brief Return model value and set analytical gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the radial spatial model for a given true photon arrival
 * direction.
 ***************************************************************************/
double GModelSpatialRadial::eval_gradients(const GSkyDir& srcDir) const
{
    // Compute distance from source (in radians)
    double theta = srcDir.dist(dir());

    // Evaluate model and set gradients
    double value = eval_gradients(theta);

    // Return result
    return value;
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
 * Read the radial source location information from an XML element. The XML
 * element is required to have at least 2 parameters.
 * The location is named either "RA" and "DEC" or "GLON" and "GLAT".
 ***************************************************************************/
void GModelSpatialRadial::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has at least 2 parameters
    if (xml.elements() < 2 || npars < 2) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Radial model requires at least 2 parameters.");
    }

    // Extract model parameters
    bool has_glon = false;
    bool has_glat = false;
    int  npar[2]  = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

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
              "Require either \"RA\"/\"DEC\" or \"GLON\"/\"GLAT\".");
    }

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"RA\"/\"DEC\" and \"GLON\"/\"GLAT\" parameters.");
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
 * Write the radial source location information into an XML element. The XML
 * element will have at least 2 parameter leafs named "RA" and "DEC".
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialRadial::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Radial model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append 2 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"RA\""));
        xml.append(GXmlElement("parameter name=\"DEC\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has at least 2 parameters
    if (xml.elements() < 2 || npars < 2) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Radial source model requires at least 2 parameters.");
    }

    // Set or update model parameter attributes
    int npar[2] = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

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

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"RA\" and \"DEC\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return position of radial spatial model
 ***************************************************************************/
GSkyDir GModelSpatialRadial::dir(void) const
{
    // Allocate sky direction
    GSkyDir srcDir;

    // Set sky direction
    srcDir.radec_deg(ra(), dec());

    // Return direction
    return srcDir;
}


/***********************************************************************//**
 * @brief Set position of radial spatial model
 ***************************************************************************/
void GModelSpatialRadial::dir(const GSkyDir& dir)
{
    // Assign Right Ascension and Declination
    m_ra.real_value(dir.ra_deg());
    m_dec.real_value(dir.dec_deg());

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialRadial::init_members(void)
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

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial spatial model.
 ***************************************************************************/
void GModelSpatialRadial::copy_members(const GModelSpatialRadial& model)
{
    // Copy members
    m_ra  = model.m_ra;
    m_dec = model.m_dec;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadial::free_members(void)
{
    // Return
    return;
}
