/***************************************************************************
 *  GModelSpatialElliptical.cpp - Abstract elliptical spatial model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialElliptical.cpp
 * @brief Abstract elliptical spatial model base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelSpatialElliptical.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                  "GModelSpatialElliptical::read(GXmlElement&)"
#define G_WRITE                "GModelSpatialElliptical::write(GXmlElement&)"

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
GModelSpatialElliptical::GModelSpatialElliptical(void) : GModelSpatial()
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
 * Constructs an elliptical spatial model component by extracting information
 * from an XML element. See read() for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialElliptical::GModelSpatialElliptical(const GXmlElement& xml) : 
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
 * @param[in] model Elliptical spatial model.
 ***************************************************************************/
GModelSpatialElliptical::GModelSpatialElliptical(const GModelSpatialElliptical& model) :
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
GModelSpatialElliptical::~GModelSpatialElliptical(void)
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
 * @param[in] model Elliptical spatial model.
 * @return Elliptical spatial model.
 ***************************************************************************/
GModelSpatialElliptical& GModelSpatialElliptical::operator=(const GModelSpatialElliptical& model)
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
 * Evaluates the elliptical spatial model for a given true photon arrival
 * direction.
 ***************************************************************************/
double GModelSpatialElliptical::eval(const GSkyDir& srcDir) const
{
    // Compute distance from source and position angle (in radians)
    double theta  = dir().dist(srcDir);
    double posang = dir().posang(srcDir);

    // Evaluate model
    double value = eval(theta, posang);

    // Return result
    return value;
}


/***********************************************************************//**
 * @brief Return model value and set analytical gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the elliptical spatial model for a given true photon arrival
 * direction.
 ***************************************************************************/
double GModelSpatialElliptical::eval_gradients(const GSkyDir& srcDir) const
{
    // Compute distance from source and position angle (in radians)
    double theta  = dir().dist(srcDir);
    double posang = dir().posang(srcDir);

    // Evaluate model and set gradients
    double value = eval_gradients(theta, posang);

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
 * Read the elliptical source location and position angle information from
 * an XML element in the following format
 *
 *     <spatialModel type="...">
 *       <parameter name="RA"  scale="1" value="83.63" min="-360" max="360" free="1"/>
 *       <parameter name="DEC" scale="1" value="22.01" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"  scale="1" value="45.0"  min="-360" max="360" free="1"/>
 *       ...
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="...">
 *       <parameter name="GLON" scale="1" value="83.63" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT" scale="1" value="22.01" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"   scale="1" value="45.0"  min="-360" max="360" free="1"/>
 *       ...
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialElliptical::read(const GXmlElement& xml)
{
    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has at least 3 parameters
    if (xml.elements() < 3 || npars < 3) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Elliptical model requires at least 3 parameters.");
    }

    // Extract model parameters
    bool has_glon = false;
    bool has_glat = false;
    int  npar[3]  = {0, 0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = static_cast<GXmlElement*>(xml.element("parameter", i));

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

        // Handle PA
        else if (par->attribute("name") == "PA") {
            m_posangle.read(*par);
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
              "Require either \"RA\"/\"DEC\" or \"GLON\"/\"GLAT\".");
    }

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"RA\"/\"DEC\", \"GLON\"/\"GLAT\" and \"PA\" parameters.");
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
 *            Existing XML element is not of requested type.
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the elliptical source location and position angle information into
 * an XML element in the following format
 *
 *     <spatialModel type="...">
 *       <parameter name="RA"  scale="1" value="83.63" min="-360" max="360" free="1"/>
 *       <parameter name="DEC" scale="1" value="22.01" min="-90"  max="90"  free="1"/>
 *       <parameter name="PA"  scale="1" value="45.0"  min="-360" max="360" free="1"/>
 *       ...
 *     </spatialModel>
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialElliptical::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Elliptical model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append 3 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"RA\""));
        xml.append(new GXmlElement("parameter name=\"DEC\""));
        xml.append(new GXmlElement("parameter name=\"PA\""));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has at least 3 parameters
    if (xml.elements() < 3 || npars < 3) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Elliptical source model requires at least 3 parameters.");
    }

    // Set or update model parameter attributes
    int npar[3] = {0, 0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = static_cast<GXmlElement*>(xml.element("parameter", i));

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

        // Handle PA
        else if (par->attribute("name") == "PA") {
            npar[2]++;
            m_posangle.write(*par);
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"RA\", \"DEC\" and \"PA\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return position of elliptical spatial model
 ***************************************************************************/
GSkyDir GModelSpatialElliptical::dir(void) const
{
    // Allocate sky direction
    GSkyDir srcDir;

    // Set sky direction
    srcDir.radec_deg(ra(), dec());

    // Return direction
    return srcDir;
}


/***********************************************************************//**
 * @brief Set position of elliptical spatial model
 ***************************************************************************/
void GModelSpatialElliptical::dir(const GSkyDir& dir)
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
void GModelSpatialElliptical::init_members(void)
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

    // Initialise Position Angle
    m_posangle.clear();
    m_posangle.name("PA");
    m_posangle.unit("deg");
    m_posangle.fix();
    m_posangle.scale(1.0);
    m_posangle.gradient(0.0);
    m_posangle.hasgrad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_posangle);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial spatial model.
 ***************************************************************************/
void GModelSpatialElliptical::copy_members(const GModelSpatialElliptical& model)
{
    // Copy members
    m_ra       = model.m_ra;
    m_dec      = model.m_dec;
    m_posangle = model.m_posangle;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);
    m_pars.push_back(&m_posangle);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialElliptical::free_members(void)
{
    // Return
    return;
}
