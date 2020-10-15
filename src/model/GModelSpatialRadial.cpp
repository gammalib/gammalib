/***************************************************************************
 *    GModelSpatialRadial.cpp - Abstract radial spatial model base class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * Constructs a radial spatial model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
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
 * @return Radial spatial model.
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
 * @param[in] photon Incident Photon.
 * @param[in] gradients Compute gradients?
 * @return Value of spatial radial model.
 *
 * Evaluates the radial spatial model value for a specific incident
 * @p photon.
 *
 * If the @p gradients flag is true the method will also compute the
 * parameter gradients for all model parameters.
 ***************************************************************************/
double GModelSpatialRadial::eval(const GPhoton& photon,
                                 const bool&    gradients) const
{
    // Compute distance from source (in radians)
    double theta = photon.dir().dist(dir());

    // Evaluate model
    double value = eval(theta, photon.energy(), photon.time(), gradients);

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
 * Reads the radial source location and position angle information from an
 * XML element in the following format
 *
 *     <spatialModel type="...">
 *       <parameter name="RA"  scale="1" value="83.63" min="-360" max="360" free="1"/>
 *       <parameter name="DEC" scale="1" value="22.01" min="-90"  max="90"  free="1"/>
 *       ...
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="...">
 *       <parameter name="GLON" scale="1" value="83.63" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT" scale="1" value="22.01" min="-90"  max="90"  free="1"/>
 *       ...
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadial::read(const GXmlElement& xml)
{
    // Read RA/DEC parameters
    if (gammalib::xml_has_par(xml, "RA") && gammalib::xml_has_par(xml, "DEC")) {

        // Get parameters
        const GXmlElement* ra  = gammalib::xml_get_par(G_READ, xml, "RA");
        const GXmlElement* dec = gammalib::xml_get_par(G_READ, xml, "DEC");

        // Read parameters
        m_ra.read(*ra);
        m_dec.read(*dec);

    }

    // ... otherwise read GLON/GLAT parameters
    else {

        // Get parameters
        const GXmlElement* glon = gammalib::xml_get_par(G_READ, xml, "GLON");
        const GXmlElement* glat = gammalib::xml_get_par(G_READ, xml, "GLAT");

        // Read parameters
        m_ra.read(*glon);
        m_dec.read(*glat);

        // Convert into RA/DEC
        GSkyDir dir;
        dir.lb_deg(ra(), dec()),
        m_ra.value(dir.ra_deg());
        m_dec.value(dir.dec_deg());

        // Set names to RA/DEC
        m_ra.name("RA");
        m_dec.name("DEC");

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
 * Writes the radial source location and position angle information into an
 * XML element in the following format
 *
 *     <spatialModel type="...">
 *       <parameter name="RA"  scale="1" value="83.63" min="-360" max="360" free="1"/>
 *       <parameter name="DEC" scale="1" value="22.01" min="-90"  max="90"  free="1"/>
 *       ...
 *     </spatialModel>
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

    // Get or create parameters
    GXmlElement* ra  = gammalib::xml_need_par(G_WRITE, xml, m_ra.name());
    GXmlElement* dec = gammalib::xml_need_par(G_WRITE, xml, m_dec.name());

    // Write parameters
    m_ra.write(*ra);
    m_dec.write(*dec);

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
    m_ra.value(dir.ra_deg());
    m_dec.value(dir.dec_deg());

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
    m_ra.has_grad(false);

    // Initialise Declination
    m_dec.clear();
    m_dec.name("DEC");
    m_dec.unit("deg");
    m_dec.fix();
    m_dec.scale(1.0);
    m_dec.gradient(0.0);
    m_dec.has_grad(false);

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
