/***************************************************************************
 *    GModelSpatialRadial.cpp - Abstract radial spatial model base class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2022 by Juergen Knoedlseder                         *
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
 *       <parameter name="RA"  scale="1" value="83.6331" min="-360" max="360" free="0" />
 *       <parameter name="DEC" scale="1" value="22.0145" min="-90"  max="90"  free="0" />
 *       ...
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="...">
 *       <parameter name="GLON" scale="1" value="184.5575" min="-360" max="360" free="0" />
 *       <parameter name="GLAT" scale="1" value="-5.7843"  min="-90"  max="90"  free="0" />
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
        m_lon.read(*ra);
        m_lat.read(*dec);

    }

    // ... otherwise read GLON/GLAT parameters
    else {

        // Get parameters
        const GXmlElement* glon = gammalib::xml_get_par(G_READ, xml, "GLON");
        const GXmlElement* glat = gammalib::xml_get_par(G_READ, xml, "GLAT");

        // Read parameters
        m_lon.read(*glon);
        m_lat.read(*glat);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Depending on the coordinate system, write the radial source information
 * into an XML element with the following format
 *
 *     <spatialModel type="...">
 *       <parameter name="RA"  scale="1" value="83.6331" min="-360" max="360" free="0" />
 *       <parameter name="DEC" scale="1" value="22.0145" min="-90"  max="90"  free="0" />
 *       ...
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="PointSource">
 *       <parameter name="GLON" scale="1" value="184.5575" min="-360" max="360" free="0" />
 *       <parameter name="GLAT" scale="1" value="-5.7843"  min="-90"  max="90"  free="0" />
 *       ...
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialRadial::write(GXmlElement& xml) const
{
    // Check model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get or create parameters
    GXmlElement* lon = gammalib::xml_need_par(G_WRITE, xml, m_lon.name());
    GXmlElement* lat = gammalib::xml_need_par(G_WRITE, xml, m_lat.name());

    // Write parameters
    m_lon.write(*lon);
    m_lat.write(*lat);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return position of radial spatial model
 *
 * @return Radial spatial model sky direction.
 *
 * Returns the sky direction of the radial spatial model.
 ***************************************************************************/
const GSkyDir& GModelSpatialRadial::dir(void) const
{
    // Get longitude and latitude values
    double lon = m_lon.value();
    double lat = m_lat.value();

    // If longitude or latitude values have changed then update sky
    // direction cache
    if ((lon != m_last_lon) || (lat != m_last_lat)) {

        // Update last values
        m_last_lon = lon;
        m_last_lat = lat;

        // Update sky direction dependent on model coordinate system
        if (is_celestial()) {
            m_dir.radec_deg(m_last_lon, m_last_lat);
        }
        else {
            m_dir.lb_deg(m_last_lon, m_last_lat);
        }

    } // endif: update of sky direction cache required

    // Return sky direction
    return (m_dir);
}


/***********************************************************************//**
 * @brief Set position of radial spatial model
 *
 * @param[in] dir Sky direction of radial spatial model.
 *
 * Sets the sky direction of the radial spatial model.
 ***************************************************************************/
void GModelSpatialRadial::dir(const GSkyDir& dir)
{
    // Assign sky direction depending on the model coordinate system
    if (is_celestial()) {
        m_lon.value(dir.ra_deg());
        m_lat.value(dir.dec_deg());
    }
    else {
        m_lon.value(dir.l_deg());
        m_lat.value(dir.b_deg());
    }

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
    m_lon.clear();
    m_lon.name("RA");
    m_lon.unit("deg");
    m_lon.fix();
    m_lon.scale(1.0);
    m_lon.gradient(0.0);
    m_lon.has_grad(false);

    // Initialise Declination
    m_lat.clear();
    m_lat.name("DEC");
    m_lat.unit("deg");
    m_lat.fix();
    m_lat.scale(1.0);
    m_lat.gradient(0.0);
    m_lat.has_grad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_lon);
    m_pars.push_back(&m_lat);

    // Initialise cache
    m_dir.clear();
    m_last_lon = -9999.0;
    m_last_lat = -9999.0;

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
    m_lon = model.m_lon;
    m_lat = model.m_lat;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_lon);
    m_pars.push_back(&m_lat);

    // Copy cache
    m_dir      = model.m_dir;
    m_last_lon = model.m_last_lon;
    m_last_lat = model.m_last_lat;

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
