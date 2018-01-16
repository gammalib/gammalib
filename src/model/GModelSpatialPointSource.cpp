/***************************************************************************
 *     GModelSpatialPointSource.cpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialPointSource.cpp
 * @brief Point source spatial model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */
const double tolerance =  0.000027777778;    // angular tolerance is 1 arcsec


/* __ Globals ____________________________________________________________ */
const GModelSpatialPointSource g_spatial_ptsrc_seed;
const GModelSpatialRegistry    g_spatial_ptsrc_registry(&g_spatial_ptsrc_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialPointSource g_spatial_ptsrc_legacy_seed(true, "SkyDirFunction");
const GModelSpatialRegistry    g_spatial_ptsrc_legacy_registry(&g_spatial_ptsrc_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GModelSpatialPointSource::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialPointSource::write(GXmlElement&)"

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
 *
 * Constructs empty point source model.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model type constructor
 *
 * @param[in] dummy Dummy flag.
 * @param[in] type Model type.
 *
 * Constructs empty point source model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const bool&        dummy,
                                                   const std::string& type) :
                          GModelSpatial()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky direction constructor
 *
 * @param[in] dir Sky direction.
 *
 * Construct a point source spatial model from a sky direction.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const GSkyDir& dir) :
                          GModelSpatial()
{
    // Initialise members
    init_members();

    // Assign direction
    this->dir(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Value constructor
 *
 * @param[in] ra Right Ascencion of model centre.
 * @param[in] dec Declination of model centre.
 *
 * Construct a point source spatial model from the Right Ascension and
 * Declination of the model centre.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const double& ra,
                                                   const double& dec) :
                          GModelSpatial()
{
    // Initialise members
    init_members();

    // Set values
    m_ra.value(ra);
    m_dec.value(dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Construct a point source spatial model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const GXmlElement& xml) :
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
 * @param[in] model Point source spatial model.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const GModelSpatialPointSource& model) :
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
GModelSpatialPointSource::~GModelSpatialPointSource(void)
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
 * @param[in] model Point source spatial model.
 * @return Point source spatial model.
 ***************************************************************************/
GModelSpatialPointSource& GModelSpatialPointSource::operator=(const GModelSpatialPointSource& model)
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
 * @brief Clear point source model
 ***************************************************************************/
void GModelSpatialPointSource::clear(void)
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
 * @brief Clone point source model
 *
 * @return Pointer to deep copy of point source model.
 ***************************************************************************/
GModelSpatialPointSource* GModelSpatialPointSource::clone(void) const
{
    // Clone point source model
    return new GModelSpatialPointSource(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] photon Incident photon.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * Evaluates the spatial part for a point source model. It implements a delta
 * function with respect to the coordinates of the source. For numerical
 * reasons, a certain tolerance is accepted (typically 0.1 arcsec, i.e. well
 * below the angular resolution of gamma-ray telescopes).
 *
 * The method will not compute analytical parameter gradients, even if the
 * @p gradients argument is set to true. Point source parameter gradients
 * need to be computed numerically.
 ***************************************************************************/
double GModelSpatialPointSource::eval(const GPhoton& photon,
                                      const bool&    gradients) const
{
    // Set value dependent on source distance
    double value = (photon.dir().dist_deg(dir()) < tolerance) ? 1.0 : 0.0;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction for the point source model. As the point
 * source is a point in the sky, the methods always returns the directon of
 * the point source.
 ***************************************************************************/
GSkyDir GModelSpatialPointSource::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Return sky direction
    return (dir());
}


/***********************************************************************//**
 * @brief Checks where model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the point source
 * model.
 ***************************************************************************/
bool GModelSpatialPointSource::contains(const GSkyDir& dir,
                                        const double&  margin) const
{
    // Compute distance to centre (radian)
    double distance = dir.dist(this->dir());

    // Return flag
    return (distance <= margin*gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing point source model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the point source information from an XML element with the following
 * format
 *
 *     <spatialModel type="PointSource">
 *       <parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="PointSource">
 *       <parameter free="0" max="360" min="-360" name="GLON" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="GLAT" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialPointSource::read(const GXmlElement& xml)
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
 *            Existing XML element is not of type 'SkyDirFunction'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the point source information into an XML element with the following
 * format
 *
 *     <spatialModel type="PointSource">
 *       <parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialPointSource::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \""+type()+"\".");
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
 * @brief Print point source information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point source information.
 ***************************************************************************/
std::string GModelSpatialPointSource::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialPointSource ===");

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


/***********************************************************************//**
 * @brief Return position of point source
 *
 * @return Point source sky direction.
 *
 * Returns the sky direction of the point source.
 ***************************************************************************/
GSkyDir GModelSpatialPointSource::dir(void) const
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
 *
 * Sets the sky direction of the point source.
 ***************************************************************************/
void GModelSpatialPointSource::dir(const GSkyDir& dir)
{
    // Assign Right Ascension and Declination
    m_ra.value(dir.ra_deg());
    m_dec.value(dir.dec_deg());

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
void GModelSpatialPointSource::init_members(void)
{
    // Initialise model type
    m_type = "PointSource";

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

    // Initialise other members
    m_region.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Point source spatial model.
 ***************************************************************************/
void GModelSpatialPointSource::copy_members(const GModelSpatialPointSource& model)
{
    // Copy members
    m_type   = model.m_type;
    m_ra     = model.m_ra;
    m_dec    = model.m_dec;
    m_region = model.m_region;

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
void GModelSpatialPointSource::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialPointSource::set_region(void) const
{
    // Set sky region centre
    m_region.centre(m_ra.value(), m_dec.value());

    // Set sky region radius to zero
    m_region.radius(0.0);

    // Return
    return;
}
