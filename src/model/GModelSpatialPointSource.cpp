/***************************************************************************
 *     GModelSpatialPointSource.cpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2022 by Juergen Knoedlseder                         *
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
#define G_CONSTRUCTOR1  "GModelSpatialPointSource::GModelSpatialPointSource("\
                                                    "GSkyDir&, std::string&)"
#define G_CONSTRUCTOR2  "GModelSpatialPointSource::GModelSpatialPointSource("\
                                            "double&, double&, std::string&)"
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
 * @param[in] coordsys Coordinate system (either "CEL" or "GAL")
 *
 * @exception GException::invalid_argument
 *            Invalid @p coordsys argument specified.
 *
 * Construct a point source spatial model from a sky direction. The
 * @p coordsys parameter specified whether the sky direction should be
 * interpreted in the celestial or Galactic coordinata system.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const GSkyDir&     dir,
                                                   const std::string& coordsys) :
                          GModelSpatial()
{
    // Throw an exception if the coordinate system is invalid
    if ((coordsys != "CEL") && (coordsys != "GAL")) {
        std::string msg = "Invalid coordinate system \""+coordsys+"\" "
                          "specified. Please specify either \"CEL\" or "
                          "\"GAL\".";
        throw GException::invalid_argument(G_CONSTRUCTOR1, msg);
    }

    // Initialise members
    init_members();

    // Set parameter names
    if (coordsys == "CEL") {
        m_lon.name("RA");
        m_lat.name("DEC");
    }
    else {
        m_lon.name("GLON");
        m_lat.name("GLAT");
    }

    // Assign direction
    this->dir(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Value constructor
 *
 * @param[in] lon Longitude of point source (deg).
 * @param[in] lat Latitude of point source (deg).
 * @param[in] coordsys Coordinate system (either "CEL" or "GAL")
 *
 * @exception GException::invalid_argument
 *            Invalid @p coordsys argument specified.
 *
 * Construct a point source spatial model from the longitude and latitude
 * of a point source. The @p coordsys parameter specifies the coordinate
 * system in which @p lon and @p lat are provided. By default a celestial
 * coordinate system is assumed, which means that @p lon and @p lat are
 * taken as Right Ascension and Declination. If "GAL" is specified, @p lon
 * and @p lat are taken as Galactic longitude and latitude.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const double&      lon,
                                                   const double&      lat,
                                                   const std::string& coordsys) :
                          GModelSpatial()
{
    // Throw an exception if the coordinate system is invalid
    if ((coordsys != "CEL") && (coordsys != "GAL")) {
        std::string msg = "Invalid coordinate system \""+coordsys+"\" "
                          "specified. Please specify either \"CEL\" or "
                          "\"GAL\".";
        throw GException::invalid_argument(G_CONSTRUCTOR2, msg);
    }

    // Initialise members
    init_members();

    // Set parameter names
    if (coordsys == "CEL") {
        m_lon.name("RA");
        m_lat.name("DEC");
    }
    else {
        m_lon.name("GLON");
        m_lat.name("GLAT");
    }

    // Set values
    m_lon.value(lon);
    m_lat.value(lat);

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
 * Read the point source information from an XML element with the following
 * format
 *
 *     <spatialModel type="PointSource">
 *       <parameter name="RA"  scale="1" value="83.6331" min="-360" max="360" free="0" />
 *       <parameter name="DEC" scale="1" value="22.0145" min="-90"  max="90"  free="0" />
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="PointSource">
 *       <parameter name="GLON" scale="1" value="184.5575" min="-360" max="360" free="0" />
 *       <parameter name="GLAT" scale="1" value="-5.7843"  min="-90"  max="90"  free="0" />
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialPointSource::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 2);

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
 * Depending on the coordinate system, write the point source information
 * into an XML element with the following format
 *
 *     <spatialModel type="PointSource">
 *       <parameter name="RA"  scale="1" value="83.6331" min="-360" max="360" free="0" />
 *       <parameter name="DEC" scale="1" value="22.0145" min="-90"  max="90"  free="0" />
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="PointSource">
 *       <parameter name="GLON" scale="1" value="184.5575" min="-360" max="360" free="0" />
 *       <parameter name="GLAT" scale="1" value="-5.7843"  min="-90"  max="90"  free="0" />
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialPointSource::write(GXmlElement& xml) const
{
    // Verify model type
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
 * @brief Returns model flux integrated in sky region
 *
 * @param[in] region Sky region.
 * @param[in] srcEng Energy.
 * @param[in] srcTime Time.
 * @return Flux (adimensional or ph/cm2/s).
 *
 * Returns point source flux within a sky region. If the point source is
 * contained within the sky region the flux will be 1, otherwise the flux
 * will be 0.
 ***************************************************************************/
double GModelSpatialPointSource::flux(const GSkyRegion& region,
                                      const GEnergy&    srcEng,
                                      const GTime&      srcTime) const
{
    // Initialise flux
    double flux = 0.0;

    // If point source region overlaps with sky region then set flux to 1
    if (this->region()->overlaps(region)) {
        flux = 1.0;
    }

    // Return flux
    return flux;
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
const GSkyDir& GModelSpatialPointSource::dir(void) const
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
 * @brief Set position of point source
 *
 * @param[in] dir Sky direction of point source.
 *
 * Sets the sky direction of the point source.
 ***************************************************************************/
void GModelSpatialPointSource::dir(const GSkyDir& dir)
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
 * @param[in] model Point source spatial model.
 ***************************************************************************/
void GModelSpatialPointSource::copy_members(const GModelSpatialPointSource& model)
{
    // Copy members
    m_type = model.m_type;   // Needed to conserve model type
    m_lon  = model.m_lon;
    m_lat  = model.m_lat;

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
    // Set sky region circle
    GSkyRegionCircle region(dir(), 0.0);

    // Set region (circumvent const correctness)
    const_cast<GModelSpatialPointSource*>(this)->m_region = region;

    // Return
    return;
}
