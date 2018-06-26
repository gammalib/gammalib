/***************************************************************************
 *           GModelSpatialDiffuseMap.cpp - Spatial map model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseMap.cpp
 * @brief Spatial map model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GHealpix.hpp"
#include "GModelSpatialDiffuseMap.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDiffuseMap g_spatial_map_seed;
const GModelSpatialRegistry   g_spatial_map_registry(&g_spatial_map_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialDiffuseMap g_spatial_map_legacy_seed(true, "SpatialMap");
const GModelSpatialRegistry   g_spatial_map_legacy_registry(&g_spatial_map_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_MC           "GModelSpatialDiffuseMap::mc(GEnergy&, GTime&, GRan&)"
#define G_READ                  "GModelSpatialDiffuseMap::read(GXmlElement&)"
#define G_WRITE                "GModelSpatialDiffuseMap::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_MC                                     //!< Debug MC method
//#define G_DEBUG_MC_CACHE                               //!< Debug MC cache


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty spatial map model.
 ***************************************************************************/
GModelSpatialDiffuseMap::GModelSpatialDiffuseMap(void) :
                         GModelSpatialDiffuse()
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
 * Constructs empty spatial map model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialDiffuseMap::GModelSpatialDiffuseMap(const bool&        dummy,
                                                 const std::string& type) :
                         GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs spatial map model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialDiffuseMap::GModelSpatialDiffuseMap(const GXmlElement& xml) :
                         GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File name constructor
 *
 * @param[in] filename File name.
 * @param[in] value Normalization factor.
 * @param[in] normalize Normalize map?
 *
 * Constructs spatial map model by loading a sky map from the file specified
 * by @p filename and by setting the @p value by which the map will be
 * multiplied (or normalized).
 ***************************************************************************/
GModelSpatialDiffuseMap::GModelSpatialDiffuseMap(const GFilename& filename,
                                                 const double&    value,
                                                 const bool&      normalize) :
                         GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set normalisation parameter
    m_value.value(value);

    // Set normalization flag
    m_normalize = normalize;

    // Load sky map
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky map constructor
 *
 * @param[in] map Sky map.
 * @param[in] value Normalization factor.
 * @param[in] normalize Normalize map.
 *
 * Constructs spatial map model by setting the sky @p map and by setting the
 * @p value by which the map will be multiplied (or normalized).
 ***************************************************************************/
GModelSpatialDiffuseMap::GModelSpatialDiffuseMap(const GSkyMap& map,
                                                 const double&  value,
                                                 const bool&    normalize) :
                         GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set normalization parameter
    m_value.value(value);

    // Set normalization flag
    m_normalize = normalize;

    // Set sky map
    m_map = map;

    // Prepare sky map
    prepare_map();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spatial map model.
 *
 * Construct a spatial map model by copying another spatial map model.
 ***************************************************************************/
GModelSpatialDiffuseMap::GModelSpatialDiffuseMap(const GModelSpatialDiffuseMap& model) :
                         GModelSpatialDiffuse(model)
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
GModelSpatialDiffuseMap::~GModelSpatialDiffuseMap(void)
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
 * @param[in] model Spatial map model.
 * @return Spatial map model.
 *
 * Assigns a spatial map model to another spatial map model.
 ***************************************************************************/
GModelSpatialDiffuseMap& GModelSpatialDiffuseMap::operator=(const GModelSpatialDiffuseMap& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialDiffuse::operator=(model);

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
 * @brief Clear spatial map model
 ***************************************************************************/
void GModelSpatialDiffuseMap::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatialDiffuse::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialDiffuse::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone spatial map model
 *
 * @return Pointer to deep copy of spatial map model.
 ***************************************************************************/
GModelSpatialDiffuseMap* GModelSpatialDiffuseMap::clone(void) const
{
    // Clone diffuse map
    return new GModelSpatialDiffuseMap(*this);
}


/***********************************************************************//**
 * @brief Return intensity of sky map
 *
 * @param[in] photon Incident photon.
 * @param[in] gradients Compute gradients?
 * @return Sky map intensity (\f$\mbox{ph cm}^{-2}\mbox{sr}^{-1}\mbox{s}^{-1}\f$)
 *
 * Returns the intensity of the sky map at the specified sky direction
 * multiplied by the normalization factor. If the sky direction falls outside
 * the sky map, an intensity of 0 is returned.
 *
 * If the @p gradients flag is true the method will also evaluate the partial
 * derivatives of the model.
 ***************************************************************************/
double GModelSpatialDiffuseMap::eval(const GPhoton& photon,
                                     const bool&    gradients) const
{
    // Get skymap intensity
    double intensity = m_map(photon.dir());

    // Optionally compute partial derivatives
    if (gradients) {

        // Compute partial derivatives of the parameter values
        double g_value = (m_value.is_free()) ? intensity * m_value.scale() : 0.0;

        // Set gradient
        m_value.factor_gradient(g_value);

    } // endif: computed partial derivatives

    // Return intensity times normalization factor
    return (intensity * m_value.value());
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy (ignored).
 * @param[in] time Photon arrival time (ignored).
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * @exception GException::invalid_value
 *            No sky map defined.
 * @exception GException::invalid_return_value
 *            Simulation cone not defined, does not overlap with map or map
 *            is empty.
 *
 * Returns a random sky direction according to the intensity distribution of
 * the model sky map. The method uses a rejection method to determine the sky
 * direction. If no sky direction could be determined, the method throws an
 * GException::invalid_return_value exception.
 ***************************************************************************/
GSkyDir GModelSpatialDiffuseMap::mc(const GEnergy& energy,
                                    const GTime&   time,
                                    GRan&          ran) const
{
    // Throw an exception if the maximum MC intensity is not positive. This
    // can be the case because the simulation cone has not been defined or
    // because it does not overlap with the sky map
    if (m_mc_max <= 0.0) {
        std::string msg = "Simulation cone has not been defined or does not "
                          "overlap with the sky map. Please specify a valid "
                          "simulation cone.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Determine number of skymap pixels
    int npix = m_map.npix();

    // Throw an exception if there are no sky map pixels
    if (npix <= 0) {
        std::string msg = "No sky map defined. Please specify a valid sky map.";
        throw GException::invalid_value(G_MC, msg);
    }

    // Allocate sky direction
    GSkyDir dir;

    // Debug option: initialise counter
    #if defined(G_DEBUG_MC)
    int num_iterations = 0;
    #endif

    // Get sky direction
    while (true) {

        // Debug option: increment counter
        #if defined(G_DEBUG_MC)
        num_iterations++;
        #endif

        // Simulate random sky direction within Monte Carlo simulation cone
        double theta = std::acos(1.0 - ran.uniform() * m_mc_one_minus_cosrad) *
                       gammalib::rad2deg;
        double phi   = 360.0 * ran.uniform();
        dir = m_mc_centre;
        dir.rotate_deg(phi, theta);

        // Get map value at simulated sky direction. If the map value is non-
        // positive then simulate a new sky direction.
        double value = m_map(dir);
        if (value <= 0.0) {
            continue;
        }

        // Get uniform random number
        double uniform = ran.uniform() * m_mc_max;

        // Exit loop if the random number is not larger than the sky map value
        if (uniform <= value) {
            break;
        }

    } // endwhile: loop until sky direction was accepted

    // Debug option: log counter
    #if defined(G_DEBUG_MC)
    std::cout << num_iterations << " ";
    #endif

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Return normalization of diffuse map for Monte Carlo simulations
 *
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (degrees).
 * @return Normalization.
 *
 * Returns the normalization of a diffuse map. The normalization is given
 * by the model value times the integrated flux in the sky map over a
 * circular region.
 ***************************************************************************/
double GModelSpatialDiffuseMap::mc_norm(const GSkyDir& dir,
                                        const double&  radius) const
{
    // Set the MC cone
    set_mc_cone(dir, radius);

    // Retrieve normalization
    double norm = m_mc_norm * value();

    // Return normalization
    return norm;
}


/***********************************************************************//**
 * @brief Signals whether model contains sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction falls within the bounding circle of
 * the diffuse map.
 ***************************************************************************/
bool GModelSpatialDiffuseMap::contains(const GSkyDir& dir,
                                       const double&  margin) const
{
    // Initialise containment flag
    bool contains = false;

    // Continue only if radius is positive
    if (m_radius > 0.0) {

        // Compute distance to centre
        double distance = m_centre.dist_deg(dir);

        // If distance is smaller than radius plus margin we consider
        // the position to be contained within the bounding circle
        if (distance < m_radius + margin) {
            contains = true;
        }

    }

    // Return containment
    return (contains);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Model parameters not found in XML element.
 *
 * Reads the spatial information for a diffuse map from an XML element. The
 * expected format of the XML element is
 *
 *     <spatialModel type="DiffuseMap" file="myfile.fits" normalize="1">
 *       <parameter name="Normalization" ../>
 *     </spatialModel>
 *
 * The @p file attribute provides the filename of the diffuse map FITS file.
 * The filename may be either an absolute filename (starting with '/') or a
 * relative filename. If no access path is given, the file is expected to
 * reside in the same location as the XML file.
 *
 * The @p normalize attribute specifies whether the sky map should be
 * normalized to unity flux or not. If the attribute is not given, the map
 * will be automatically normalized. To prevent normalization,
 * @p normalize="0" needs to be specified.
 ***************************************************************************/
void GModelSpatialDiffuseMap::read(const GXmlElement& xml)
{
    // If "Normalization" parameter exists then read parameter from this
    // XML element
    if (gammalib::xml_has_par(xml, "Normalization")) {
        const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Normalization");
        m_value.read(*par);
    }

    // ... otherwise try reading parameter from "Prefactor" parameter
    #if defined(G_LEGACY_XML_FORMAT)
    else if (gammalib::xml_has_par(xml, "Prefactor")) {
        const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Prefactor");
        m_value.read(*par);
    }
    #endif

    // ... otherwise throw an exception
    else {
        #if defined(G_LEGACY_XML_FORMAT)
        std::string msg = "Diffuse map model requires either "
                          "\"Normalization\" or \"Prefactor\" parameter.";
        #else
        std::string msg = "Diffuse map model requires \"Normalization\" "
                          "parameter.";
        #endif
        throw GException::invalid_value(G_READ, msg);
    }

    // Get optional normalization attribute
    m_normalize     = true;
    m_has_normalize = false;
    if (xml.has_attribute("normalize")) {
        m_has_normalize = true;
        std::string arg = xml.attribute("normalize");
        if (arg == "0" || gammalib::tolower(arg) == "false") {
            m_normalize = false;
        }
    }

    // Load sky map.
    load(gammalib::xml_file_expand(xml, xml.attribute("file")));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type "SpatialMap"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spatial information for a diffuse map into an XML element. The
 * format of the XML element is
 *
 *     <spatialModel type="DiffuseMap" file="myfile.fits" normalize="1">
 *       <parameter name="Prefactor" value="1" min="0.1" max="10" scale="1" free="0"/>
 *     </spatialModel>
 *
 * The @p file attribute provides the filename of the diffuse map FITS file.
 * The filename may be either an absolute filename (starting with '/') or a
 * relative filename. If no access path is given, the file is expected to
 * reside in the same location as the XML file.
 *
 * The @p normalize attribute specifies whether the sky map should be
 * normalized to unity flux or not. The attribute will only be written if the
 * normalization is disabled.
 ***************************************************************************/
void GModelSpatialDiffuseMap::write(GXmlElement& xml) const
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

    // Set sky map file name
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // If XML element has 0 nodes then append parameter node. The name
    // of the node is "Prefactor" as this is the Fermi/LAT standard.
    // We thus assure that the XML files will be compatible with
    // Fermi/LAT.
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Prefactor\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spatial map model requires exactly 1 parameter.");
    }

    // Get pointer on model parameter
    GXmlElement* par = xml.element("parameter", 0);

    // Set or update parameter
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Prefactor") {
        m_value.write(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Spatial map model requires either \"Prefactor\" or"
              " \"Normalization\" parameter.");
    }

    // Set optional normalization attribute
    if (m_has_normalize || !m_normalize) {
        if (m_normalize) {
            xml.attribute("normalize", "1");
        }
        else {
            xml.attribute("normalize", "0");
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set Monte Carlo simulation cone
 *
 * @param[in] centre Simulation cone centre.
 * @param[in] radius Simulation cone radius (degrees).
 *
 * Sets the simulation cone centre and radius that defines the directions
 * that will be simulated using the mc() method and pre-computes the maximum
 * intensity and the spatially integrated flux of the map within the
 * simulation cone region.
 ***************************************************************************/
void GModelSpatialDiffuseMap::set_mc_cone(const GSkyDir& centre,
                                          const double&  radius) const
{
    // Continue only if the simulation cone has changed
    if ((centre != m_mc_centre) || (radius != m_mc_radius)) {

        // Save simulation cone definition
        m_mc_centre = centre;
        m_mc_radius = radius;

        // Pre-compute 1 - cosine of radius
        m_mc_one_minus_cosrad = 1.0 - std::cos(m_mc_radius*gammalib::deg2rad);

        // Initialise map maximum and normalisation
        m_mc_max  = 0.0;
        m_mc_norm = 0.0;

        // Determine number of skymap pixels
        int npix = m_map.npix();

        // Continue only if there are pixels
        if (npix > 0) {

            // Compute flux and maximum map intensity within the simulation cone
            double sum     = 0.0;
            double sum_map = 0.0;
            for (int i = 0; i < npix; ++i) {

                // Get map flux, intensity and distance from MC cone centre
                double flux      = m_map.flux(i);
                double intensity = m_map(i);
                double distance  = centre.dist_deg(m_map.pix2dir(i));

                // Add flux if positive
                if (flux > 0.0) {
                    if (distance <= radius) {
                        sum += flux; // flux within simulation cone
                    }
                    sum_map += flux; // total flux
                }
    
                // Update maximum intensity
                if (distance <= radius) {
                    if (intensity > m_mc_max) {
                        m_mc_max = intensity;
                    }
                }

            } // endfor: looped over all map pixels

            // Set the normalization factor for the MC simulations. In case
            // that the map is normalised, this is the fraction of the flux
            // that is comprised within the simulation cone. For non-normalised
            // maps, this is simply the flux comprised within the simulation
            // cone.
            if (sum_map > 0.0) {
                if (m_normalize) {
                    m_mc_norm = sum / sum_map;
                }
                else {
                    m_mc_norm = sum;
                }
            }

            // Log maximum intensity and total flux for debugging
            #if defined(G_DEBUG_MC_CACHE)
            std::cout << "GModelSpatialDiffuseMap::set_mc_cone:" << std::endl;
            std::cout << "  Maximum map intensity:";
            std::cout << m_mc_max << " ph/cm2/s/sr" << std::endl;
            std::cout << "  Spatially integrated flux:" << std::endl;
            std::cout << sum << " ph/cm2/s" << std::endl;
            std::cout << "  Map normalisation:" << std::endl;
            std::cout << m_mc_norm << " ph/cm2/s" << std::endl;
            #endif

        } // endif: there were map pixels

    } // endif: simulation cone has changed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print map information
 *
 * @param[in] chatter Chattiness.
 * @return String with diffuse map model information.
 ***************************************************************************/
std::string GModelSpatialDiffuseMap::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialDiffuseMap ===");

        // Append parameters
        result.append("\n"+gammalib::parformat("Sky map file"));
        result.append(m_filename);
        result.append("\n"+gammalib::parformat("Map normalization"));
        result.append(gammalib::str(m_mc_norm)+" ph/cm2/s");
        if (normalize()) {
            result.append(" [normalized]");
        }
        result.append("\n"+gammalib::parformat("Map centre"));
        result.append(m_centre.print());
        result.append("\n"+gammalib::parformat("Map radius"));
        result.append(gammalib::str(m_radius)+" deg");
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
 * @brief Load skymap into the model class
 *
 * @param[in] filename Sky map file.
 *
 * Loads skymap into the model class. The method calls the protected method
 * prepare_map() that prepares the map for usage by the class.
 *
 * If the @p filename is empty, no map will be loaded.
 ***************************************************************************/
void GModelSpatialDiffuseMap::load(const GFilename& filename)
{
    // Initialise skymap
    m_map.clear();

    // Continue only if filename is not empty
    if (!filename.is_empty()) {

        // Store filename of skymap (for XML writing). Note that we do not
        // expand any environment variable at this level, so that if we
        // write back the XML element we write the filepath with the
        // environment variables
        m_filename = filename;

        // Load sky map
        m_map.load(m_filename);

        // Prepare sky map
        prepare_map();

    } // endif: filename was not empty

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
void GModelSpatialDiffuseMap::init_members(void)
{
    // Initialise model type
    m_type = "DiffuseMap";

    // Initialise Value
    m_value.clear();
    m_value.name("Normalization");
    m_value.value(1.0);
    m_value.scale(1.0);
    m_value.range(0.001, 1000.0);
    m_value.gradient(0.0);
    m_value.fix();
    m_value.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Initialise other members
    m_map.clear();
    m_filename.clear();
    m_normalize     = true;
    m_has_normalize = false;
    m_centre.clear();
    m_radius        = 0.0;
    m_region.clear();

    // Initialise MC cache
    m_mc_centre.clear();
    m_mc_radius           = -1.0;   //!< Signal initialisation
    m_mc_one_minus_cosrad =  1.0;
    m_mc_norm             =  0.0;
    m_mc_max              =  0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spatial map cube model.
 ***************************************************************************/
void GModelSpatialDiffuseMap::copy_members(const GModelSpatialDiffuseMap& model)
{
    // Copy members
    m_type          = model.m_type;
    m_value         = model.m_value;
    m_map           = model.m_map;
    m_filename      = model.m_filename;
    m_normalize     = model.m_normalize;
    m_has_normalize = model.m_has_normalize;
    m_centre        = model.m_centre;
    m_radius        = model.m_radius;
    m_region        = model.m_region;

    // Copy MC cache
    m_mc_centre           = model.m_mc_centre;
    m_mc_radius           = model.m_mc_radius;
    m_mc_one_minus_cosrad = model.m_mc_one_minus_cosrad;
    m_mc_norm             = model.m_mc_norm;
    m_mc_max              = model.m_mc_max;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialDiffuseMap::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Prepare sky map after loading
 *
 * Prepares a sky map after loading. Negative, infinite or undefined skymap
 * pixels are set to zero intensity. The method also determine the centre
 * and radius of a circle enclosing the map, and the Monte Carlo simulation
 * cone is set to this circle.
 *
 * If normalize() is true, the map is furthermore normalised so that the
 * total flux in the map amounts to 1 ph/cm2/s. Negative skymap pixels are
 * set to zero intensity.
 ***************************************************************************/
void GModelSpatialDiffuseMap::prepare_map(void)
{
    // Initialise centre and radius
    m_centre.clear();
    m_radius = 0.0;

    // Determine number of skymap pixels
    int npix = m_map.npix();

    // Continue only if there are skymap pixels
    if (npix > 0) {

        // Initialise flux sum
        double sum = 0.0;

        // Compute flux sum and set negative or invalid pixels to zero
        // intensity.
        for (int i = 0; i < npix; ++i) {
            double flux = m_map.flux(i);
            if (flux < 0.0 ||
                gammalib::is_notanumber(flux) ||
                gammalib::is_infinite(flux)) {
                m_map(i) = 0.0;
                flux     = 0.0;
            }
            sum += flux;
        }

        // Optionally normalize the sky map.
        if (sum > 0.0) {
            if (normalize()) {
                for (int i = 0; i < npix; ++i) {
                    m_map(i) /= sum;
                }
            }
        }

        // If we have a HealPix map then set radius to 180 deg
        if (m_map.projection()->code() == "HPX") {
            m_radius = 180.0;
        }

        // ... otherwise compute map centre and radius
        else {

            // Get map centre
            GSkyPixel centre(m_map.nx()/2.0, m_map.ny()/2.0);
            m_centre = m_map.pix2dir(centre);

            // Determine map radius
            for (int i = 0; i < npix; ++i) {
                double radius = m_map.inx2dir(i).dist_deg(m_centre);
                if (radius > m_radius) {
                    m_radius = radius;
                }
            }

        } // endelse: computed map centre and radius

        // Set simulation cone
        set_mc_cone(m_centre, m_radius);

    } // endif: there were skymap pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialDiffuseMap::set_region(void) const
{
    // Set sky region centre to bounding circle centre
    m_region.centre(m_centre);

    // Set sky region radius to bounding circle radius
    m_region.radius(m_radius);

    // Return
    return;
}
