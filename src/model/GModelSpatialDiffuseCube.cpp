/***************************************************************************
 *       GModelSpatialDiffuseCube.cpp - Spatial map cube model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseCube.cpp
 * @brief Spatial map cube model class implementation
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
#include "GEnergies.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"
#include "GXmlElement.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpatialDiffuseCube.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDiffuseCube g_spatial_cube_seed;
const GModelSpatialRegistry    g_spatial_cube_registry(&g_spatial_cube_seed);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpatialDiffuseCube g_spatial_cube_legacy_seed(true, "MapCubeFunction");
const GModelSpatialRegistry    g_spatial_cube_legacy_registry(&g_spatial_cube_legacy_seed);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_MC          "GModelSpatialDiffuseCube::mc(GEnergy&, GTime&, GRan&)"
#define G_READ                 "GModelSpatialDiffuseCube::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialDiffuseCube::write(GXmlElement&)"
#define G_ENERGIES           "GModelSpatialDiffuseCube::energies(GEnergies&)"
#define G_LOAD_CUBE         "GModelSpatialDiffuseCube::load_cube(GFilename&)"

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
 * Constructs empty map cube model.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(void) :
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
 * Constructs empty map cube model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const bool&        dummy,
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
 * Constructs map cube model by extracting information from an XML element.
 * See the read() method for more information about the expected structure
 * of the XML element.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const GXmlElement& xml) :
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
 * @brief Filename constructor
 *
 * @param[in] filename File name.
 * @param[in] value Normalization factor (defaults to 1).
 *
 * Constructs map cube model by loading a map cube from @p filename and by
 * assigning the normalization @p value.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const GFilename& filename,
                                                   const double&    value) :
                          GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set parameter
    m_value.value(value);

    // Perform autoscaling of parameter
    autoscale();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky map constructor
 *
 * @param[in] cube Sky map cube.
 * @param[in] energies Sky map energies.
 * @param[in] value Normalization factor (defaults to 1).
 *
 * Constructs map cube model by extracting a @p cube from a sky map. The
 * constructor also assigns the energy values for all maps and sets the
 * scaling @p value. The filename will remain blank.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const GSkyMap&   cube,
                                                   const GEnergies& energies,
                                                   const double&    value) :
                          GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set parameter
    m_value.value(value);

    // Perform autoscaling of parameter
    autoscale();

    // Set map cube
    this->cube(cube);

    // Set energies
    this->energies(energies);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Map cube model.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const GModelSpatialDiffuseCube& model) :
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
GModelSpatialDiffuseCube::~GModelSpatialDiffuseCube(void)
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
 * @param[in] model Map cube model.
 * @return Map cube model.
 ***************************************************************************/
GModelSpatialDiffuseCube& GModelSpatialDiffuseCube::operator=(const GModelSpatialDiffuseCube& model)
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
 * @brief Clear map cube model
 ***************************************************************************/
void GModelSpatialDiffuseCube::clear(void)
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
 * @brief Clone map cube model
 *
 * @return Pointer to deep copy of map cube model.
 ***************************************************************************/
GModelSpatialDiffuseCube* GModelSpatialDiffuseCube::clone(void) const
{
    // Clone map cube model
    return new GModelSpatialDiffuseCube(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] photon Incident photon.
 * @param[in] gradients Compute gradients?
 * @return Sky map intensity (\f$\mbox{ph cm}^{-2}\mbox{sr}^{-1}\mbox{s}^{-1}\f$)
 *
 * Computes the spatial diffuse model as function of photon parameters.
 *
 * If the @p gradients flag is true the method will also evaluate the partial
 * derivatives of the model.
 ***************************************************************************/
double GModelSpatialDiffuseCube::eval(const GPhoton& photon,
                                      const bool&    gradients) const
{
    // Get log-log interpolated cube intensity
    double intensity = cube_intensity(photon);

    // Set the intensity times the scaling factor as model value
    double value = intensity * m_value.value();

    // Optionally compute partial derivatives
    if (gradients) {

        // Compute partial derivatives of the parameter value. In case that
        // the value is negative set the gradient to zero.
        double g_value = (m_value.is_free()) ? intensity * m_value.scale() : 0.0;
        if (value < 0.0) {
            g_value = 0.0;
        }

        // Set gradient
        m_value.factor_gradient(g_value);

    } // endif: computed partial derivatives

    // Make sure that value is not negative
    if (value < 0.0) {
        value = 0.0;
    }

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
 * @exception GException::invalid_value
 *            No map cube defined.
 *            No energy boundaries defined.
 * @exception GException::invalid_return_value
 *            Simulation cone not defined, does not overlap with map cube or
 *            map cube is empty for the specified energy.
 *
 * Returns a random sky direction according to the intensity distribution of
 * the model sky map and the specified energy. The method uses a rejection
 * method to determine the sky direction. If no sky direction could be
 * determined, the method throws an GException::invalid_return_value
 * exception.
 ***************************************************************************/
GSkyDir GModelSpatialDiffuseCube::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Fetch cube
    fetch_cube();

    // Determine number of skymap pixels
    int npix = pixels();

    // Throw an exception if there are no sky map pixels
    if (npix <= 0) {
        std::string msg = "No map cube defined. Please specify a valid map cube.";
        throw GException::invalid_value(G_MC, msg);
    }

    // Throw an exception if no energy boundaries are defined
    if (m_ebounds.size() < 1) {
        std::string msg = "The energy boundaries of the maps in the cube have "
                          "not been defined. Maybe the map cube file is missing "
                          "the \"ENERGIES\" extension which defines the energy "
                          "of each map in the cube. Please provide the energy "
                          "information.";
        throw GException::invalid_value(G_MC, msg);
    }

    // Set energy for interpolation
    m_logE.set_value(energy.log10MeV());

    // Compute maximum map value within simulation cone by log-log interpolation
    // of the maximum values that are stored for each map
    double max       = 0.0;
    double max_left  = m_mc_max[m_logE.inx_left()];
    double max_right = m_mc_max[m_logE.inx_right()];
    if (max_left > 0.0 && max_right > 0.0) {
        double max_log = m_logE.wgt_left()  * std::log(max_left) +
                         m_logE.wgt_right() * std::log(max_right);
        max = std::exp(max_log);
    }
    else if (max_left > 0.0) {
        max = max_left;
    }
    else if (max_right > 0.0) {
        max = max_right;
    }

    // Throw an exception if the maximum MC intensity is not positive. This
    // can happen if the simulation cone has not been defined or if there is
    // no overlap with the sky map or if the sky map is empty
    if (max <= 0.0) {
        std::string msg;
        if (m_mc_max.size() == 0) {
            msg = "Simulation cone has not been defined. Please specify a "
                  "valid simulation cone before calling the method.";
        }
        else {
            msg = "The map cube is empty at "+energy.print()+" within the "
                  "simulation cone. Please specify a valid map cube or "
                  "restrict the simulated energies to values where the map "
                  "cube is non-zero within the simulation cone.";
        }
        throw GException::invalid_return_value(G_MC, msg);
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

        // Get map value at simulated sky direction
        double value       = 0.0;
        double value_left  = m_cube(dir, m_logE.inx_left());
        double value_right = m_cube(dir, m_logE.inx_right());
        if (value_left > 0.0 && value_right > 0.0) {
            double value_log = m_logE.wgt_left()  * std::log(value_left) +
                               m_logE.wgt_right() * std::log(value_right);
            value = std::exp(value_log);
        }
        else if (value_left > 0.0) {
            value = value_left;
        }
        else if (value_right > 0.0) {
            value = value_right;
        }

        // If map value is non-positive then simulate a new sky direction
        if (value <= 0.0) {
            continue;
        }

        // Get uniform random number
        double uniform = ran.uniform() * max;

        // Exit loop if the random number is not larger than the map cube value
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
 * @brief Signals whether model contains sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction falls within the bounding circle of
 * the diffuse cube.
 *
 * @todo To be implemented.
 ***************************************************************************/
bool GModelSpatialDiffuseCube::contains(const GSkyDir& dir,
                                        const double&  margin) const
{
    return (true);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Model parameters not found in XML element.
 *
 * Read the map cube information from an XML element. The XML element should
 * have the format
 *
 *     <spatialModel type="DiffuseMapCube" file="test_file.fits">
 *       <parameter name="Normalization" ../>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialDiffuseCube::read(const GXmlElement& xml)
{
    // If "Normalization" parameter exists then read parameter from this
    // XML element
    if (gammalib::xml_has_par(xml, "Normalization")) {
        const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Normalization");
        m_value.read(*par);
    }

    // ... otherwise try reading parameter from "Value" parameter
    #if defined(G_LEGACY_XML_FORMAT)
    else if (gammalib::xml_has_par(xml, "Value")) {
        const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Value");
        m_value.read(*par);
    }
    #endif

    // ... otherwise throw an exception
    else {
        #if defined(G_LEGACY_XML_FORMAT)
        std::string msg = "Diffuse map cube model requires either "
                          "\"Normalization\" or \"Value\" parameter.";
        #else
        std::string msg = "Diffuse map cube model requires \"Normalization\" "
                          "parameter.";
        #endif
        throw GException::invalid_value(G_READ, msg);
    }

    // Save filename
    m_filename = gammalib::xml_file_expand(xml, xml.attribute("file"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type "MapCubeFunction"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the map cube information into an XML element. The XML element will
 * have either the format
 *
 *     <spatialModel type="DiffuseMapCube" file="test_file.fits">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *     </spatialModel>
 *
 * or alternatively
 *
 *     <spatialModel type="DiffuseMapCube" file="test_file.fits">
 *       <parameter name="Value" scale="1" value="1" min="0.1" max="10" free="0"/>
 *     </spatialModel>
 *
 * The latter format is the default for newly written XML elements. 
 ***************************************************************************/
void GModelSpatialDiffuseCube::write(GXmlElement& xml) const
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

    // If XML element has 0 nodes then append parameter node. The name
    // of the node is "Normalization" as this is the Fermi-LAT standard.
    // We thus assure that the XML files will be compatible with
    // Fermi-LAT.
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Normalization\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Map cube spatial model requires exactly 1 parameter.");
    }

    // Get pointers on model parameter
    GXmlElement* par = xml.element("parameter", 0);

    // Set or update parameter
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value") {
        m_value.write(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Map cube spatial model requires either \"Value\" or"
              " \"Normalization\" parameter.");
    }

    // Set filename
    //xml.attribute("file", m_filename);
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set map cube
 *
 * @param[in] cube Sky map.
 *
 * Set the map cube of the spatial map cube model.
 ***************************************************************************/
void GModelSpatialDiffuseCube::cube(const GSkyMap& cube)
{
    // Clear filename and signal map as unloaded
    m_filename.clear();
    m_loaded = false;

    // Assign map
    m_cube = cube;

    // Update Monte Carlo cache
    update_mc_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energies for map cube
 *
 * @param[in] energies Sky map energies.
 *
 * @exception GException::invalid_argument
 *            Specified sky map energies incompatible with map cube.
 *
 * Sets the energies for the map cube.
 ***************************************************************************/
void GModelSpatialDiffuseCube::energies(const GEnergies& energies)
{
    // Initialise energies
    m_logE.clear();

    // Fetch cube
    fetch_cube();

    // Extract number of energies in vector
    int num = energies.size();

    // Check if energy binning is consistent with number of maps in the cube
    if (num != m_cube.nmaps() ) {
        std::string msg = "Number of specified energies ("+gammalib::str(num)+")"
                          " does not match the number of maps ("
                          ""+gammalib::str(m_cube.nmaps())+" in the map cube.\n"
                          "The energies argument shall provide a vector of length"
                          " "+gammalib::str(m_cube.nmaps())+".";
        throw GException::invalid_argument(G_ENERGIES, msg);
    }

    // Set log10(energy) nodes, where energy is in units of MeV
    for (int i = 0; i < num; ++i) {
        m_logE.append(energies[i].log10MeV());
    }

    // Set energy boundaries
    set_energy_boundaries();

    // Update MC cache
    update_mc_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns map cube energies
 *
 * @return Map cube energies.
 *
 * Returns the energies for the map cube in a vector.
 ***************************************************************************/
GEnergies GModelSpatialDiffuseCube::energies(void)
{
    // Initialise energies container
    GEnergies energies;

    // Fetch cube
    fetch_cube();

    // Get number of map energies
    int num = m_logE.size();

    // Continue only if there are maps in the cube
    if (num > 0) {

        // Reserve space for all energies
        energies.reserve(num);

        // Set log10(energy) nodes, where energy is in units of MeV
        for (int i = 0; i < num; ++i) {
            GEnergy energy;
            energy.log10MeV(m_logE[i]);
            energies.append(energy);
        }

    } // endif: there were maps in the cube

    // Return energies
    return energies;
}


/***********************************************************************//**
 * @brief Set Monte Carlo simulation cone
 *
 * @param[in] centre Simulation cone centre.
 * @param[in] radius Simulation cone radius (degrees).
 *
 * Sets the simulation cone centre and radius that defines the directions
 * that will be simulated using the mc() method and pre-computes the maximum
 * intensity and the spatially integrated flux of each map within the
 * simulation cone region.
 ***************************************************************************/
void GModelSpatialDiffuseCube::set_mc_cone(const GSkyDir& centre,
                                           const double&  radius) const
{
    // Continue only if the simulation cone has changed
    if ((centre != m_mc_centre) || (radius != m_mc_radius)) {

        // Save simulation cone definition
        m_mc_centre = centre;
        m_mc_radius = radius;

        // Pre-compute 1 - cosine of radius
        m_mc_one_minus_cosrad = 1.0 - std::cos(m_mc_radius*gammalib::deg2rad);

        // Initialise Monte Carlo cache
        m_mc_max.clear();
        m_mc_spectrum.clear();

        // Fetch cube
        fetch_cube();

        // Determine number of cube pixels and maps
        int npix  = pixels();
        int nmaps = maps();

        // Continue only if there are pixels and maps
        if (npix > 0 && nmaps > 0) {

            // Reserve space in cache
            m_mc_max.reserve(nmaps);
            m_mc_spectrum.reserve(nmaps);

            // Loop over all maps
            for (int i = 0; i < nmaps; ++i) {

                // Compute flux and maximum map intensity within the
                // simulation cone
                double total_flux    = 0.0;
                double max_intensity = 0.0;
                for (int k = 0; k < npix; ++k) {
                    double distance = centre.dist_deg(m_cube.pix2dir(k));
                    if (distance <= radius) {
                        double flux = m_cube.flux(k,i);
                        if (flux > 0.0) {
                            total_flux += flux;
                        }
                        double value = m_cube(k,i);
                        if (value > max_intensity) {
                            max_intensity = value;
                        }
                    }
                }

                // Store maximum map intensity
                m_mc_max.push_back(max_intensity);

                // Store flux as spectral node
                if (m_logE.size() == nmaps) {

                    // Set map energy
                    GEnergy energy;
                    energy.log10MeV(m_logE[i]);

                    // Only append node if the integrated flux is positive
                    // (as we do a log-log interpolation)
                    if (total_flux > 0.0) {
                        m_mc_spectrum.append(energy, total_flux);
                    }

                }

            } // endfor: looped over all maps

            // Log maximum intensity and total flux for debugging
            #if defined(G_DEBUG_MC_CACHE)
            std::cout << "GModelSpatialDiffuseCube::set_mc_cone:" << std::endl;
            std::cout << "  Maximum map intensity:" << std::endl;
            for (int i = 0; i < m_mc_max.size(); ++i) {
                GEnergy energy;
                energy.log10MeV(m_logE[i]);
                std::cout << "  " << i;
                std::cout << " " << energy;
                std::cout << " " << m_mc_max[i] << " ph/cm2/s/sr/MeV";
                std::cout << std::endl;
            }
            std::cout << "  Spatially integrated flux:" << std::endl;
            for (int i = 0; i < m_mc_spectrum.nodes(); ++i) {
                std::cout << "  " << i;
                std::cout << " " << m_mc_spectrum.energy(i);
                std::cout << " " << m_mc_spectrum.intensity(i);
                std::cout << " ph/cm2/s/MeV" << std::endl;
            }
            #endif

        } // endif: there were cube pixels and maps

    } // endif: simulation cone has changed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load cube into the model class
 *
 * @param[in] filename cube file.
 *
 * Loads cube into the model class.
 ***************************************************************************/
void GModelSpatialDiffuseCube::load(const GFilename& filename)
{
    // Load cube
    load_cube(filename);

    // Update Monte Carlo cache
    update_mc_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save cube into FITS file
 *
 * @param[in] filename Cube FITS file name.
 * @param[in] clobber Overwrite existing FITS file (default: false).
 *
 * Saves spatial cube model into a FITS file.
 ***************************************************************************/
void GModelSpatialDiffuseCube::save(const GFilename& filename,
                                    const bool&      clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write skymap into FITS file
 *
 * @param[in] file FITS file pointer.
 ***************************************************************************/
void GModelSpatialDiffuseCube::write(GFits& file) const
{
    // Write cube
    m_cube.write(file);

    // Create energy intervals
    GEnergies energies;
    for (int i = 0; i < m_logE.size(); ++i) {
        double energy = std::pow(10.0, m_logE[i]);
        energies.append(GEnergy(energy, "MeV"));
    }

    // Write energy intervals
    energies.write(file);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print map cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialDiffuseCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialDiffuseCube ===");

        // Append parameters
        result.append("\n"+gammalib::parformat("Map cube file"));
        result.append(m_filename);
        if (m_loaded) {
            result.append(" [loaded]");
        }
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

        // Append detailed information only if a map cube exists
        if (m_cube.npix() > 0) {

            // NORMAL: Append sky map
            if (chatter >= NORMAL) {
                result.append("\n"+m_cube.print(chatter));
            }

            // EXPLICIT: Append energy nodes
            if (chatter >= EXPLICIT && m_logE.size() > 0) {
                result.append("\n"+gammalib::parformat("Map energy values"));
                if (m_logE.size() > 0) {
                    for (int i = 0; i < m_logE.size(); ++i) {
                        result.append("\n"+gammalib::parformat("  Map "+gammalib::str(i+1)));
                        result.append(gammalib::str(std::pow(10.0, m_logE[i])));
                        result.append(" MeV (log10E=");
                        result.append(gammalib::str(m_logE[i]));
                        result.append(")");
                        if (m_ebounds.size() == m_logE.size()) {
                            result.append(" [");
                            result.append(m_ebounds.emin(i).print());
                            result.append(", ");
                            result.append(m_ebounds.emax(i).print());
                            result.append("]");
                        }
                    }
                }
                else {
                    result.append("not specified");
                }
            }

            // VERBOSE: Append MC cache
            if (chatter >= VERBOSE) {
                result.append("\n"+gammalib::parformat("Map flux"));
                if (m_mc_spectrum.nodes() > 0) {
                    for (int i = 0; i < m_mc_spectrum.nodes(); ++i) {
                        result.append("\n"+gammalib::parformat("  Map "+gammalib::str(i+1)));
                        result.append(gammalib::str(m_mc_spectrum.intensity(i)));
                    }
                }
                else {
                    result.append("not specified");
                }
            }

        } // endif: map cube exists

    } // endif: chatter was not silent

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
void GModelSpatialDiffuseCube::init_members(void)
{
    // Initialise model type
    m_type = "DiffuseMapCube";

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
    m_filename.clear();
    m_cube.clear();
    m_logE.clear();
    m_ebounds.clear();
    m_loaded = false;
    m_region.clear();

    // Initialise MC cache
    m_mc_centre.clear();
    m_mc_radius           = -1.0;    //!< Signal that initialisation is needed
    m_mc_one_minus_cosrad =  1.0;
    m_mc_max.clear();
    m_mc_spectrum.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spatial map cube model.
 ***************************************************************************/
void GModelSpatialDiffuseCube::copy_members(const GModelSpatialDiffuseCube& model)
{
    // Copy members
    m_type     = model.m_type;
    m_value    = model.m_value;
    m_filename = model.m_filename;
    m_cube     = model.m_cube;
    m_logE     = model.m_logE;
    m_ebounds  = model.m_ebounds;
    m_loaded   = model.m_loaded;
    m_region   = model.m_region;

    // Copy MC cache
    m_mc_centre           = model.m_mc_centre;
    m_mc_radius           = model.m_mc_radius;
    m_mc_one_minus_cosrad = model.m_mc_one_minus_cosrad;
    m_mc_max              = model.m_mc_max;
    m_mc_spectrum         = model.m_mc_spectrum;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialDiffuseCube::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch cube
 *
 * Load diffuse cube if it is not yet loaded. The loading is thread save.
 ***************************************************************************/
void GModelSpatialDiffuseCube::fetch_cube(void) const
{
    // Load cube if it is not yet loaded
    if (!m_loaded && !m_filename.is_empty()) {

        // Put in a OMP critical zone
        #pragma omp critical(GModelSpatialDiffuseCube_fetch_cube)
        {
            const_cast<GModelSpatialDiffuseCube*>(this)->load_cube(m_filename);
        } // end of pragma

    } // endif: file not

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load cube
 *
 * @param[in] filename Cube file.
 *
 * @exception GException::invalid_value
 *            Number of maps in cube mismatches number of energy bins.
 *
 * Load diffuse cube.
 ***************************************************************************/
void GModelSpatialDiffuseCube::load_cube(const GFilename& filename)
{
    // Initialise skymap
    m_cube.clear();
    m_logE.clear();

    // Store filename of cube (for XML writing). Note that we do
    // not expand any environment variable at this level, so that
    // if we write back the XML element we write the filepath with
    // the environment variables
    m_filename = filename;

    // Load cube
    m_cube.load(filename);

    // Load energies
    GEnergies energies(filename);

    // Extract number of energy bins
    int num = energies.size();

    // Check if energy binning is consistent with primary image hdu
    if (num != m_cube.nmaps() ) {
        std::string msg = "Number of energies in \"ENERGIES\""
                          " extension ("+gammalib::str(num)+")"
                          " does not match the number of maps ("+
                          gammalib::str(m_cube.nmaps())+" in the"
                          " map cube.\n"
                          "The \"ENERGIES\" extension table shall"
                          " provide one enegy value for each map"
                          " in the cube.";
        throw GException::invalid_value(G_LOAD_CUBE, msg);
    }

    // Set log10(energy) nodes, where energy is in units of MeV
    for (int i = 0; i < num; ++i) {
        m_logE.append(energies[i].log10MeV());
    }

    // Signal that cube has been loaded
    m_loaded = true;

    // Set energy boundaries
    set_energy_boundaries();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy boundaries
 *
 * Computes energy boundaries from the energy nodes. The boundaries will be
 * computed as the centre in log10(energy) between the nodes. If only a
 * single map exists, no boundaries will be computed.
 ***************************************************************************/
void GModelSpatialDiffuseCube::set_energy_boundaries(void)
{
    // Initialise energy boundaries
    m_ebounds.clear();

    // Determine number of energy bins
    int num = m_logE.size();

    // Continue only if there are at least two energy nodes
    if (num > 1) {

        // Loop over all nodes
        for (int i = 0; i < num; ++i) {

            // Calculate minimum energy
            double e_min = (i == 0) ? m_logE[i] - 0.5 * (m_logE[i+1] - m_logE[i])
                                    : m_logE[i] - 0.5 * (m_logE[i] - m_logE[i-1]);

            // Calculate maximum energy
            double e_max = (i < num-1) ? m_logE[i] + 0.5 * (m_logE[i+1] - m_logE[i])
                                       : m_logE[i] + 0.5 * (m_logE[i] - m_logE[i-1]);

            // Set energy boundaries
            GEnergy emin;
            GEnergy emax;
            emin.log10MeV(e_min);
            emax.log10MeV(e_max);

            // Append energy bin to energy boundary arra
            m_ebounds.append(emin, emax);

        } // endfor: looped over energy nodes

    } // endif: there were at least two energy nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update Monte Carlo cache
 *
 * Initialise the cache for Monte Carlo sampling of the map cube. The Monte
 * Carlo cache consists of a linear array that maps a value between 0 and 1
 * into the skymap pixel for all maps in the cube.
 ***************************************************************************/
void GModelSpatialDiffuseCube::update_mc_cache(void)
{
    // Set centre and radius to all sky
    GSkyDir centre;
    double  radius = 360.0;

    // Compute cache
    set_mc_cone(centre, radius);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute cube intensity by log-log interpolation
 *
 * @param[in] photon Incident photon.
 * @return Cube intensity value.
 ***************************************************************************/
double GModelSpatialDiffuseCube::cube_intensity(const GPhoton& photon) const
{
    // Fetch cube
    fetch_cube();

    // Initialise intensity value
    double intensity = 0.0;

    // Continue only if there is energy information for the map cube
    if (m_logE.size() > 0) {

        // Set energy for interpolation in log-energy
        m_logE.set_value(photon.energy().log10MeV());

        // Compute map cube intensity for the left and right map
        double left_intensity  = m_cube(photon.dir(), m_logE.inx_left());
        double right_intensity = m_cube(photon.dir(), m_logE.inx_right());

        // Perform log-log interpolation
        if (left_intensity > 0.0 && right_intensity > 0.0) {
            double log_intensity = m_logE.wgt_left()  * std::log(left_intensity) +
                                   m_logE.wgt_right() * std::log(right_intensity);
            intensity = std::exp(log_intensity);
        }
        else if (left_intensity > 0.0) {
            intensity = left_intensity;
        }
        else if (right_intensity > 0.0) {
            intensity = right_intensity;
        }

    } // endif: energy information was available

    // Return intensity
    return intensity;
}


/***********************************************************************//**
 * @brief Set boundary sky region
 *
 * @todo Implement determination of the cube boundary circle
 ***************************************************************************/
void GModelSpatialDiffuseCube::set_region(void) const
{
    // Set sky region centre to (0,0)
    m_region.centre(0.0, 0.0);

    // Set sky region radius to 180 degrees (all points included)
    m_region.radius(180.0);

    // Return
    return;
}
