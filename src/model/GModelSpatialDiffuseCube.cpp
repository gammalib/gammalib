/***************************************************************************
 *       GModelSpatialDiffuseCube.cpp - Spatial map cube model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDiffuseCube g_spatial_cube_seed;
const GModelSpatialRegistry    g_spatial_cube_registry(&g_spatial_cube_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL                     "GModelSpatialDiffuseCube::eval(GSkyDir&)"
#define G_EVAL_GRADIENTS "GModelSpatialDiffuseCube::eval_gradients(GSkyDir&)"
#define G_MC          "GModelSpatialDiffuseCube::mc(GEnergy&, GTime&, GRan&)"
#define G_READ                 "GModelSpatialDiffuseCube::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialDiffuseCube::write(GXmlElement&)"
#define G_LOAD                 "GModelSpatialDiffuseCube::load(std::string&)"
#define G_ENERGIES           "GModelSpatialDiffuseCube::energies(std::vector"\
                                                                "<GEnergy>&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_CACHE                             //!< Dump MC cache values


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
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
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const std::string& filename,
                                                   const double&      value) :
                          GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set parameter
    m_value.value(value);

    // Load the file
    load(filename);

    // Perform autoscaling of parameter
    autoscale();

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
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const GSkymap&              cube,
                                                   const std::vector<GEnergy>& energies,
                                                   const double&               value) :
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
 * @return Model value.
 *
 * Computes the spatial diffuse model as function of photon parameters.
 ***************************************************************************/
double GModelSpatialDiffuseCube::eval(const GPhoton& photon) const
{
    // Initialise value
    double value = 0.0;

    // Continue only if there is energy information for the map cube
    if (m_logE.size() > 0) {

        // Compute diffuse model value by interpolation in log10(energy)
        m_logE.set_value(photon.energy().log10MeV());
        double intensity = m_logE.wgt_left() *
                           m_cube(photon.dir(), m_logE.inx_left()) +
                           m_logE.wgt_right() *
                           m_cube(photon.dir(), m_logE.inx_right());

        // Set the intensity times the scaling factor as model value
        value = intensity * m_value.value();

    } // endif: energy information was available

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] photon Incident photon.
 * @return Model value.
 *
 * Computes the spatial diffuse model as function of photon parameters and
 * sets the value gradient.
 ***************************************************************************/
double GModelSpatialDiffuseCube::eval_gradients(const GPhoton& photon) const
{
    // Initialise intensity
    double intensity = 0.0;

    // Continue only if there is energy information for the map cube
    if (m_logE.size() > 0) {

        // Compute diffuse model value by interpolation in log10(energy)
        m_logE.set_value(photon.energy().log10MeV());
        intensity = m_logE.wgt_left() *
                    m_cube(photon.dir(), m_logE.inx_left()) +
                    m_logE.wgt_right() *
                    m_cube(photon.dir(), m_logE.inx_right());


    } // endif: energy information was available

    // Compute the model value
    double value = intensity * m_value.value();

	// Compute partial derivatives of the parameter value
	double g_value = (m_value.isfree()) ? intensity * m_value.scale() : 0.0;

	// Set gradient to 0 (circumvent const correctness)
	const_cast<GModelSpatialDiffuseCube*>(this)->m_value.factor_gradient(g_value);

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
 *            No energy boundaries specified, or energy boundaries do not
 *            cover the specified @p energy.
 *
 * Returns a random sky direction according to the intensity distribution of
 * the model sky map and the specified energy. The method makes use of a
 * cache array that contains the normalised cumulative flux values for each
 * of the sky maps in the cube. The specified energy is used to select the
 * appropriate cache array from the cube. Using a uniform random number, the
 * selected cache array is scanned using a bi-section method to determine
 * the skymap pixel for which the position should be returned. To avoid
 * binning problems, the exact position within the pixel is set by a uniform
 * random number generator (neglecting thus pixel distortions). The
 * fractional skymap pixel is then converted into a sky direction.
 ***************************************************************************/
GSkyDir GModelSpatialDiffuseCube::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Allocate sky direction
    GSkyDir dir;

    // Determine number of skymap pixels
    int npix = pixels();

    // Continue only if there are skymap pixels
    if (npix > 0) {

        // If no energy boundaries are defined, throw an exception
        if (m_ebounds.size() < 1) {
            std::string msg = "The energy boundaries of the maps in the cube"
                              " have not been defined. Maybe the map cube file"
                              " is missing the \"ENERGIES\" extension which"
                              " defines the energy of each map in the cube.\n"
                              "Please provide the energy information."; 
            throw GException::invalid_value(G_MC, msg);
        }

        // Determine the map that corresponds best to the specified energy.
        // This is not 100% clean, as ideally some map interpolation should
        // be done to the exact energy specified. However, as long as the map
        // does not change drastically with energy, taking the closest map
        // seems to be fine.
        int i = m_ebounds.index(energy);
        if (i < 0) {
            if (energy <= m_ebounds.emin()) {
                i = 0;
            }
            else if (energy >= m_ebounds.emax()) {
                i = m_ebounds.size()-1;
            }
            else {
                std::string msg = "The specified energy "+energy.print()+" does"
                                  " not fall in any of the energy boundaries of"
                                  " the map cube.\n"
                                  "Please make sure that the map cube energies"
                                  " are properly defined.";
                throw GException::invalid_value(G_MC, msg);
            }
        }
        
        // Get uniform random number
        double u = ran.uniform();

        // Get pixel index according to random number. We use a bi-section
        // method to find the corresponding skymap pixel
        int offset = i * (npix+1);
        int low    = offset;
        int high   = offset + npix;
        while ((high - low) > 1) {
            int mid = (low+high) / 2;
            if (u < m_mc_cache[mid]) {
                high = mid;
            }
            else if (m_mc_cache[mid] <= u) {
                low = mid;
            }
        }

        // Convert 1D pixel index to 2D pixel index
        GSkyPixel pixel = m_cube.pix2xy(low-offset);

        // Randomize pixel
        pixel.x(pixel.x() + ran.uniform() - 0.5);
        pixel.y(pixel.y() + ran.uniform() - 0.5);

        // Get sky direction
        dir = m_cube.xy2dir(pixel);
    
    } // endif: there were pixels in sky map

    // Return sky direction
    return dir;
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
 * Read the map cube information from an XML element. The XML element should
 * have either the format
 *
 *     <spatialModel type="MapCubeFunction" file="test_file.fits">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *     </spatialModel>
 *
 * or alternatively
 *
 *     <spatialModel type="MapCubeFunction" file="test_file.fits">
 *       <parameter name="Value" scale="1" value="1" min="0.1" max="10" free="0"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialDiffuseCube::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Verify that XML element has exactly 1 parameters
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Map cube spatial model requires exactly 1 parameter.");
    }

    // Get pointer on model parameter
    const GXmlElement* par = xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value") {
        m_value.read(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Map cube spatial model requires either \"Value\" or"
              " \"Normalization\" parameter.");
    }

    // Save filename
    m_filename = gammalib::expand_env(xml.attribute("file"));

    // Load the cube
    load(xml.attribute("file"));

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
 *     <spatialModel type="MapCubeFunction" file="test_file.fits">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *     </spatialModel>
 *
 * or alternatively
 *
 *     <spatialModel type="MapCubeFunction" file="test_file.fits">
 *       <parameter name="Value" scale="1" value="1" min="0.1" max="10" free="0"/>
 *     </spatialModel>
 *
 * The latter format is the default for newly written XML elements. 
 ***************************************************************************/
void GModelSpatialDiffuseCube::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "MapCubeFunction");
    }

    // Verify model type
    if (xml.attribute("type") != "MapCubeFunction") {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \"MapCubeFunction\".");
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
    xml.attribute("file", m_filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load cube into the model class
 *
 * @param[in] filename cube file.
 *
 * @exception GException::invalid_value
 *            Number of maps in cube mismatches number of energy bins.
 *
 * Loads cube into the model class.
 ***************************************************************************/
void GModelSpatialDiffuseCube::load(const std::string& filename)
{
    // Initialise skymap
    m_cube.clear();
    m_logE.clear();

    // Store filename of cube (for XML writing). Note that we do not
    // expand any environment variable at this level, so that if we write
    // back the XML element we write the filepath with the environment
    // variables
    m_filename = filename;

    // Load cube
    m_cube.load(gammalib::expand_env(m_filename));

    // Load the energy binning definition
    GFits file(gammalib::expand_env(m_filename));

    // Get the table extension "ENERGIES"
    GFitsTable* energies = file.table("ENERGIES");

    // Read only if extension "ENERGIES" exists
    if (energies != NULL) {

        // Extract number of energy bins in FITS file
        int num = energies->integer("NAXIS2");

        // Check if energy binning is consistent with primary image hdu
        if (num != m_cube.nmaps() ) {
            std::string msg = "Number of energies in \"ENERGIES\" extension"
                              " ("+gammalib::str(num)+") does not match the"
                              " number of maps ("+gammalib::str(m_cube.nmaps())+""
                              " in the map cube.\n"
                              "The \"ENERGIES\" extension table shall provide"
                              " one enegy value for each map in the cube.";
            throw GException::invalid_value(G_LOAD, msg);
        }

        // Get the unit of the energies. If no TUNIT1 header keyword is found
        // then use MeV
        std::string unit = "MeV";
        if (energies->header()->hascard("TUNIT1")) {
            unit = energies->string("TUNIT1");
        }

        // Get the column with the name "Energy"
        const GFitsTableCol& col_energy = (*energies)["Energy"];

        // Set log10(energy) nodes, where energy is in units of MeV
        for (int i = 0; i < num; ++i) {
            GEnergy energy(col_energy.real(i), unit);
            m_logE.append(energy.log10MeV());
        }

    } // endif: ENERGIES extension existed

    // Set energy boundaries
    set_energy_boundaries();

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
void GModelSpatialDiffuseCube::energies(const std::vector<GEnergy>& energies)
{
    // Initialise energies
    m_logE.clear();

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
std::vector<GEnergy> GModelSpatialDiffuseCube::energies(void)
{
    // Initialise vector of energies
    std::vector<GEnergy> energies;

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
            energies.push_back(energy);
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
 * that will be simulated using the mc() method.
 ***************************************************************************/
void GModelSpatialDiffuseCube::set_mc_cone(const GSkyDir& centre,
                                           const double&  radius)
{
    // Initialise cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();

    // Determine number of cube pixels and maps
    int npix  = pixels();
    int nmaps = maps();

    // Continue only if there are pixels and maps
    if (npix > 0 && nmaps > 0) {

        // Reserve space for all pixels in cache
        m_mc_cache.reserve((npix+1)*nmaps);

        // Loop over all maps
        for (int i = 0; i < nmaps; ++i) {

            // Compute pixel offset
            int offset = i * (npix+1);

            // Set first cache value to 0
            m_mc_cache.push_back(0.0);

            // Initialise cache with cumulative pixel fluxes and compute
            // total flux in skymap for normalization. Negative pixels are
            // excluded from the cumulative map.
            double flux_centre = 0.0;
            double flux_circle = 0.0;
        	for (int k = 0; k < npix; ++k) {

                // Derive effective pixel radius from half opening angle
                // that corresponds to the pixel's solid angle
                double pixel_radius =
                       std::acos(1.0 - m_cube.omega(k)/gammalib::twopi) *
                       gammalib::rad2deg * 1.5;

                // Add up flux. First we check if the pixel is within the
                // safety radius (simulation cone + effective pixel radius),
                // then we check if the pixels is within the simulation cone
                double distance = centre.dist_deg(m_cube.pix2dir(k));
                if (distance <= radius+pixel_radius) {
                    double flux = m_cube(k,i) * m_cube.omega(k);
                    if (flux > 0.0) {
                        flux_circle += flux;
                        if (distance <= radius) {
                            flux_centre += flux;
                        }
                    }
                }

                // Push back flux derived using safety radius
        		m_mc_cache.push_back(flux_circle); // units: ph/cm2/s/MeV
        	}

            // Normalize cumulative pixel fluxes so that the values in the
            // cache run from 0 to 1
            if (flux_circle > 0.0) {
        		for (int k = 0; k < npix; ++k) {
        			m_mc_cache[k+offset] /= flux_circle;
        		}
        	}

            // Make sure that last pixel in the cache is >1
            m_mc_cache[npix+offset] = 1.0001;

            // Store centre flux in node array
            if (m_logE.size() == nmaps) {
                GEnergy energy;
                energy.log10MeV(m_logE[i]);

                // Only append new node if flux > 0
                if (flux_centre>0.0) {
                	m_mc_spectrum.append(energy, flux_centre);
                }
            }

        } // endfor: looped over all maps

        // Dump cache values for debugging
        #if defined(G_DEBUG_CACHE)
        for (int i = 0; i < m_mc_cache.size(); ++i) {
            std::cout << "i=" << i;
            std::cout << " c=" << m_mc_cache[i] << std::endl;
        }
        #endif

    } // endif: there were cube pixels and maps

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
        result.append("\n"+gammalib::parformat("Map cube file")+m_filename);
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

        // NORMAL: Append sky map
        if (chatter >= NORMAL) {
            result.append("\n"+m_cube.print(chatter));
        }

        // EXPLICIT: Append energy nodes
        if (chatter >= EXPLICIT) {
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
    // Initialise Value
    m_value.clear();
    m_value.name("Normalization");
    m_value.value(1.0);
    m_value.scale(1.0);
    m_value.range(0.001, 1000.0);
    m_value.gradient(0.0);
    m_value.fix();
    m_value.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Initialise other members
    m_filename.clear();
    m_cube.clear();
    m_logE.clear();
    m_ebounds.clear();

    // Initialise MC cache
    m_mc_cache.clear();
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
    m_value    = model.m_value;
    m_filename = model.m_filename;
    m_cube     = model.m_cube;
    m_logE     = model.m_logE;
    m_ebounds  = model.m_ebounds;

    // Copy MC cache
    m_mc_cache    = model.m_mc_cache;
    m_mc_spectrum = model.m_mc_spectrum;

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
