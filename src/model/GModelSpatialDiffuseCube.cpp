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
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GFITSTableDoubleCol.hpp"


/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDiffuseCube g_spatial_cube_seed;
const GModelSpatialRegistry    g_spatial_cube_registry(&g_spatial_cube_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL                     "GModelSpatialDiffuseCube::eval(GSkyDir&)"
#define G_EVAL_GRADIENTS "GModelSpatialDiffuseCube::eval_gradients(GSkyDir&)"
#define G_MC                            "GModelSpatialDiffuseCube::mc(GRan&)"
#define G_READ                 "GModelSpatialDiffuseCube::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialDiffuseCube::write(GXmlElement&)"
#define G_LOAD               "GModelSpatialDiffuseCube::load(std::string& filename)"

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
 * Constructs map cube model by assigning the normalization @p value and the
 * @p filename of the map cube. Note that the map cube file is not opened,
 * only the filename is stored for use when required.
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
 * @param[in] map Sky map.
 * @param[in] value Normalization factor (defaults to 1).
 *
 * Constructs map cube model by assigning the normalisation @p value and by
 * loading a @p map cube from a sky map. The filename will remain blank.
 ***************************************************************************/
GModelSpatialDiffuseCube::GModelSpatialDiffuseCube(const GSkymap& map,
                                                   const double&  value) :
                          GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set parameter
    m_value.value(value);

    // Perform autoscaling of parameter
    autoscale();

    // Set map cube
    cube(map);

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
GModelSpatialDiffuseCube& GModelSpatialDiffuseCube::operator= (const GModelSpatialDiffuseCube& model)
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
 * @exception GException::feature_not_implemented
 *            Method not yet implemented
 *
 * @todo Implement method.
 ***************************************************************************/
double GModelSpatialDiffuseCube::eval(const GPhoton& photon) const
{
	// Get skymap bin
	int imap = m_ebounds.index(photon.energy());

	// Get skymap intensity
	double intensity = m_cube(photon.dir(), imap);

	// Return intensity times normalization factor
	return (intensity * m_value.value());

}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] photon Incident photon.
 * @return Model value.
 *
 * @todo Implement method.
 ***************************************************************************/
double GModelSpatialDiffuseCube::eval_gradients(const GPhoton& photon) const
{
	// Get skymap bin
	int imap = m_ebounds.index(photon.energy());

	// Get skymap intensity
	double intensity = m_cube(photon.dir(), imap);

	// Compute partial derivatives of the parameter value
	double g_value = (m_value.isfree()) ? intensity * m_value.scale() : 0.0;

	// Set gradient to 0 (circumvent const correctness)
	const_cast<GModelSpatialDiffuseCube*>(this)->m_value.factor_gradient(g_value);

	// Return intensity times normalization factor
	return (intensity * m_value.value());

}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented
 *
 * @todo Implement method.
 ***************************************************************************/
GSkyDir GModelSpatialDiffuseCube::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Allocate sky direction
    GSkyDir dir;

    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_MC);
    
    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief _ model from XML element
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

        // Append sky map
        if (m_loaded) {
            result.append("\n"+m_cube.print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}

/***********************************************************************//**
 * @brief Load cube into the model class
 *
 * @param[in] filename cube file.
 *
 * Loads cube into the model class. The method calls the protected method
 * prepare_cube() that prepares the cube for usage by the class.
 ***************************************************************************/
void GModelSpatialDiffuseCube::load(const std::string& filename)
{
    // Initialise skymap
    m_cube.clear();
    m_ebounds.clear();
    m_loaded = false;

    // Store filename of cube (for XML writing). Note that we do not
    // expand any environment variable at this level, so that if we write
    // back the XML element we write the filepath with the environment
    // variables
    m_filename = filename;

    // Load cube
    m_cube.load(gammalib::expand_env(m_filename));

    // Loading the energy binning definition
    GFits file(filename);

    // Get the table extension "ENERGIES"
    GFitsTable* energies = file.table("ENERGIES");

    // Read only if extension "ENERGIES" exists
    if (energies != NULL) {

        // Extract number of energy bins in FITS file
        int num = energies->integer("NAXIS2");

        // Check if energy binning is consistent with primary image hdu
        if (num != m_cube.nmaps() ) {

            // if inconsistent throw an exception
            // ToDo: Excepton has to be adjusted
            throw GException::skymap_bad_size(G_LOAD,num,m_cube.nmaps());
        }

        // Check if there are enough energy nodes to construct binning
        if (num < 2) {
            // if inconsistent throw an exception
            // ToDo: Excepton has to be adjusted
            throw GException::skymap_bad_size(G_LOAD,num,2);
        }

        // Get the column with the name "Energy"
        GFitsTableDoubleCol* ptr_energy = (GFitsTableDoubleCol*)&(*energies)["Energy"];
        ptr_energy->print();
        
        // Get the unit of the energies
        // Default for Fermi Galactic diffuse model
        std::string unit = "MeV";
        if (energies->header()->hascard("TUNIT1")) {
            unit = energies->string("TUNIT1");
        }

        // Read the energy binning
        // We construct energy bins from the given bin centers
        // bin edges will be computed in the logarithimic center
        // between the nodes
        for (int i = 0; i < num; ++i) {

            double emin = 0.0;
            double emax = 0.0;
            double step = 0.0;
            // calculate bin width
            // if first bin, use second bin to calculate stepsize
            if (i == 0) {
                step = std::log10((*ptr_energy)(i + 1)) - std::log10((*ptr_energy)(i));
            }

            // else use previous bin to get stepsize
            else {
                step = std::log10((*ptr_energy)(i)) - std::log10((*ptr_energy)(i - 1));
            }

            // calculate bin minimum and bin maximum
            emin = std::log10((*ptr_energy)(i)) - step / 2.0;
            emax = std::log10((*ptr_energy)(i)) + step / 2.0;

            // Convert it to a GEnergy without logscale
            GEnergy bin_min = GEnergy(std::pow(10.0,emin),unit);
            GEnergy bin_max = GEnergy(std::pow(10.0,emax),unit);

            // Append energy bin to Ebounds
            m_ebounds.append(bin_min,bin_max);

        } // endfor: looped over energy bins

    } // endif: ENERGIES extension existed

    // Prepare cube
    prepare_cube();

    // Cube is loaded
    m_loaded = true;

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
    m_ebounds.clear();
    m_loaded = false;

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
    m_loaded   = model.m_loaded;
    m_ebounds = model.m_ebounds;

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
 * @brief Prepare cube after loading
 *
 * Prepares a cube after loading. The cube is normalised so that the total
 * flux of the cube amounts to 1 ph/cm2/s. Negative skymap pixels are set to
 * zero intensity.
 *
 * The method also initialises a cache for Monte Carlo sampling of the
 * cube. This Monte Carlo cache consists of a linear array that maps a
 * value between 0 and 1 into the skymap pixel.
 *
 * Note that if the GSkymap object contains multiple maps, only the first
 * map is used.
 ***************************************************************************/
void GModelSpatialDiffuseCube::prepare_cube(void)
{
    // Initialise cache
    m_mc_cache.clear();

    // Determine number of cube pixels
    int npix = m_cube.npix();
    int nmaps =  m_cube.nmaps();

    // Continue only if there are pixels
    if (npix > 0) {

        // Reserve space for all pixels in cache
        m_mc_cache.reserve(npix+1);

        // Set first cache value to 0
        m_mc_cache.push_back(0.0);

        // Initialise cache with cumulative pixel fluxes and compute total
        // flux in skymap for normalization. Negative pixels are set to
        // zero intensity in the skymap.
        double sum = 0.0;
        for (int i=0;i< nmaps; ++i) {
        	for (int j = 0; j < npix; ++j) {
        		double flux = m_cube(j,i) * m_cube.omega(i);
        		if (flux < 0.0) {
        		    m_cube(j,i) = 0.0;
        		    flux     = 0.0;
        		}
        		sum += flux;
        		m_mc_cache.push_back(sum);
        	}
        }

        // Normalize skymap and pixel fluxes in the cache so that the values
        // in the cache run from 0 to 1
        if (sum > 0.0) {
        	for (int i = 0; i < nmaps; ++i) {
        		for (int j = 0; j < npix; ++j) {
        			m_cube(j, i)      /= sum;
        			m_mc_cache[i] /= sum;
        		}
        	}
        }

        // Make sure that last pixel in the cache is >1
        m_mc_cache[npix] = 1.0001;

        // Dump preparation results
        #if defined(G_DEBUG_PREPARE)
        double sum_control = 0.0;
        for (int i = 0; i < nmaps; ++i) {
        	for (int j = 0; j < npix; ++j) {
        		double flux = m_cube(j,i) * m_cube.omega(i);
        		if (flux >= 0.0) {
        			sum_control += flux;
        		}
        	}
        }
        std::cout << "Total flux before normalization: " << sum << std::endl;
        std::cout << "Total flux after normalization : " << sum_control << std::endl;
        #endif

        // Dump cache values for debugging
        #if defined(G_DEBUG_CACHE)
        for (int i = 0; i < m_mc_cache.size(); ++i) {
            std::cout << "i=" << i;
            std::cout << " c=" << m_mc_cache[i] << std::endl;
        }
        #endif

    } // endif: there were cube pixels

    // Return
    return;
}
