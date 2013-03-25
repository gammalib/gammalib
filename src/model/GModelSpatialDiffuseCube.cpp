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

    // Set filename
    m_filename = filename;

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
    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_EVAL);

    // Return value
    return 1.0;
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
    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_EVAL_GRADIENTS);

    // Set gradient to 0 (circumvent const correctness)
    const_cast<GModelSpatialDiffuseCube*>(this)->m_value.factor_gradient(0.0);

    // Return value
    return 1.0;
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
    m_filename = expand_env(xml.attribute("file"));

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
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialDiffuseCube::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialDiffuseCube ===");

    // Append parameters
    result.append("\n"+parformat("Map cube file")+m_filename);
    result.append("\n"+parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+m_pars[i]->print());
    }

    // Append sky map
    if (m_loaded) {
        result.append("\n"+m_cube.print());
    }

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
