/***************************************************************************
 *        GModelSpectralFunc.cpp  -  Spectral function model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GModelSpectralFunc.cpp
 * @brief Spectral function model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GCsv.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralFunc     g_spectral_func_seed;
const GModelSpectralRegistry g_spectral_func_registry(&g_spectral_func_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                 "GModelSpectralFunc::flux(GEnergy&, GEnergy&)"
#define G_EFLUX               "GModelSpectralFunc::eflux(GEnergy&, GEnergy&)"
#define G_MC              "GModelSpectralFunc::mc(GEnergy&, GEnergy&, GRan&)"
#define G_READ                       "GModelSpectralFunc::read(GXmlElement&)"
#define G_WRITE                     "GModelSpectralFunc::write(GXmlElement&)"
#define G_LOAD_NODES           "GModelSpectralFunc::load_nodes(std::string&)"

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
GModelSpectralFunc::GModelSpectralFunc(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name of nodes.
 *
 * Creates instance of a spectral function model from a list of nodes that is
 * found in the specified file. See GModelSpectralFunc::load_nodes() for more
 * information about the expected structure of the file.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const std::string& filename) :
                    GModelSpectral()
{
    // Initialise members
    init_members();

    // Load nodes
    load_nodes(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a spectral function model by extracting information
 * from an XML element. See GModelSpectralFunc::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const GXmlElement& xml) :
                    GModelSpectral()
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
 * @param[in] model Spectral function model.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const GModelSpectralFunc& model) :
                    GModelSpectral(model)
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
GModelSpectralFunc::~GModelSpectralFunc(void)
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
 * @param[in] model Spectral function model.
 ***************************************************************************/
GModelSpectralFunc& GModelSpectralFunc::operator= (const GModelSpectralFunc& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GModelSpectralFunc::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpectralFunc* GModelSpectralFunc::clone(void) const
{
    return new GModelSpectralFunc(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 * Note that the node energies are stored as log10 of energy in units of
 * MeV.
 ***************************************************************************/
double GModelSpectralFunc::eval(const GEnergy& srcEng) const
{
    // Interpolate function
    double func = m_nodes.interpolate(srcEng.log10MeV(), m_values);

    // Compute function value
    double value  = norm() * func;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralFunc::eval";
        std::cout << "(srcEng=" << srcEng << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", norm=" << norm();
        std::cout << ", func=" << func;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization of the function.
 * Note that the normalization is factorised into a scaling factor and a
 * value and that the method is expected to return the gradient with respect
 * to the parameter value \f$n_v\f$.
 *
 * The partial derivative of the normalization value is given by
 * \f[dI/dn_v=n_s f(E)\f]
 ***************************************************************************/
double GModelSpectralFunc::eval_gradients(const GEnergy& srcEng) const
{
    // Interpolate function
    double func = m_nodes.interpolate(srcEng.log10MeV(), m_values);

    // Compute function value
    double value  = norm() * func;

    // Compute partial derivatives of the parameter values
    double g_norm  = (m_norm.isfree())  ? m_norm.scale() * func : 0.0;

    // Set gradients (circumvent const correctness)
    ((GModelSpectralFunc*)this)->m_norm.gradient(g_norm);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralFunc::eval_gradients";
        std::cout << "(srcEng=" << srcEng << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", norm=" << norm();
        std::cout << ", func=" << func;
        std::cout << ", g_norm=" << g_norm;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 *
 * @todo Implement method
 ***************************************************************************/
double GModelSpectralFunc::flux(const GEnergy& emin, const GEnergy& emax) const
{
    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_FLUX);

    // Return
    return 0.0;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) E dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 *
 * @todo Implement method
 ***************************************************************************/
double GModelSpectralFunc::eflux(const GEnergy& emin, const GEnergy& emax) const
{
    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_EFLUX);

    // Return
    return 0.0;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] ran Random number generator.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented
 *
 * @todo To be implemented
 ***************************************************************************/
GEnergy GModelSpectralFunc::mc(const GEnergy& emin, const GEnergy& emax,
                               GRan& ran) const
{
    // Allocate energy
    GEnergy energy;

    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_MC);

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Read the function information from an XML element and load the nodes
 * from the associated file. The XML element is required to have an
 * attribute "file" that specifies the function nodes and a parameter
 * named "Normalization".
 ***************************************************************************/
void GModelSpectralFunc::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Spectral function requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization") {
        m_norm.read(*par);
    }
    else
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require \"Normalization\" parameter.");

    // Load nodes from file
    load_nodes(xml.attribute("file"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "FileFunction"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the spectral function information into an XML element. The XML
 * element has to be of type "FileFunction" and will have 1 parameter leaf
 * named "Normalization". Note that the function nodes will not be written
 * since they will not be altered by any method.
 ***************************************************************************/
void GModelSpectralFunc::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "FileFunction");

    // Verify model type
    if (xml.attribute("type") != "FileFunction")
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"FileFunction\".");

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"Normalization\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral function requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Set parameyter
    if (par->attribute("name") == "Normalization") {
        m_norm.write(*par);
    }
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Require \"Normalization\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelSpectralFunc::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralFunc ===\n");
    result.append(parformat("Function file")+m_filename);
    result.append(parformat("Number of nodes")+str(m_nodes.size()));
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralFunc::init_members(void)
{
    // Initialise powerlaw normalisation
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.scale(1.0);
    m_norm.value(1.0);
    m_norm.range(0.0,1000.0);
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Initialise other members
    m_nodes.clear();
    m_values.clear();
    m_filename.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral function model.
 ***************************************************************************/
void GModelSpectralFunc::copy_members(const GModelSpectralFunc& model)
{
    // Copy members
    m_norm     = model.m_norm;
    m_nodes    = model.m_nodes;
    m_values   = model.m_values;
    m_filename = model.m_filename;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralFunc::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load nodes from file
 *
 * @param[in] filename File name.
 *
 * @exception GException::not_enough_data
 *            File contains less than 2 nodes
 * @exception GException::not_enough_columns
 *            File contains less than 2 columns
 *
 * The file function is stored as a column separated value table (CSV) in an
 * ASCII file with (at least) 2 columns. The first column specifies the
 * energy in MeV while the second column specifies the intensity at this
 * energy in units of ph/cm2/s/MeV.
 * The node energies will be stored as log10 of the energy given in units
 * of MeV. At least 2 nodes and 2 columns are required.
 ***************************************************************************/
void GModelSpectralFunc::load_nodes(const std::string& filename)
{
    // Clear nodes and values
    m_nodes.clear();
    m_values.clear();

    // Set filename
    m_filename = filename;

    // Load file
    GCsv csv = GCsv(filename);

    // Check if there are at least 2 nodes
    if (csv.nrows() < 2)
        throw GException::not_enough_data(G_LOAD_NODES, filename,
                                          csv.nrows());

    // Check if there are at least 2 columns
    if (csv.ncols() < 2)
        throw GException::not_enough_columns(G_LOAD_NODES, filename,
                                             csv.ncols());

    // Setup nodes
    for (int i = 0; i < csv.nrows(); ++i) {
        m_nodes.append(std::log10(csv.real(i,0)));
        m_values.push_back(csv.real(i,1));
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
