/***************************************************************************
 *         GModelSpectralFunc.cpp - Spectral function model class          *
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
 * @file GModelSpectralFunc.cpp
 * @brief Spectral function model class implementation
 * @author Juergen Knoedlseder
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
#define G_MC      "GModelSpectralFunc::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
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
 * @param[in] norm Normalization factor.
 *
 * Constructs spectral file function model from a list of nodes that is found
 * in the specified ASCII file. See the load_nodes() method for more
 * information about the expected structure of the file.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const std::string& filename,
                                       const double&      norm) :
                    GModelSpectral()
{
    // Initialise members
    init_members();

    // Load nodes
    load_nodes(filename);

    // Set normalization
    m_norm.value(norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs spectral file function model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
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
 * @param[in] model File function model.
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
 * @param[in] model File function model.
 * @return File function model.
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
 * @brief Clear file function
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
 * @brief Clone file function
***************************************************************************/
GModelSpectralFunc* GModelSpectralFunc::clone(void) const
{
    // Clone file function
    return new GModelSpectralFunc(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization factor.
 ***************************************************************************/
double GModelSpectralFunc::eval(const GEnergy& srcEng,
                                const GTime&   srcTime) const
{
    // Interpolate function. This is done in log10-log10 space, but the
    // linear value is returned.
    double arg  = m_log_nodes.interpolate(srcEng.log10MeV(), m_log_values);
    double func = std::pow(10.0, arg);

    // Compute function value
    double value  = m_norm.value() * func;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralFunc::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
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
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization factor.
 *
 * The method also evaluates the partial derivatives of the model with
 * respect to the normalization parameter using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_norm}} =
 *      \frac{S_{\rm E}(E | t)}{{\tt m\_norm}}
 * \f]
 ***************************************************************************/
double GModelSpectralFunc::eval_gradients(const GEnergy& srcEng,
                                          const GTime&   srcTime)
{
    // Interpolate function. This is done in log10-log10 space, but the
    // linear value is returned.
    double arg  = m_log_nodes.interpolate(srcEng.log10MeV(), m_log_values);
    double func = std::pow(10.0, arg);

    // Compute function value
    double value  = m_norm.value() * func;

    // Compute partial derivatives of the parameter values
    double g_norm  = (m_norm.isfree())  ? m_norm.scale() * func : 0.0;

    // Set gradients
    m_norm.factor_gradient(g_norm);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralFunc::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
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
 * @return Photon flux (ph/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralFunc::flux(const GEnergy& emin, const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();
    
        // Determine left node index for minimum energy
        m_lin_nodes.set_value(e_min);
        int inx_emin = m_lin_nodes.inx_left();

        // Determine left node index for maximum energy
        m_lin_nodes.set_value(e_max);
        int inx_emax = m_lin_nodes.inx_left();
    
        // If both energies are within the same nodes then simply
        // integrate over the energy interval using the appropriate power
        // law parameters
        if (inx_emin == inx_emax) {
            flux = m_prefactor[inx_emin] * 
                   plaw_photon_flux(e_min,
                                    e_max, 
                                    m_epivot[inx_emin],
                                    m_gamma[inx_emin]);
        }
    
        // ... otherwise integrate over the nodes where emin and emax
        // resides and all the remaining nodes
        else {

            // If we are requested to extrapolate beyond first node,
            // the use the first nodes lower energy and upper energy
            // boundary for the initial integration.
            int i_start = (e_min < m_lin_nodes[0]) ? inx_emin : inx_emin+1;

            // Integrate from emin to the node boundary
            flux = m_prefactor[inx_emin] *
                   plaw_photon_flux(e_min,
                                    m_lin_nodes[i_start],
                                    m_epivot[inx_emin],
                                    m_gamma[inx_emin]);

            // Integrate over all nodes between
            for (int i = i_start; i < inx_emax; ++i) {
                flux += m_flux[i];
            }

            // Integrate from node boundary to emax
            flux += m_prefactor[inx_emax] *
                    plaw_photon_flux(m_lin_nodes[inx_emax],
                                     e_max,
                                     m_epivot[inx_emax],
                                     m_gamma[inx_emax]);
        
        } // endelse: emin and emax not between same nodes
    
        // Multiply flux by normalisation factor
        flux *= m_norm.value();
    
    } // endif: integration range was valid

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Energy flux (erg/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) E \, dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralFunc::eflux(const GEnergy& emin, const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {
    
        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();
    
        // Determine left node index for minimum energy
        m_lin_nodes.set_value(e_min);
        int inx_emin = m_lin_nodes.inx_left();

        // Determine left node index for maximum energy
        m_lin_nodes.set_value(e_max);
        int inx_emax = m_lin_nodes.inx_left();
    
        // If both energies are within the same nodes then simply
        // integrate over the energy interval using the appropriate power
        // law parameters
        if (inx_emin == inx_emax) {
            eflux = m_prefactor[inx_emin] * 
                    plaw_energy_flux(e_min,
                                     e_max, 
                                     m_epivot[inx_emin],
                                     m_gamma[inx_emin]) * gammalib::MeV2erg;
        }

        // ... otherwise integrate over the nodes where emin and emax
        // resides and all the remaining nodes
        else {

            // If we are requested to extrapolate beyond first node,
            // the use the first nodes lower energy and upper energy
            // boundary for the initial integration.
            int i_start = (e_min < m_lin_nodes[0]) ? inx_emin : inx_emin+1;

            // Integrate from emin to the node boundary
            eflux = m_prefactor[inx_emin] *
                    plaw_energy_flux(e_min,
                                     m_lin_nodes[i_start],
                                     m_epivot[inx_emin],
                                     m_gamma[inx_emin]) * gammalib::MeV2erg;

            // Integrate over all nodes between
            for (int i = i_start; i < inx_emax; ++i) {
                eflux += m_eflux[i];
            }

            // Integrate from node boundary to emax
            eflux += m_prefactor[inx_emax] *
                     plaw_energy_flux(m_lin_nodes[inx_emax],
                                      e_max,
                                      m_epivot[inx_emax],
                                      m_gamma[inx_emax]) * gammalib::MeV2erg;
        
        } // endelse: emin and emax not between same nodes
    
        // Multiply flux by normalisation factor
        eflux *= norm();
    
    } // endif: integration range was valid

    // Return flux
    return eflux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] time True photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Energy.
 *
 * @exception GException::erange_invalid
 *            Energy range is invalid (emin < emax required).
 *
 * Returns Monte Carlo energy by randomly drawing from a broken power law
 * defined by the file function.
 ***************************************************************************/
GEnergy GModelSpectralFunc::mc(const GEnergy& emin,
                               const GEnergy& emax,
                               const GTime&   time,
                               GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
              "Minimum energy < maximum energy required.");
    }

    // Allocate energy
    GEnergy energy;
    
    // Update cache
    mc_update(emin, emax);

    // Determine in which bin we reside
    int inx = 0;
    if (m_mc_cum.size() > 1) {
        double u = ran.uniform();
        for (inx = m_mc_cum.size()-1; inx > 0; --inx) {
            if (m_mc_cum[inx-1] <= u) {
                break;
            }
        }
    }

    // Get random energy for specific bin
    if (m_mc_exp[inx] != 0.0) {
        double e_min = m_mc_min[inx];
        double e_max = m_mc_max[inx];
        double u     = ran.uniform();
        double eng   = (u > 0.0) 
                        ? std::exp(std::log(u * (e_max - e_min) + e_min) / m_mc_exp[inx])
                        : 0.0;
        energy.MeV(eng);
    }
    else {
        double e_min = m_mc_min[inx];
        double e_max = m_mc_max[inx];
        double u     = ran.uniform();
        double eng   = std::exp(u * (e_max - e_min) + e_min);
        energy.MeV(eng);
    }
    
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
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *     <spectrum type="FileFunction" file="..">
 *       <parameter name="Normalization" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralFunc::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Spectral function requires exactly 1 parameter.");
    }

    // Get parameter element
    const GXmlElement* par = xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization") {
        m_norm.read(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require \"Normalization\" parameter.");
    }

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
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="FileFunction" file="..">
 *       <parameter name="Normalization" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 * 
 * Note that the function nodes will not be written since they will not be
 * altered by any method.
 ***************************************************************************/
void GModelSpectralFunc::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "FileFunction");
    }

    // Verify model type
    if (xml.attribute("type") != "FileFunction") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"FileFunction\".");
    }

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Normalization\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral function requires exactly 1 parameter.");
    }

    // Get parameter element
    GXmlElement* par = xml.element("parameter", 0);

    // Set parameter
    if (par->attribute("name") == "Normalization") {
        m_norm.write(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Require \"Normalization\" parameter.");
    }

    // Set file attribute
    xml.attribute("file", m_filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print file function information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing file function information.
 ***************************************************************************/
std::string GModelSpectralFunc::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralFunc ===");

        // Append information
        result.append("\n"+parformat("Function file")+m_filename);
        result.append("\n"+parformat("Number of nodes")+str(m_lin_nodes.size()));
        result.append("\n"+parformat("Number of parameters")+str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

        // Append node information
        for (int i = 0; i < m_prefactor.size(); ++i) {
            result.append("\n"+parformat("Node "+str(i+1)));
            result.append("Epivot="+str(m_epivot[i]));
            result.append(" Prefactor="+str(m_prefactor[i]));
            result.append(" Gamma="+str(m_gamma[i]));
            result.append(" Flux="+str(m_flux[i]));
            result.append(" EFlux="+str(m_eflux[i]));
        }

    } // endif: chatter was not silent

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

    // Initialise members
    m_lin_nodes.clear();
    m_log_nodes.clear();
    m_lin_values.clear();
    m_log_values.clear();
    m_filename.clear();
    
    // Initialise cache
    m_mc_emin.clear();
    m_mc_emax.clear();
    m_mc_cum.clear();
    m_mc_min.clear();
    m_mc_max.clear();
    m_mc_exp.clear();

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
    m_norm       = model.m_norm;
    m_lin_nodes  = model.m_lin_nodes;
    m_log_nodes  = model.m_log_nodes;
    m_lin_values = model.m_lin_values;
    m_log_values = model.m_log_values;
    m_filename   = model.m_filename;
    
    // Copy pre-computation cache
    m_prefactor  = model.m_prefactor;
    m_gamma      = model.m_gamma;
    m_epivot     = model.m_epivot;
    m_flux       = model.m_flux;
    m_eflux      = model.m_eflux;

    // Copy MC cache
    m_mc_emin    = model.m_mc_emin;
    m_mc_emax    = model.m_mc_emax;
    m_mc_cum     = model.m_mc_cum;
    m_mc_min     = model.m_mc_min;
    m_mc_max     = model.m_mc_max;
    m_mc_exp     = model.m_mc_exp;

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
 * @exception GException::file_function_data
 *            File contains less than 2 nodes
 * @exception GException::file_function_columns
 *            File contains less than 2 columns
 * @exception GException::file_function_value
 *            File contains invalid value
 *
 * The file function is stored as a column separated value table (CSV) in an
 * ASCII file with (at least) 2 columns. The first column specifies the
 * energy in MeV while the second column specifies the intensity at this
 * energy in units of ph/cm2/s/MeV.
 * The node energies and values will be stored both linearly and as log10.
 * The log10 storing requires that node energies and node values are
 * positive. Also, at least 2 nodes and 2 columns are required in the file
 * function.
 ***************************************************************************/
void GModelSpectralFunc::load_nodes(const std::string& filename)
{
    // Clear nodes and values
    m_lin_nodes.clear();
    m_log_nodes.clear();
    m_lin_values.clear();
    m_log_values.clear();

    // Set filename
    m_filename = filename;

    // Load file
    GCsv csv = GCsv(filename);

    // Check if there are at least 2 nodes
    if (csv.nrows() < 2) {
        throw GException::file_function_data(G_LOAD_NODES, filename,
                                             csv.nrows());
    }

    // Check if there are at least 2 columns
    if (csv.ncols() < 2) {
        throw GException::file_function_columns(G_LOAD_NODES, filename,
                                                csv.ncols());
    }

    // Setup nodes
    double last_energy = 0.0;
    for (int i = 0; i < csv.nrows(); ++i) {
    
        // Get log10 of node energy and value. Make sure they are valid.
        double log10energy;
        double log10value;
        if (csv.real(i,0) > 0) {
            log10energy = std::log10(csv.real(i,0));
        }
        else {
            throw GException::file_function_value(G_LOAD_NODES, filename,
                  csv.real(i,0), "Energy value must be positive.");
        }
        if (csv.real(i,1) > 0) {
            log10value = std::log10(csv.real(i,1));
        }
        else {
            throw GException::file_function_value(G_LOAD_NODES, filename,
                  csv.real(i,1), "Intensity value must be positive.");
        }
        
        // Make sure that energies are increasing
        if (csv.real(i,0) <= last_energy) {
            throw GException::file_function_value(G_LOAD_NODES, filename,
                  csv.real(i,0), "Energy values must be monotonically increasing.");
        }
        
        // Append log10 of node energy and value    
        m_lin_nodes.append(csv.real(i,0));
        m_log_nodes.append(log10energy);
        m_lin_values.push_back(csv.real(i,1));
        m_log_values.push_back(log10value);
        
        // Store last energy for monotonically increasing check
        last_energy = csv.real(i,0);
        
    } // endfor: looped over nodes
    
    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pre-computation cache
 ***************************************************************************/
void GModelSpectralFunc::set_cache(void) const
{
    // Clear any existing values
    m_prefactor.clear();
    m_gamma.clear();
    m_epivot.clear();
    m_flux.clear();
    m_eflux.clear();
    
    // Loop over all nodes-1
    for (int i = 0; i < m_lin_nodes.size()-1; ++i) {
    
        // Get energies and function values
        double emin = m_lin_nodes[i];
        double emax = m_lin_nodes[i+1];
        double fmin = m_lin_values[i];
        double fmax = m_lin_values[i+1];
    
        // Compute pivot energy (MeV). We use here the geometric mean of
        // the node boundaries.
        double epivot = std::sqrt(emin*emax);
        
        // Compute spectral index
        double gamma = std::log(fmin/fmax) / std::log(emin/emax);
        
        // Compute power law normalisation
        double prefactor = fmin / std::pow(emin/epivot, gamma);
        
        // Compute photon flux between nodes
        double flux = prefactor*plaw_photon_flux(emin, emax, epivot, gamma);

        // Compute energy flux between nodes
        double eflux = prefactor*plaw_energy_flux(emin, emax, epivot, gamma);

        // Convert energy flux from MeV/cm2/s to erg/cm2/s
        eflux *= gammalib::MeV2erg;
        
        // Push back values on pre-computation cache
        m_prefactor.push_back(prefactor);
        m_gamma.push_back(gamma);
        m_epivot.push_back(epivot);
        m_flux.push_back(flux);
        m_eflux.push_back(eflux);
    
    } // endfor: looped over all nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set MC pre-computation cache
 *
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * This method sets up an array of indices and the cumulative distribution
 * function needed for MC simulations.
 ***************************************************************************/
void GModelSpectralFunc::mc_update(const GEnergy& emin,
                                   const GEnergy& emax) const
{
    // Check if we need to update the cache
    if (emin != m_mc_emin || emax != m_mc_emax) {
    
        // Store new energy interval
        m_mc_emin = emin;
        m_mc_emax = emax;
        
        // Initialise cache
        m_mc_cum.clear();
        m_mc_min.clear();
        m_mc_max.clear();
        m_mc_exp.clear();

        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();
    
        // Continue only if e_max > e_min
        if (e_max > e_min) {
        
            // Allocate flux
            double flux;
    
            // Determine left node index for minimum energy
            m_lin_nodes.set_value(e_min);
            int inx_emin = m_lin_nodes.inx_left();

            // Determine left node index for maximum energy
            m_lin_nodes.set_value(e_max);
            int inx_emax = m_lin_nodes.inx_left();
    
            // If both energies are within the same node then just
            // add this one node on the stack
            if (inx_emin == inx_emax) {
                flux = m_prefactor[inx_emin] * 
                       plaw_photon_flux(e_min,
                                        e_max, 
                                        m_epivot[inx_emin],
                                        m_gamma[inx_emin]);
                m_mc_cum.push_back(flux);
                m_mc_min.push_back(e_min);
                m_mc_max.push_back(e_max);
                m_mc_exp.push_back(m_gamma[inx_emin]);
            }
    
            // ... otherwise integrate over the nodes where emin and emax
            // resides and all the remaining nodes
            else {

                // If we are requested to extrapolate beyond first node,
                // the use the first nodes lower energy and upper energy
                // boundary for the initial integration.
                int i_start = (e_min < m_lin_nodes[0]) ? inx_emin : inx_emin+1;

                // Add emin to the node boundary
                flux = m_prefactor[inx_emin] *
                       plaw_photon_flux(e_min,
                                        m_lin_nodes[i_start],
                                        m_epivot[inx_emin],
                                        m_gamma[inx_emin]);
                m_mc_cum.push_back(flux);
                m_mc_min.push_back(e_min);
                m_mc_max.push_back(m_lin_nodes[i_start]);
                m_mc_exp.push_back(m_gamma[inx_emin]);

                // Add all nodes between
                for (int i = i_start; i < inx_emax; ++i) {
                    flux = m_flux[i];
                    m_mc_cum.push_back(flux);
                    m_mc_min.push_back(m_lin_nodes[i]);
                    m_mc_max.push_back(m_lin_nodes[i+1]);
                    m_mc_exp.push_back(m_gamma[i]);
                }

                // Add node boundary to emax
                flux = m_prefactor[inx_emax] *
                       plaw_photon_flux(m_lin_nodes[inx_emax],
                                        e_max,
                                        m_epivot[inx_emax],
                                        m_gamma[inx_emax]);
                m_mc_cum.push_back(flux);
                m_mc_min.push_back(m_lin_nodes[inx_emax]);
                m_mc_max.push_back(e_max);
                m_mc_exp.push_back(m_gamma[inx_emax]);
        
            } // endelse: emin and emax not between same nodes

            // Build cumulative distribution
            for (int i = 1; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] += m_mc_cum[i-1];
            }
            double norm = m_mc_cum[m_mc_cum.size()-1];
            for (int i = 0; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] /= norm;
            }

            // Set MC values
            for (int i = 0; i < m_mc_cum.size(); ++i) {

                // Compute exponent
                double exponent = m_mc_exp[i] + 1.0;
                
                // Exponent dependend computation
                if (std::abs(exponent) > 1.0e-11) {
                    m_mc_exp[i] = exponent;
                    m_mc_min[i] = std::pow(m_mc_min[i], exponent);
                    m_mc_max[i] = std::pow(m_mc_max[i], exponent);
                }
                else {
                    m_mc_exp[i] = 0.0;
                    m_mc_min[i] = std::log(m_mc_min[i]);
                    m_mc_max[i] = std::log(m_mc_max[i]);
                }
            
            } // endfor: set MC values
            
        } // endif: e_max > e_min
        
    } // endif: Update was required


    // Return
    return;
}
