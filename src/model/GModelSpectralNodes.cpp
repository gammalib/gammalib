/***************************************************************************
 *          GModelSpectralNodes.cpp - Spectral nodes model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralNodes.cpp
 * @brief Spectral nodes model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralNodes    g_spectral_nodes_seed;
const GModelSpectralRegistry g_spectral_nodes_registry(&g_spectral_nodes_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                "GModelSpectralNodes::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralNodes::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralNodes::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                      "GModelSpectralNodes::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralNodes::write(GXmlElement&)"
#define G_INSERT       "GModelSpectralNodes::insert(int&, GEnergy&, double&)"
#define G_REMOVE                          "GModelSpectralNodes::remove(int&)"
#define G_ENERGY_GET                      "GModelSpectralNodes::energy(int&)"
#define G_ENERGY_SET            "GModelSpectralNodes::energy(int&, GEnergy&)"
#define G_INTENSITY_GET                "GModelSpectralNodes::intensity(int&)"
#define G_INTENSITY_SET       "GModelSpectralNodes::intensity(int&, double&)"

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
GModelSpectralNodes::GModelSpectralNodes(void) : GModelSpectral()
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
 * Construct spectral nodes model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralNodes::GModelSpectralNodes(const GXmlElement& xml) :
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
 * @param[in] model Spectral nodes model.
 ***************************************************************************/
GModelSpectralNodes::GModelSpectralNodes(const GModelSpectralNodes& model) :
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
GModelSpectralNodes::~GModelSpectralNodes(void)
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
 * @param[in] model Spectral nodes model.
 * @return Spectral nodes model.
 ***************************************************************************/
GModelSpectralNodes& GModelSpectralNodes::operator=(const GModelSpectralNodes& model)
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
 * @brief Clear spectral nodes model
 ***************************************************************************/
void GModelSpectralNodes::clear(void)
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
 * @brief Clone spectral nodes model
***************************************************************************/
GModelSpectralNodes* GModelSpectralNodes::clone(void) const
{
    // Clone spectral nodes model
    return new GModelSpectralNodes(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Computes
 *
 * \f[
 *    S_{\rm E}(E | t) =
 *    10^{(\log {\tt m\_values[i]}) w_{i} + 
 *        (\log {\tt m\_values[i+1]}) w_{i+1}}
 * \f]
 *
 * where
 * - \f${\tt m\_values[i]}\f$ is the intensity of node \f$i\f$,
 * - \f${\tt m\_values[i+1]}\f$ is the intensity of node \f$i+1\f$,
 * - \f$w_{i}\f$ is the weighting of node \f$i\f$,
 * - \f$w_{i+1}\f$ is the weighting of node \f$i+1\f$, and
 * - \f${\tt m\_energies[i]} <= E <= {\tt m\_energies[i+1]}\f$.
 *
 * The weightings \f$w_{i}\f$ and \f$w_{i+1}\f$ are computed by linear
 * interpolation (in the log-log plane) between the nodes
 * \f$(\log {\tt m\_energies[i]}, \log{\tt m\_values[i]})\f$
 * and
 * \f$(\log {\tt m\_energies[i+1]}, \log{\tt m\_values[i+1]})\f$
 * to the requested energy \f$\log E\f$.
 ***************************************************************************/
double GModelSpectralNodes::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime) const
{
    // Update evaluation cache
    update_eval_cache();

    // Interpolate function. This is done in log10-log10 space, but the
    // linear value is returned.
    double exponent = m_log_energies.interpolate(srcEng.log10MeV(), m_log_values);
    double value    = std::pow(10.0, exponent);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralNodes::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", exponent=" << exponent;
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
 * Computes
 *
 * \f[
 *    S_{\rm E}(E | t) =
 *    10^{(\log {\tt m\_values[i]}) w_{i} + 
 *        (\log {\tt m\_values[i+1]}) w_{i+1}}
 * \f]
 *
 * where
 * - \f${\tt m\_values[i]}\f$ is the intensity of node \f$i\f$,
 * - \f${\tt m\_values[i+1]}\f$ is the intensity of node \f$i+1\f$,
 * - \f$w_{i}\f$ is the weighting of node \f$i\f$,
 * - \f$w_{i+1}\f$ is the weighting of node \f$i+1\f$, and
 * - \f${\tt m\_energies[i]} <= E <= {\tt m\_energies[i+1]}\f$.
 *
 * The weightings \f$w_{i}\f$ and \f$w_{i+1}\f$ are computed by linear
 * interpolation (in the log-log plane) between the nodes
 * \f$(\log {\tt m\_energies[i]}, \log{\tt m\_values[i]})\f$
 * and
 * \f$(\log {\tt m\_energies[i+1]}, \log{\tt m\_values[i+1]})\f$
 * to the requested energy \f$\log E\f$.
 *
 * The method also evaluates the partial derivatives of the model with
 * respect to the parameters using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_values[i]}} =
 *      \frac{S_{\rm E}(E | t) \, w_i}{{\tt m\_values[i]}}
 * \f]
 ***************************************************************************/
double GModelSpectralNodes::eval_gradients(const GEnergy& srcEng,
                                           const GTime&   srcTime)
{
    // Update evaluation cache
    update_eval_cache();

    // Set interpolation value
    m_log_energies.set_value(srcEng.log10MeV());

    // Get indices and weights for interpolation
    int    inx_left  = m_log_energies.inx_left();
    int    inx_right = m_log_energies.inx_right();
    double wgt_left  = m_log_energies.wgt_left();
    double wgt_right = m_log_energies.wgt_right();

    // Interpolate function
    double exponent = m_log_values[inx_left]  * wgt_left +
                      m_log_values[inx_right] * wgt_right;

    // Compute linear value
    double value = std::pow(10.0, exponent);

    // Initialise gradients
    for (int i = 0; i < m_values.size(); ++i) {
        m_values[i].factor_gradient(0.0);
    }

    // Gradient for left node
    if (m_values[inx_left].isfree()) {
        double grad = value * wgt_left / m_values[inx_left].factor_value();
        m_values[inx_left].factor_gradient(grad);
    }

    // Gradient for right node
    if (m_values[inx_right].isfree()) {
        double grad = value * wgt_right / m_values[inx_right].factor_value();
        m_values[inx_right].factor_gradient(grad);
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralNodes::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", exponent=" << exponent;
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
double GModelSpectralNodes::flux(const GEnergy& emin, const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Update flux cache
        update_flux_cache();

        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();
    
        // Determine left node index for minimum energy
        m_lin_energies.set_value(e_min);
        int inx_emin = m_lin_energies.inx_left();

        // Determine left node index for maximum energy
        m_lin_energies.set_value(e_max);
        int inx_emax = m_lin_energies.inx_left();
    
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
            int i_start = (e_min < m_lin_energies[0]) ? inx_emin : inx_emin+1;

            // Integrate from emin to the node boundary
            flux = m_prefactor[inx_emin] *
                   plaw_photon_flux(e_min,
                                    m_lin_energies[i_start],
                                    m_epivot[inx_emin],
                                    m_gamma[inx_emin]);

            // Integrate over all nodes between
            for (int i = i_start; i < inx_emax; ++i) {
                flux += m_flux[i];
            }

            // Integrate from node boundary to emax
            flux += m_prefactor[inx_emax] *
                    plaw_photon_flux(m_lin_energies[inx_emax],
                                     e_max,
                                     m_epivot[inx_emax],
                                     m_gamma[inx_emax]);
        }

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
double GModelSpectralNodes::eflux(const GEnergy& emin, const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Update flux cache
        update_flux_cache();
    
        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();
    
        // Determine left node index for minimum energy
        m_lin_energies.set_value(e_min);
        int inx_emin = m_lin_energies.inx_left();

        // Determine left node index for maximum energy
        m_lin_energies.set_value(e_max);
        int inx_emax = m_lin_energies.inx_left();
    
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
            int i_start = (e_min < m_lin_energies[0]) ? inx_emin : inx_emin+1;

            // Integrate from emin to the node boundary
            eflux = m_prefactor[inx_emin] *
                    plaw_energy_flux(e_min,
                                     m_lin_energies[i_start],
                                     m_epivot[inx_emin],
                                     m_gamma[inx_emin]) * gammalib::MeV2erg;

            // Integrate over all nodes between
            for (int i = i_start; i < inx_emax; ++i) {
                eflux += m_eflux[i];
            }

            // Integrate from node boundary to emax
            eflux += m_prefactor[inx_emax] *
                     plaw_energy_flux(m_lin_energies[inx_emax],
                                      e_max,
                                      m_epivot[inx_emax],
                                      m_gamma[inx_emax]) * gammalib::MeV2erg;
        
        } // endelse: emin and emax not between same nodes
    
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
 * Returns Monte Carlo energy by randomly drawing from node function.
 ***************************************************************************/
GEnergy GModelSpectralNodes::mc(const GEnergy& emin,
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
            if (m_mc_cum[inx-1] <= u)
                break;
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
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *     <spectrum type="NodeFunction">
 *       <node>
 *         <parameter name="Energy"    scale=".." value=".." min=".." max=".." free=".."/>
 *         <parameter name="Intensity" scale=".." value=".." min=".." max=".." free=".."/>
 *       </node>
 *       ...
 *       <node>
 *         <parameter name="Energy"    scale=".." value=".." min=".." max=".." free=".."/>
 *         <parameter name="Intensity" scale=".." value=".." min=".." max=".." free=".."/>
 *       </node>
 *     </spectrum>
 *
 * @todo Check that all energies and intensities are positive
 * @todo Check that nodes are ordered
 * @todo Check that energy boundaries are not overlapping
 ***************************************************************************/
void GModelSpectralNodes::read(const GXmlElement& xml)
{
    // Free space for nodes
    m_energies.clear();
    m_values.clear();

    // Get number of nodes from XML file
    int nodes = xml.elements("node");

    // Throw an error if there not at least two nodes
    if (nodes < 1) {
        std::string message = "Node function model requires at least two"
                              " nodes.";
        throw GException::model_invalid_parnum(G_READ, xml, message);
    }

    // Loop over all nodes
    for (int i = 0; i < nodes; ++i) {

        // Allocate node parameters
        GModelPar energy;
        GModelPar intensity;
            
        // Get node
        const GXmlElement* node = xml.element("node", i);

        // Verify that node XML element has exactly 2 parameters
        if (node->elements() != 2 || node->elements("parameter") != 2) {
            throw GException::model_invalid_parnum(G_READ, xml,
                  "Node requires exactly 2 parameters.");
        }

        // Extract node parameters
        int npar[] = {0, 0};
        for (int k = 0; k < 2; ++k) {

            // Get parameter element
            const GXmlElement* par = node->element("parameter", k);

            // Handle energy
            if (par->attribute("name") == "Energy") {
                energy.read(*par);
                npar[0]++;
            }

            // Handle intensity
            else if (par->attribute("name") == "Intensity") {
                intensity.read(*par);
                npar[1]++;
            }

        } // endfor: looped over parameters

        // Verify that all parameters were found
        if (npar[0] != 1 || npar[1] != 1) {
            throw GException::model_invalid_parnames(G_READ, xml,
                  "Require \"Energy\" and \"Intensity\" parameters.");
        }

        // Set parameter names
        std::string energy_name    = "Energy"+str(i);
        std::string intensity_name = "Intensity"+str(i);

        // Set energy attributes
        energy.name(energy_name);
        energy.unit("MeV");
        energy.hasgrad(false);

        // Set intensity attributes
        intensity.name(intensity_name);
        intensity.unit("ph/cm2/s/MeV");
        intensity.hasgrad(true);

        // Push node parameters on list
        m_energies.push_back(energy);
        m_values.push_back(intensity);

    } // endfor: looped over nodes

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of required type
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="NodeFunction">
 *       <node>
 *         <parameter name="Energy"    scale=".." value=".." min=".." max=".." free=".."/>
 *         <parameter name="Intensity" scale=".." value=".." min=".." max=".." free=".."/>
 *       </node>
 *       ...
 *       <node>
 *         <parameter name="Energy"    scale=".." value=".." min=".." max=".." free=".."/>
 *         <parameter name="Intensity" scale=".." value=".." min=".." max=".." free=".."/>
 *       </node>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralNodes::write(GXmlElement& xml) const
{
    // Determine number of nodes
    int nodes = m_energies.size();

    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append nodes
    if (xml.elements() == 0) {
        for (int i = 0; i < nodes; ++i) {
            xml.append(GXmlElement("node"));
        }
    }

    // Verify that XML element has the required number of nodes
    if (xml.elements() != nodes || xml.elements("node") != nodes) {
        std::string message = "Spectral model requires exactly " +
                              str(nodes) + " nodes.";
        throw GException::model_invalid_parnum(G_WRITE, xml, message);
    }

    // Loop over all nodes
    for (int i = 0; i < nodes; ++i) {

        // Get node
        GXmlElement* node = xml.element("node", i);

        // If XML element has 0 leafes then append energy and intensity
        // element
        if (node->elements() == 0) {
            node->append(GXmlElement("parameter name=\"Energy\""));
            node->append(GXmlElement("parameter name=\"Intensity\""));
        }

        // Verify that node XML element has exactly 2 parameters
        if (node->elements() != 2 || node->elements("parameter") != 2) {
            throw GException::model_invalid_parnum(G_WRITE, xml,
                  "Node requires exactly 2 parameters.");
        }

        // Set or update model parameter attributes
        int npar[] = {0, 0};
        for (int k = 0; k < 2; ++k) {

            // Get parameter element
            GXmlElement* par = node->element("parameter", k);

            // Handle prefactor
            if (par->attribute("name") == "Energy") {
                npar[0]++;
                m_energies[i].write(*par);
            }

            // Handle index
            else if (par->attribute("name") == "Intensity") {
                npar[1]++;
                m_values[i].write(*par);
            }

        } // endfor: looped over parameters

        // Check of all required parameters are present
        if (npar[0] != 1 || npar[1] != 1) {
            throw GException::model_invalid_parnames(G_WRITE, xml,
                  "Require \"Energy\" and \"Intensity\" parameters.");
        }

    } // endfor: looped over nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append node
 *
 * @param[in] energy Node energy.
 * @param[in] intensity Node intensity.
 *
 * Appends one node to the node function.
 ***************************************************************************/
void GModelSpectralNodes::append(const GEnergy& energy,
                                 const double&  intensity)
{
    // Allocate node parameters
    GModelPar e_par;
    GModelPar i_par;

    // Set energy attributes
    e_par.name("Energy");
    e_par.value(energy.MeV());
    e_par.unit("MeV");
    e_par.hasgrad(false);

    // Set intensity attributes
    i_par.name("Intensity");
    i_par.value(intensity);
    i_par.unit("ph/cm2/s/MeV");
    i_par.hasgrad(true);

    // Append to nodes
    m_energies.push_back(e_par);
    m_values.push_back(i_par);

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert node
 *
 * @param[in] index Node index [0,...,nodes()-1].
 * @param[in] energy Node energy.
 * @param[in] intensity Node intensity.
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Inserts a node into the node function before the node with the specified
 * @p index.
 ***************************************************************************/
void GModelSpectralNodes::insert(const int&     index,
                                 const GEnergy& energy,
                                 const double&  intensity)
{
    // Raise exception if index is outside boundary
    #if defined(G_RANGE_CHECK)
    if (m_energies.empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, index, 0, nodes()-1);
        }
    }
    else {
        if (index < 0 || index >= nodes()) {
            throw GException::out_of_range(G_INSERT, index, 0, nodes()-1);
        }
    }
    #endif

    // Allocate node parameters
    GModelPar e_par;
    GModelPar i_par;

    // Set energy attributes
    e_par.name("Energy");
    e_par.value(energy.MeV());
    e_par.unit("MeV");
    e_par.hasgrad(false);

    // Set intensity attributes
    i_par.name("Intensity");
    i_par.value(intensity);
    i_par.unit("ph/cm2/s/MeV");
    i_par.hasgrad(true);

    // Insert node
    m_energies.insert(m_energies.begin()+index, e_par);
    m_values.insert(m_values.begin()+index, i_par);

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove node
 *
 * @param[in] index Node index [0,...,nodes()-1].
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Removes node of specified @p index from the node function.
 ***************************************************************************/
void GModelSpectralNodes::remove(const int& index)
{
    // Raise exception if index is outside boundary
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= nodes()) {
        throw GException::out_of_range(G_REMOVE, index, 0, nodes()-1);
    }
    #endif

    // Erase energy and intensity
    m_energies.erase(m_energies.begin() + index);
    m_values.erase(m_values.begin() + index);

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve space for nodes
 *
 * @param[in] num Number of reseved nodes.
 *
 * Reserves space for @p num nodes.
 ***************************************************************************/
void GModelSpectralNodes::reserve(const int& num)
{
    // Reserve space
    m_energies.reserve(num);
    m_values.reserve(num);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append nodes from node function
 *
 * @param[in] nodes Node function.
 *
 * Appends all nodes from a node function to current object.
 ***************************************************************************/
void GModelSpectralNodes::extend(const GModelSpectralNodes& nodes)
{
    // Get number of nodes in node function. Note that we extract the size
    // first to avoid an endless loop that arises when a container is
    // appended to itself.
    int num = nodes.nodes();

    // Continue only if node function is not empty
    if (num > 0) {

        // Reserve enough space
        reserve(this->nodes() + num);

        // Loop over all nodes and append them to the node function
        for (int i = 0; i < num; ++i) {
            m_energies.push_back(nodes.m_energies[i]);
            m_values.push_back(nodes.m_values[i]);
        }

        // Update parameter mapping
        update_pars();

        // Set pre-computation cache
        set_cache();

    } // endif: node function was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return node energy
 *
 * @param[in] index Node index [0,...,nodes()-1].
 * @return Energy of node @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns the energy of node @p index.
 ***************************************************************************/
GEnergy GModelSpectralNodes::energy(const int& index) const
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= nodes()) {
        throw GException::out_of_range(G_ENERGY_GET, index, nodes()-1);
    }
    #endif

    // Retrieve energy
    GEnergy energy;
    energy.MeV(m_energies[index].value());

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Set node energy
 *
 * @param[in] index Node index [0,...,nodes()-1].
 * @param[in] energy Node energy.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Sets the energy of node @p index.
 ***************************************************************************/
void GModelSpectralNodes::energy(const int& index, const GEnergy& energy)
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= nodes()) {
        throw GException::out_of_range(G_ENERGY_SET, index, nodes()-1);
    }
    #endif

    // Set energy
    m_energies[index].value(energy.MeV());

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return node intensity
 *
 * @param[in] index Node index [0,...,nodes()-1].
 * @return Intensity of node @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns the intensity of node @p index.
 ***************************************************************************/
double GModelSpectralNodes::intensity(const int& index) const
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= nodes()) {
        throw GException::out_of_range(G_INTENSITY_GET, index, nodes()-1);
    }
    #endif

    // Return intensity
    return (m_values[index].value());
}


/***********************************************************************//**
 * @brief Set node intensity
 *
 * @param[in] index Node index [0,...,nodes()-1].
 * @param[in] intensity Node Intensity.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Set the intensity of node @p index.
 ***************************************************************************/
void GModelSpectralNodes::intensity(const int& index, const double& intensity)
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= nodes()) {
        throw GException::out_of_range(G_INTENSITY_SET, index, nodes()-1);
    }
    #endif

    // Set intensity
    m_values[index].value(intensity);

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print node function information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing node function information.
 ***************************************************************************/
std::string GModelSpectralNodes::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralNodes ===");

        // Append information
        result.append("\n"+parformat("Number of nodes")+str(m_energies.size()));
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
void GModelSpectralNodes::init_members(void)
{
    // Initialise node energies and values
    m_energies.clear();
    m_values.clear();

    // Initialise evaluation cache
    m_old_energies.clear();
    m_old_values.clear();
    m_log_energies.clear();
    m_log_values.clear();

    // Initialise flux computation cache
    m_lin_energies.clear();
    m_lin_values.clear();
    m_prefactor.clear();
    m_gamma.clear();
    m_epivot.clear();
    m_flux.clear();
    m_eflux.clear();
    
    // Initialise MC cache
    m_mc_emin.clear();
    m_mc_emax.clear();
    m_mc_cum.clear();
    m_mc_min.clear();
    m_mc_max.clear();
    m_mc_exp.clear();

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral function model.
 ***************************************************************************/
void GModelSpectralNodes::copy_members(const GModelSpectralNodes& model)
{
    // Copy node energies and values
    m_energies     = model.m_energies;
    m_values       = model.m_values;
    
    // Copy evaluation cache
    m_old_energies = model.m_old_energies;
    m_old_values   = model.m_old_values;
    m_log_energies = model.m_log_energies;
    m_log_values   = model.m_log_values;

    // Copy flux computation cache
    m_lin_energies = model.m_lin_energies;
    m_lin_values   = model.m_lin_values;
    m_prefactor    = model.m_prefactor;
    m_gamma        = model.m_gamma;
    m_epivot       = model.m_epivot;
    m_flux         = model.m_flux;
    m_eflux        = model.m_eflux;

    // Copy MC cache
    m_mc_emin      = model.m_mc_emin;
    m_mc_emax      = model.m_mc_emax;
    m_mc_cum       = model.m_mc_cum;
    m_mc_min       = model.m_mc_min;
    m_mc_max       = model.m_mc_max;
    m_mc_exp       = model.m_mc_exp;

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralNodes::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update parameter mapping
 *
 * Sets the parameter pointers in the m_pars array, enabling iterating over
 * all model parameters. This method needs to be called after changing the
 * number of nodes in the spectral model. The method needs not to be called
 * after value update.
 ***************************************************************************/
void GModelSpectralNodes::update_pars(void)
{
    // Clear parameter pointer(s)
    m_pars.clear();

    // Get number of nodes
    int nodes = m_energies.size();

    // Set parameter pointers for all nodes
    for (int i = 0; i < nodes; ++i) {

        // Set parameter names
        std::string energy_name    = "Energy"+str(i);
        std::string intensity_name = "Intensity"+str(i);

        // Set energy attributes
        m_energies[i].name(energy_name);
        m_energies[i].unit("MeV");
        m_energies[i].hasgrad(false);

        // Set intensity attributes
        m_values[i].name(intensity_name);
        m_values[i].unit("ph/cm2/s/MeV");
        m_values[i].hasgrad(true);

        // Set pointer
        m_pars.push_back(&(m_energies[i]));
        m_pars.push_back(&(m_values[i]));

    } // endfor: looped over nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pre-computation cache
 ***************************************************************************/
void GModelSpectralNodes::set_cache(void) const
{
    // Set evaluation cache
    set_eval_cache();

    // Set flux computation cache
    set_flux_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set evaluation cache
 *
 * The evaluation cache contains pre-computed results that are needed for
 * fast function evaluation.
 *
 * @todo Check that all energies and intensities are > 0
 ***************************************************************************/
void GModelSpectralNodes::set_eval_cache(void) const
{
    // Clear any existing values
    m_old_energies.clear();
    m_old_values.clear();
    m_log_energies.clear();
    m_log_values.clear();

    // Compute log10 of energies and intensities
    for (int i = 0; i < m_energies.size(); ++i) {
    
        // Set log10(energy)
        double log10energy = std::log10(m_energies[i].value());
        m_old_energies.push_back(log10energy);
        m_log_energies.append(log10energy);

        // Set log10(intensity)
        double log10value = std::log10(m_values[i].value());
        m_old_values.push_back(log10value);
        m_log_values.push_back(log10value);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set flux computation cache
 *
 * Pre-computes some values that are needed for flux computation.
 *
 * @todo Handle special case emin=emax and fmin=fmax
 ***************************************************************************/
void GModelSpectralNodes::set_flux_cache(void) const
{
    // Clear any existing values
    m_lin_energies.clear();
    m_lin_values.clear();
    m_prefactor.clear();
    m_gamma.clear();
    m_epivot.clear();
    m_flux.clear();
    m_eflux.clear();

    // Store linear energies and values
    for (int i = 0; i < m_energies.size(); ++i) {
        m_lin_energies.append(m_energies[i].value());
        m_lin_values.push_back(m_values[i].value());
    }

    // Loop over all nodes-1
    for (int i = 0; i < m_energies.size()-1; ++i) {
    
        // Get energies and function values
        double emin = m_lin_energies[i];
        double emax = m_lin_energies[i+1];
        double fmin = m_values[i].value();
        double fmax = m_values[i+1].value();
    
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
 * @brief Update evaluation cache
 *
 * Updates the evaluation cache by computing only results for values that
 * changed.
 *
 * @todo Check that all energies and intensities are > 0
 ***************************************************************************/
void GModelSpectralNodes::update_eval_cache(void) const
{
    // Update energies
    for (int i = 0; i < m_energies.size(); ++i) {
        double energy = m_energies[i].value();
        if (energy != m_old_energies[i]) {
            m_log_energies[i] = std::log10(energy);
            m_old_energies[i] = energy;
        }
    }

    // Update intensities
    for (int i = 0; i < m_values.size(); ++i) {
        double value = m_values[i].value();
        if (value != m_old_values[i]) {
            m_log_values[i] = std::log10(value);
            m_old_values[i] = value;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update flux computation cache
 *
 * Updates the flux computation cache if either the energy boundaries or the
 * intensity values have changed.
 *
 * @todo Handle special case emin=emax and fmin=fmax
 ***************************************************************************/
void GModelSpectralNodes::update_flux_cache(void) const
{
    // Loop over all nodes-1
    for (int i = 0; i < m_energies.size()-1; ++i) {
    
        // Get energies and function values
        double emin = m_lin_energies[i];
        double emax = m_lin_energies[i+1];
        double fmin = m_values[i].value();
        double fmax = m_values[i+1].value();

        // Update values only if energies or function values have changed
        if (emin != m_energies[i].value() ||
            emax != m_energies[i+1].value() ||
            fmin != m_values[i].value() ||
            fmax != m_values[i+1].value()) {
    
            // Compute pivot energy (MeV). We use here the geometric mean
            // of the node boundaries.
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
        
            // Store values on pre-computation cache
            m_prefactor[i] = prefactor;
            m_gamma[i]     = gamma;
            m_epivot[i]    = epivot;
            m_flux[i]      = flux;
            m_eflux[i]     = eflux;

        } // endif: update was required
    
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
void GModelSpectralNodes::mc_update(const GEnergy& emin, const GEnergy& emax) const
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
            m_lin_energies.set_value(e_min);
            int inx_emin = m_lin_energies.inx_left();

            // Determine left node index for maximum energy
            m_lin_energies.set_value(e_max);
            int inx_emax = m_lin_energies.inx_left();
    
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
                int i_start = (e_min < m_lin_energies[0]) ? inx_emin : inx_emin+1;

                // Add emin to the node boundary
                flux = m_prefactor[inx_emin] *
                       plaw_photon_flux(e_min,
                                        m_lin_energies[i_start],
                                        m_epivot[inx_emin],
                                        m_gamma[inx_emin]);
                m_mc_cum.push_back(flux);
                m_mc_min.push_back(e_min);
                m_mc_max.push_back(m_lin_energies[i_start]);
                m_mc_exp.push_back(m_gamma[inx_emin]);

                // Add all nodes between
                for (int i = i_start; i < inx_emax; ++i) {
                    flux = m_flux[i];
                    m_mc_cum.push_back(flux);
                    m_mc_min.push_back(m_lin_energies[i]);
                    m_mc_max.push_back(m_lin_energies[i+1]);
                    m_mc_exp.push_back(m_gamma[i]);
                }

                // Add node boundary to emax
                flux = m_prefactor[inx_emax] *
                       plaw_photon_flux(m_lin_energies[inx_emax],
                                        e_max,
                                        m_epivot[inx_emax],
                                        m_gamma[inx_emax]);
                m_mc_cum.push_back(flux);
                m_mc_min.push_back(m_lin_energies[inx_emax]);
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
