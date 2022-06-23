/***************************************************************************
 *          GModelSpectralTable.cpp - Spectral table model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTable.cpp
 * @brief Spectral table model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GEnergy.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GModelSpectralTable.hpp"
#include "GModelSpectralTablePar.hpp"
#include "GModelSpectralTablePars.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralTable    g_spectral_table_seed;
const GModelSpectralRegistry g_spectral_table_registry(&g_spectral_table_seed);

/* __ Method name definitions ____________________________________________ */
#define G_CONST   "GModelSpectralTable(GEbounds&, GModelSpectralTablePars&, "\
                                                                 "GNdarray&)"
#define G_FLUX                "GModelSpectralTable::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralTable::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralTable::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                      "GModelSpectralTable::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralTable::write(GXmlElement&)"
#define G_LOAD                        "GModelSpectralTable::load(GFilename&)"
#define G_TABLE_PAR                    "GModelSpectralTable::table_par(int&)"
#define G_ENERGY_SCALE      "GModelSpectralTable::energy_scale(std::string&)"
#define G_LOAD_PAR                    "GModelSpectralTable::load_par(GFits&)"
#define G_PAR_INDEX            "GModelSpectralTable::par_index(std::string&)"
#define G_UPDATE                              "GModelSpectralTable::update()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_LOAD_SPEC                     //!< Debug load_spec() method
//#define G_DEBUG_ESCALE                            //!< Debug energy scaling
//#define G_DEBUG_UPDATE                           //!< Debug update() method
//#define G_DEBUG_CALLTREE                               //!< Debug call tree


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name of table model.
 * @param[in] norm Normalization factor.
 *
 * Constructs a spectral table model from a FITS file. See the load() method
 * for more information about the expected structure of the FITS file.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GFilename& filename,
                                         const double&    norm) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Load table
    load(filename);

    // Set normalization as the scale factor and limit change in
    // normalisation to a factor of 1000
    m_norm.remove_range();
    m_norm.scale(norm);
    m_norm.value(norm);
    m_norm.range(0.0,1.0e3*norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Table model constructor
 *
 * @param[in] ebounds Energy boundaries.
 * @param[in] pars Table model parameters.
 * @param[in] spectra Spectra.
 *
 * Constructs a spectral table model from energy boundaries, table model
 * parameters, and spectra.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GEbounds&                ebounds,
                                         const GModelSpectralTablePars& pars,
                                         const GNdarray&                spectra) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Check consistency of spectra dimension
    if (spectra.dim() != pars.size()+1) {
        std::string msg = "Spectra array has "+gammalib::str(spectra.dim())+
                          " dimensions but an array with "+
                          gammalib::str(pars.size()+1)+" dimensions is "
                          "expected. Please specify a spectra array with the "
                          "correct dimension.";
        throw GException::invalid_argument(G_CONST, msg);
    }
    
    // Check consistency of parameter dimensions
    for (int i = 0; i < pars.size(); ++i) {
        int npars = pars[i]->size();
        int nspec = spectra.shape()[i];
        if (npars != nspec) {
            std::string msg = "Parameter \""+pars[i]->par().name()+"\" has "+
                              gammalib::str(npars)+" values but there are "+
                              gammalib::str(nspec)+" spectra for table axis "+
                              gammalib::str(i)+". Please specify a spectra "
                              "array with the correct number of spectra.";
            throw GException::invalid_argument(G_CONST, msg);
        }
    }
 
    // Check consistency of energy bins
    int nebins = ebounds.size();
    int nspec  = spectra.shape()[pars.size()];
    if (nebins != nspec) {
        std::string msg = "Spectra array has "+gammalib::str(nspec)+" energy "
                          "bins but there are "+gammalib::str(nebins)+
                          " energy boundaries. Please specify a spectra "
                          "array with the correct number of energy bins.";
        throw GException::invalid_argument(G_CONST, msg);
    }

    // Set members
    m_ebounds    = ebounds;
    m_table_pars = pars;
    m_spectra    = spectra;

    // Set energy nodes
    set_energy_nodes();

    // Set parameter pointers
    set_par_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a spectral table model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GXmlElement& xml) :
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
 * @param[in] model Table model.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GModelSpectralTable& model) :
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
GModelSpectralTable::~GModelSpectralTable(void)
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
 * @param[in] model Table model.
 * @return Table model.
 ***************************************************************************/
GModelSpectralTable& GModelSpectralTable::operator=(const GModelSpectralTable& model)
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
 * @brief Clear table model
***************************************************************************/
void GModelSpectralTable::clear(void)
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
 * @brief Clone table model
***************************************************************************/
GModelSpectralTable* GModelSpectralTable::clone(void) const
{
    // Clone table model
    return new GModelSpectralTable(*this);
}


/***********************************************************************//**
 * @brief Evaluate spectral table model
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm} \times
 *                       \left( w_l F_l(p) + w_r F_r(p) \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization factor of the spectral table
 *   model,
 * - \f$w_l\f$ is the weight of the spectral vector \f$F_l(p)\f$ with the
 *   largest energy \f$E_l\f$ that satisfies \f$E<E_l\f$, and
 * - \f$w_r\f$ is the weight of the spectral vector \f$F_r(p)\f$ with the
 *   smallest energy \f$E_r\f$ that satisfies \f$E>E_r\f$.
 *
 * The weights are computed using
 *
 * \f[
 *    w_r = \frac{\log_{10} E - \log_{10} E_l}{\log_{10} E_r - \log_{10} E_l}
 * \f]
 *
 * and \f$w_l = 1 - w_r\f$.
 *
 * If @p gradient is true, the method also computes the parameter gradients
 * using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_norm}} =
 *      \frac{S_{\rm E}(E | t)}{{\tt m\_norm}}
 * \f]
 *
 * and
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta p} =
 *    {\tt m\_norm} \times
 *    \left( w_l \frac{\delta F_l(p)}{\delta p} +
 *           w_r \frac{\delta F_r(p)}{\delta p} \right)
 * \f]
 *
 * for all other parameters.
 *
 * For the computation of \f$F_l(p)\f$, \f$F_r(p)\f$,
 * \f$\frac{\delta F_l(p)}{\delta p}\f$, and
 * \f$\frac{\delta F_r(p)}{\delta p}\f$ see the update() method.
 ***************************************************************************/
double GModelSpectralTable::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime,
                                 const bool&    gradients) const
{
    // Update spectral function
    update();

    // Get energy to interpolate. In case that the spectrum was scaled
    // then use now to actual parameter value of the energy scale to
    // get scaling factor
    double log10energy = srcEng.log10MeV();
    if (m_log10escale != 0.0) {
        log10energy -= std::log10(table_par(m_escale_par).par().value());
        #if defined(G_DEBUG_ESCALE)
        std::cout << "GModelSpectralTable::eval";
        std::cout << " " << srcEng.log10MeV();
        std::cout << " " << log10energy << std::endl;
        #endif
    }

    // Interpolate function in log10 energy
    m_log_nodes.set_value(log10energy);
    double wgt_left  = m_log_nodes.wgt_left();
    double wgt_right = m_log_nodes.wgt_right();
    int    inx_left  = m_log_nodes.inx_left();
    int    inx_right = m_log_nodes.inx_right();
    double val_left  = m_lin_values(inx_left, 0);
    double val_right = m_lin_values(inx_right,0);
    double func      = wgt_left  * val_left + wgt_right * val_right;

    // Compute function value
    double value = m_norm.value() * func;

    // Optionally compute gradients
    if (gradients) {

        // Compute partial derivative of function normalisation
        double g_norm  = (m_norm.is_free())  ? m_norm.scale() * func : 0.0;

        // Set normalisation gradient
        m_norm.factor_gradient(g_norm);

        // Compute partial derivatives of all other free parameters
        int dim = m_table_pars.size();
        for (int i = 0; i < dim; ++i) {

            // Initialise gradient
            double grad = 0.0;

            // Get reference to model parameter (circumvent const correctness)
            GModelPar& par =
                const_cast<GModelSpectralTablePar*>(m_table_pars[i])->par();

            // If parameter is free then compute gradient
            if (par.is_free()) {
                double grad_left  = m_lin_values(inx_left, i+1);
                double grad_right = m_lin_values(inx_right,i+1);
                if ((m_log10escale != 0.0) && (par.name() == m_escale_par)) {
                    double dxdp = -1.0 / (gammalib::ln10 * par.factor_value());
                    grad        = ((m_log_nodes.wgt_grad_left()  * dxdp * val_left)  +
                                   (m_log_nodes.wgt_grad_right() * dxdp * val_right) +
                                   (wgt_left  * grad_left + wgt_right * grad_right) *
                                    par.scale()) *
                                   m_norm.value();
                }
                else {
                    grad = (wgt_left  * grad_left + wgt_right * grad_right) *
                           par.scale() * m_norm.value();
                }
            }

            // Set gradient
            par.factor_gradient(grad);

        } // endfor: looped over all parameters

    } // endif: gradient computation was requested

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralTable::eval";
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
double GModelSpectralTable::flux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // Update spectral function and flux cache
        update();
        update_flux();

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
                   gammalib::plaw_photon_flux(e_min,
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
                   gammalib::plaw_photon_flux(e_min,
                                              m_lin_nodes[i_start],
                                              m_epivot[inx_emin],
                                              m_gamma[inx_emin]);

            // Integrate over all nodes between
            for (int i = i_start; i < inx_emax; ++i) {
                flux += m_flux[i];
            }

            // Integrate from node boundary to emax
            flux += m_prefactor[inx_emax] *
                    gammalib::plaw_photon_flux(m_lin_nodes[inx_emax],
                                               e_max,
                                               m_epivot[inx_emax],
                                               m_gamma[inx_emax]);

        } // endelse: emin and emax not between same nodes

        // Multiply flux by normalisation factor
        flux *= norm();

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
double GModelSpectralTable::eflux(const GEnergy& emin,
                                  const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // Update spectral function and flux cache
        update();
        update_flux();

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
                    gammalib::plaw_energy_flux(e_min,
                                               e_max, 
                                               m_epivot[inx_emin],
                                               m_gamma[inx_emin]) *
                                               gammalib::MeV2erg;
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
                    gammalib::plaw_energy_flux(e_min,
                                               m_lin_nodes[i_start],
                                               m_epivot[inx_emin],
                                               m_gamma[inx_emin]) *
                                               gammalib::MeV2erg;

            // Integrate over all nodes between
            for (int i = i_start; i < inx_emax; ++i) {
                eflux += m_eflux[i];
            }

            // Integrate from node boundary to emax
            eflux += m_prefactor[inx_emax] *
                     gammalib::plaw_energy_flux(m_lin_nodes[inx_emax],
                                                e_max,
                                                m_epivot[inx_emax],
                                                m_gamma[inx_emax]) *
                                                gammalib::MeV2erg;

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
 * Returns Monte Carlo energy by randomly drawing from a broken power law
 * defined by the file function.
 ***************************************************************************/
GEnergy GModelSpectralTable::mc(const GEnergy& emin,
                                const GEnergy& emax,
                                const GTime&   time,
                                GRan&          ran) const
{
    // Check energy interval
    gammalib::check_energy_interval(G_MC, emin, emax);

    // Allocate energy
    GEnergy energy;

    // Update spectral function and flux cache
    update();
    update_flux();
    update_mc(emin, emax);

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
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *     <spectrum type="TableModel" file="..">
 *       <parameter name="Normalization" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 *
 * Optionally, values for the model table parameters can also be provided.
 ***************************************************************************/
void GModelSpectralTable::read(const GXmlElement& xml)
{
    // Load table model from file (do this first since method calls clear())
    load(gammalib::xml_file_expand(xml, xml.attribute("file")));

    // Get normalisation parameter pointer
    const GXmlElement* norm = gammalib::xml_get_par(G_READ, xml, m_norm.name());

    // Read normalisation parameter
    m_norm.read(*norm);

    // Optionally read all table parameters
    for (int i = 0; i < m_table_pars.size(); ++i) {
        GModelPar& par = m_table_pars[i]->par();
        if (gammalib::xml_has_par(xml, par.name())) {
            const GXmlElement* element =
                               gammalib::xml_get_par(G_READ, xml, par.name());
            par.read(*element);
        }
    }

    // Optionally set energy scale
    if (xml.has_attribute("energyscale")) {
        energy_scale(xml.attribute("energyscale"));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="FileFunction" file="..">
 *       <parameter name="Normalization" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 * 
 * In addition, the method writes the model table parameters into the XML
 * file.
 *
 * Note that the function nodes will not be written since they will not be
 * altered by any method.
 ***************************************************************************/
void GModelSpectralTable::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get normalisation parameter
    GXmlElement* norm  = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());

    // Write normalisation parameter
    m_norm.write(*norm);

    // Write all table parameters
    for (int i = 0; i < m_table_pars.size(); ++i) {
        const GModelPar& par     = m_table_pars[i]->par();
        GXmlElement*     element = gammalib::xml_need_par(G_WRITE, xml, par.name());
        par.write(*element);
    }

    // Set file attribute
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Optionally write energy scale attribute
    if (!m_escale_par.empty()) {
        xml.attribute("energyscale", m_escale_par);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load table from file
 *
 * @param[in] filename File name.
 *
 * Loads table model from FITS file. The format of the FITS file complies
 * with https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html
 ***************************************************************************/
void GModelSpectralTable::load(const GFilename& filename)
{
    // Clear table model
    clear();

    // Set filename
    m_filename = filename;

    // Continue only if filename is not empty
    if (!filename.is_empty()) {

        // Open FITS file
        GFits fits(filename);

        // Load data from extensions
        load_par(fits);
        load_eng(fits);
        load_spec(fits);

        // Close FITS file
        fits.close();

        // Set energy nodes
        set_energy_nodes();

        // Set parameter pointers
        set_par_pointers();

    } // endif: filename was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table into file
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file?
 *
 * Save the table model into a FITS file. The FITS file will contain three
 * binary table extensions:
 *
 *       * PARAMETERS - Table model parameters
 *       * ENERGIES - Table model energies
 *       * SPECTRA - Table model spectra
 *
 * The format of the FITS file complies with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html
 ***************************************************************************/
void GModelSpectralTable::save(const GFilename& filename,
                               const bool&      clobber) const
{
    // Create FITS file
    GFits fits;

    // Create binary tables
    GFitsBinTable par_table  = create_par_table();
    GFitsBinTable eng_table  = create_eng_table();
    GFitsBinTable spec_table = create_spec_table();

    // Append binary tables
    fits.append(par_table);
    fits.append(eng_table);
    fits.append(spec_table);

    // Set keywords in primary header
    GFitsHDU* primary = fits[0];
    primary->card("CONTENT",  "MODEL", "Spectrum file");
    primary->card("FILENAME", filename.url(), "FITS file name");
    primary->card("ORIGIN", PACKAGE_NAME, "Origin of FITS file");
    primary->card("MODLNAME", "model", "Model name");
    primary->card("MODLUNIT", "photons/cm^2/s/MeV", "Model units");
    primary->card("REDSHIFT", false, "If true then redshift will be included as a par");
    primary->card("ADDMODEL", true, "If true then this is an additive table model");
    primary->card("HDUCLASS", "OGIP", "Format conforms to OGIP standard");
    primary->card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    primary->card("HDUVERS", "1.0.0", "Version of format");

    // Save to FITS file
    fits.saveto(filename, clobber);

    // Set filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of spectra in table model
 *
 * @return Number of spectra.
 *
 * If no table model is empty the method returns zero.
 ***************************************************************************/
int GModelSpectralTable::nspectra(void) const
{
    // Initialise number of spectra
    int nspectra = 0;

    // If there are spectra then compute their number
    if (m_table_pars.size() > 0) {

        // Initialise number of spectra
        nspectra = 1;

        // Multiply number of values for all parameters
        for (int i = 0; i < m_table_pars.size(); ++i) {
            nspectra *= m_table_pars[i]->values().size();
        }

    } // endif: there were spectra

    // Return number of spectra
    return nspectra;
}


/***********************************************************************//**
 * @brief Return reference to table parameter
 *
 * @param[in] index Table parameter index [0,...,size()[.
 * @return Reference to table parameter.
 ***************************************************************************/
GModelSpectralTablePar& GModelSpectralTable::table_par(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_TABLE_PAR, "Table parameter index",
                                       index, size());
    }

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Return const reference to table parameter
 *
 * @param[in] index Table parameter index [0,...,size()[.
 * @return Const reference to table parameter.
 ***************************************************************************/
const GModelSpectralTablePar& GModelSpectralTable::table_par(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_TABLE_PAR, "Table parameter index",
                                       index, size());
    }

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Return reference to table parameter
 *
 * @param[in] name Table parameter name.
 * @return Reference to table parameter.
 ***************************************************************************/
GModelSpectralTablePar& GModelSpectralTable::table_par(const std::string& name)
{
    // Get index from name
    int index = par_index(name);

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Return const reference to table parameter
 *
 * @param[in] name Table parameter name.
 * @return Const reference to table parameter.
 ***************************************************************************/
const GModelSpectralTablePar& GModelSpectralTable::table_par(const std::string& name) const
{
    // Get index from name
    int index = par_index(name);

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Set energy scale
 *
 * @param[in] name Table parameter name for energy scaling.
 *
 * @exception GException::invalid_argument
 *            Table model parameter name not found.s
 *
 * Sets a given table parameter @p name as energy scale for table model
 * interpolation.
 ***************************************************************************/
void GModelSpectralTable::energy_scale(const std::string& name)
{
    // Continue only if energy scale is not yet set
    if ((m_escale_par.empty()) && (m_log10escale == 0)) {

        // Find parameter name in list of table parameters
        for (int i = 0; i < m_table_pars.size(); ++i) {
            GModelPar& par = m_table_pars[i]->par();
            if (par.name() == name) {
                m_escale_par = name;
                break;
            }
        }

        // Throw an exception if parameter name is not found
        if (m_escale_par.empty()) {
            std::string msg = "\"energyscale\" attribute specified but "
                              "parameter \""+name+"\" not found in the "
                              "table model. Please specify one of the "
                              "table model parameters as energy scale.";
            throw GException::invalid_argument(G_READ, msg);
        }

        // Optionally log parameter name
        #if defined(G_DEBUG_ESCALE)
        std::cout << "GModelSpectralTable::energy_scale(): ";
        std::cout << "energyscale=" << m_escale_par << std::endl;
        #endif

        // Scale energy
        scale_energy();

    } // endif: energy scale is not yet set

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print table model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing table model information.
 ***************************************************************************/
std::string GModelSpectralTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralTable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Table file"));
        result.append(m_filename.url());
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(this->size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

        // Append parameter value intervals
        for (int i = 0; i < m_table_pars.size(); ++i) {
            std::string name = m_table_pars[i]->par().name();
            int         size = m_table_pars[i]->values().size();
            result.append("\n"+gammalib::parformat(name+" values"));
            result.append(gammalib::str(size));
            if (size > 0) {
                double min = m_table_pars[i]->values()[0];
                double max = m_table_pars[i]->values()[0];
                for (int k = 0; k < size; ++k) {
                    double value = m_table_pars[i]->values()[k];
                    if (value < min) {
                        min = value;
                    }
                    if (value > max) {
                        max = value;
                    }
                }
                result.append(" ["+gammalib::str(min)+", "+gammalib::str(max)+"]");
            }
            if (m_table_pars[i]->method() == 0) {
                result.append(" (linear)");
            }
            else if (m_table_pars[i]->method() == 1) {
                result.append(" (logarithmic)");
            }
        }

        // Append energy boundaries
        result.append("\n"+gammalib::parformat("Energies"));
        result.append(gammalib::str(m_ebounds.size()));
        result.append(" [");
        result.append(m_ebounds.emin().print());
        result.append(", ");
        result.append(m_ebounds.emax().print());
        result.append("]");

        // Append shape of spectra array
        int nspectra = 0;
        int nebins   = 0;
        if (m_spectra.dim() > 0) {
            nspectra = 1;
            for (int i = 0; i < m_spectra.dim()-1; ++i) {
                nspectra *= m_spectra.shape()[i];
            }
            nebins = m_spectra.shape()[m_spectra.dim()-1];
        }
        result.append("\n"+gammalib::parformat("Spectra array dimension"));
        result.append(gammalib::str(m_spectra.dim()));
        result.append("\n"+gammalib::parformat("Number of spectra"));
        result.append(gammalib::str(nspectra));
        result.append("\n"+gammalib::parformat("Number of spectral bins"));
        result.append(gammalib::str(nebins));

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
void GModelSpectralTable::init_members(void)
{
    // Initialise normalisation
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.scale(1.0);
    m_norm.value(1.0);
    m_norm.range(0.0,1000.0);
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

    // Initialize other members
    m_table_pars.clear();
    m_spectra.clear();
    m_ebounds.clear();
    m_filename.clear();
    m_escale_par.clear();
    m_escale      = 0.0;
    m_log10escale = 0.0;

    // Initialize cache
    m_npars  = 0;
    m_nebins = 0;
    m_last_values.clear();
    m_lin_nodes.clear();
    m_log_nodes.clear();
    m_lin_values.clear();
    m_log_values.clear();

    // Initialise flux cache
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

    // Set parameter pointer(s)
    set_par_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Table model.
 ***************************************************************************/
void GModelSpectralTable::copy_members(const GModelSpectralTable& model)
{
    // Copy members
    m_norm        = model.m_norm;
    m_table_pars  = model.m_table_pars;
    m_spectra     = model.m_spectra;
    m_ebounds     = model.m_ebounds;
    m_filename    = model.m_filename;
    m_escale_par  = model.m_escale_par;
    m_escale      = model.m_escale;
    m_log10escale = model.m_log10escale;

    // Copy cache
    m_npars       = model.m_npars;
    m_nebins      = model.m_nebins;
    m_last_values = model.m_last_values;
    m_lin_nodes   = model.m_lin_nodes;
    m_log_nodes   = model.m_log_nodes;
    m_lin_values  = model.m_lin_values;
    m_log_values  = model.m_log_values;

    // Copy flux cache
    m_prefactor   = model.m_prefactor;
    m_gamma       = model.m_gamma;
    m_epivot      = model.m_epivot;
    m_flux        = model.m_flux;
    m_eflux       = model.m_eflux;

    // Copy MC cache
    m_mc_emin     = model.m_mc_emin;
    m_mc_emax     = model.m_mc_emax;
    m_mc_cum      = model.m_mc_cum;
    m_mc_min      = model.m_mc_min;
    m_mc_max      = model.m_mc_max;
    m_mc_exp      = model.m_mc_exp;

    // Set energy nodes
    set_energy_nodes();

    // Set parameter pointer(s)
    set_par_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter pointers
 ***************************************************************************/
void GModelSpectralTable::set_par_pointers(void)
{
    // Debug option: write header
    #if defined(G_DEBUG_CALLTREE)
    std::cout << "GModelSpectralTable::set_par_pointers()" << std::endl;
    #endif

    // Clear pointers
    m_pars.clear();

    // Put normalisation parameter on stack
    m_pars.push_back(&m_norm);

    // Put table model parameters on stack
    for (int i = 0; i < m_table_pars.size(); ++i) {
        m_pars.push_back(&(m_table_pars[i]->par()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy nodes from energy boundaries
 ***************************************************************************/
void GModelSpectralTable::set_energy_nodes(void)
{
    // Debug option: write header
    #if defined(G_DEBUG_CALLTREE)
    std::cout << "GModelSpectralTable::set_energy_nodes()";
    std::cout << " " << m_ebounds.size() << std::endl;
    #endif

    // Determine number of energy bins
    int nebins = m_ebounds.size();

    // Continue only if there are energy bins
    if (nebins > 0) {

        // Initialise vectors for values
        m_lin_nodes  = GNodeArray();
        m_log_nodes  = GNodeArray();
        m_lin_values = GNdarray(nebins, 1);
        m_log_values = GNdarray(nebins, 1);

        // Compute node values
        for (int i = 0; i < nebins; ++i) {
            double energy_MeV       = m_ebounds.elogmean(i).MeV();
            double log10_energy_MeV = std::log10(energy_MeV);
            m_lin_nodes.append(energy_MeV);
            m_log_nodes.append(log10_energy_MeV);
        }

        // If energy scaling is enabled then divide energy nodes by energy
        // scale
        if (has_energy_scale()) {
            for (int i = 0; i < nebins; ++i) {
                m_lin_nodes[i] /= m_escale;
                m_log_nodes[i] -= m_log10escale;
            }
        }

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Scale energy
 *
 * Scales the energy axis so that all spectra are relative to the energy
 * scaling parameter. This scaling is useful for better interpolation of
 * spectra that have a feature that depends of an energy parameter.
 ***************************************************************************/
void GModelSpectralTable::scale_energy(void)
{
    // Debug options: write header
    #if defined(G_DEBUG_CALLTREE)
    std::cout << "GModelSpectralTable::scale_energy()";
    std::cout << " \"" << m_escale_par << "\"";
    std::cout << " " << m_log10escale << std::endl;
    #elif defined(G_DEBUG_ESCALE)
    std::cout << "GModelSpectralTable::scale_energy(): entry" << std::endl;
    #endif

    // Continue only if energy scale parameter is not empty and if table
    // model was not yet scaled
    if (!(m_escale_par.empty()) && (m_escale == 0.0)) {

        // Get energy scale values in table model
        const GNodeArray& escales = m_table_pars[m_escale_par]->values();

        // Continue only if there are energy scales
        if (escales.size() > 0) {

            // Get table model dimensions, number of spectra and energy
            // scale parameter index
            int npars    = m_table_pars.size();
            int nebins   = m_ebounds.size();
            int nspectra = this->nspectra();
            int ipar     = par_index(m_escale_par);

            // Get largest energy scale in table model and its base 10
            // logarithm
            m_escale      = escales[escales.size()-1];
            m_log10escale = std::log10(m_escale);

            // Debug option: write largest energy scale
            #if defined(G_DEBUG_ESCALE)
            std::cout << " Largest energy scale: " << escale0;
            std::cout << " MeV" << std::endl;
            #endif

            // Adjust spectra according to energy scale of each spectrum
            for (int ispectrum = 0; ispectrum < nspectra; ++ispectrum) {

                // Get parameter index vector
                std::vector<int> ispec = m_spectra.index(ispectrum);

                // Get energy scale for spectrum
                double escale      = table_par(m_escale_par).values()[ispec[ipar]];
                double log10escale = std::log10(escale);

                // Debug option: write information for spectrum
                #if defined(G_DEBUG_ESCALE)
                std::cout << " " << ispectrum;
                std::cout << " " << ispec[ipar];
                std::cout << " " << escale;
                std::cout << " MeV" << std::endl;
                #endif

                // Initialise scaled spectrum
                std::vector<double> spectrum(nebins,0.0);

                // Initialise vector for spectrum access
                std::vector<int> index(npars+1,0);
                for (int k = 0; k < npars; ++k) {
                    index[k] = ispec[k];
                }

                // Setup scale spectrum
                for (int iebin = 0; iebin < nebins; ++iebin) {

                    // Get scaled log10 energy for energy bin
                    double log10energy = m_log_nodes[iebin] - m_log10escale + log10escale;

                    // Interpolate function in log10 energy
                    m_log_nodes.set_value(log10energy);
                    double wgt_left  = m_log_nodes.wgt_left();
                    double wgt_right = m_log_nodes.wgt_right();
                    index[npars]     = m_log_nodes.inx_left();
                    double val_left  = m_spectra(index);
                    index[npars]     = m_log_nodes.inx_right();
                    double val_right = m_spectra(index);
                    spectrum[iebin]  = wgt_left * val_left + wgt_right * val_right;

                } // endfor: scaled spectrum

                // Save scaled spectrum
                index[npars] = 0;
                for (int iebin = 0; iebin < nebins; ++iebin, ++index[npars]) {
                    m_spectra(index) = spectrum[iebin];
                }

            } // endfor: looped over all spectra

            // Divide energy nodes by largest energy scale value
            for (int iebin = 0; iebin < nebins; ++iebin) {
                m_lin_nodes[iebin] /= m_escale;
                m_log_nodes[iebin] -= m_log10escale;
            }

        } // endif: there were energy scales

        // Force update
        m_last_values.clear();

    } // endif: energy scale parameter was not empty

    // Debug option: write trailer
    #if defined(G_DEBUG_ESCALE)
    std::cout << "GModelSpectralTable::update(): exit" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Signal that energy scale was set
 *
 * @return True if energy scale was set
 ***************************************************************************/
bool GModelSpectralTable::has_energy_scale(void) const
{
    // Set energy scale flag
    bool flag = (m_escale != 0.0);

    // Return flag
    return flag;
}


/***********************************************************************//**
 * @brief Create PARAMETERS FITS table
 *
 * @return Binary FITS table containing table model parameters
 *
 * The method creates a binary FITS table that is compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node6.html
 ***************************************************************************/
GFitsBinTable GModelSpectralTable::create_par_table(void) const
{
    // Determine number of parameters
    int nrows = m_table_pars.size();

    // Determine maximum number of parameters values
    int ncols = 0;
    for (int i = 0; i < nrows; ++i) {
        int n = m_table_pars[i]->size();
        if (n > ncols) {
            ncols = n;
        }
    }

    // Create table columns
    GFitsTableStringCol col_name("NAME", nrows, 12);
    GFitsTableLongCol   col_method("METHOD", nrows);
    GFitsTableFloatCol  col_initial("INITIAL", nrows);
    GFitsTableFloatCol  col_delta("DELTA", nrows);
    GFitsTableFloatCol  col_minimum("MINIMUM", nrows);
    GFitsTableFloatCol  col_bottom("BOTTOM", nrows);
    GFitsTableFloatCol  col_top("TOP", nrows);
    GFitsTableFloatCol  col_maximum("MAXIMUM", nrows);
    GFitsTableLongCol   col_numbvals("NUMBVALS", nrows);
    GFitsTableFloatCol  col_value("VALUE", nrows, ncols);

    // Fill columns
    for (int i = 0; i < nrows; ++i) {
    
        // Get model parameter
        const GModelPar& par = m_table_pars[i]->par();

        // Set parameter name and initial value
        col_name(i)    = par.name();
        col_method(i)  = m_table_pars[i]->method();
        col_initial(i) = (float)par.value();

        // Handle free/fixed parameter attribute
        if (par.is_fixed()) {
            col_delta(i) = -1.0;
        }
        else {
            col_delta(i) =  1.0;  // Dummy value
        }

        // Set number of parameters
        col_numbvals(i) = m_table_pars[i]->size();

        // Set parameter values
        double min_value = 0.0;
        double max_value = 0.0;
        if (m_table_pars[i]->size() > 0) {
            min_value = m_table_pars[i]->values()[0];
            max_value = m_table_pars[i]->values()[0];
            for (int k = 0; k < m_table_pars[i]->size(); ++k) {
                double value = m_table_pars[i]->values()[k];
                if (value < min_value) {
                    min_value = value;
                }
                if (value > max_value) {
                    max_value = value;
                }
                col_value(i,k) = value;
            }
        }

        // Handle parameter minimum
        if (par.has_min()) {
            if (par.min() < min_value) {
                col_minimum(i) = (float)par.min();  // Hard minimum limit
                col_bottom(i)  = min_value;         // Soft minimum limit
            }
            else {
                col_minimum(i) = (float)par.min();  // Hard minimum limit
                col_bottom(i)  = (float)par.min();  // Soft minimum limit
            }
        }
        else {
            col_minimum(i) = min_value;             // Hard minimum limit
            col_bottom(i)  = min_value;             // Soft minimum limit
        }

        // Handle parameter maximum
        if (par.has_max()) {
            if (par.max() > max_value) {
                col_top(i)     = max_value;         // Soft maximum limit
                col_maximum(i) = (float)par.max();  // Hard maximum limit
            }
            else {
                col_top(i)     = (float)par.max();  // Soft maximum limit
                col_maximum(i) = (float)par.max();  // Hard maximum limit
            }
        }
        else {
            col_top(i)     = max_value;             // Soft maximum limit
            col_maximum(i) = max_value;             // Soft maximum limit
        }

    } // endfor: looped over parameters

    // Create binary table
    GFitsBinTable table;

    // Append columns to FITS table
    table.append(col_name);
    table.append(col_method);
    table.append(col_initial);
    table.append(col_delta);
    table.append(col_minimum);
    table.append(col_bottom);
    table.append(col_top);
    table.append(col_maximum);
    table.append(col_numbvals);
    table.append(col_value);

    // Set table extension name
    table.extname("PARAMETERS");

    // Set table keywords
    table.card("HDUCLASS", "OGIP", "Format conforms to OGIP standard");
    table.card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    table.card("HDUCLAS2", "PARAMETERS", "Extension containing parameter info");
    table.card("HDUVERS", "1.0.0", "Version of format");
    table.card("NINTPARM", nrows, "Number of interpolation parameters");
    table.card("NADDPARM", 0, "Number of additional parameters");

    // Return table
    return table;
}


/***********************************************************************//**
 * @brief Create ENERGIES FITS table
 *
 * @return Binary FITS table containing table model energies
 *
 * The method creates a binary FITS table that is compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node7.html
 ***************************************************************************/
GFitsBinTable GModelSpectralTable::create_eng_table(void) const
{
    // Create table columns
    GFitsTableFloatCol col_lo("ENERG_LO", m_ebounds.size());
    GFitsTableFloatCol col_hi("ENERG_HI", m_ebounds.size());

    // Fill energy boundary columns
    for (int i = 0; i < m_ebounds.size(); ++i) {
        col_lo(i) = m_ebounds.emin(i).keV();
        col_hi(i) = m_ebounds.emax(i).keV();
    }

    // Set energy units
    col_lo.unit("keV");
    col_hi.unit("keV");

    // Create binary table
    GFitsBinTable table;

    // Append columns to FITS table
    table.append(col_lo);
    table.append(col_hi);

    // Set table extension name
    table.extname("ENERGIES");

    // Set table keywords
    table.card("HDUCLASS", "OGIP", "Format conforms to OGIP standard");
    table.card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    table.card("HDUCLAS2", "ENERGIES", "Extension containing energy bin info");
    table.card("HDUVERS", "1.0.0", "Version of format");

    // Return table
    return table;
}


/***********************************************************************//**
 * @brief Create SPECTRA FITS table
 *
 * @return Binary FITS table containing table model spectra
 *
 * The method creates a binary FITS table that is compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node8.html
 ***************************************************************************/
GFitsBinTable GModelSpectralTable::create_spec_table(void) const
{
    // Compute number of rows
    int nrows = 1;
    for (int i = 0; i < m_spectra.dim()-1; ++i) {
        nrows *= m_spectra.shape()[i];
    }

    // Compute number of parameters
    int npars = m_table_pars.size();

    // Compute number of energy bins
    int nebins = m_ebounds.size();

    // Create columns
    GFitsTableFloatCol col_pars("PARAMVAL", nrows, npars);
    GFitsTableFloatCol col_spec("INTPSPEC", nrows, nebins);

    // Set units
    col_spec.unit("ph cm-2 s-1 MeV-1");

    // Fill columns
    std::vector<int> inx(npars+1,0);
    for (int i = 0; i < nrows; ++i) {

        // Set parameter values
        for (int k = 0; k < npars; ++k) {
            col_pars(i,k) = m_table_pars[k]->values()[inx[k]];
        }

        // Set current index
        std::vector<int> index = inx;

        // Loop over energy bins
        for (int k = 0; k < nebins; ++k, ++index[npars]) {
            col_spec(i,k) = m_spectra(index);
        }

        // Increment parameter index. Last parameter index is changing fastest
        int ipar = npars-1;
        do {
            inx[ipar] += 1;
            if (inx[ipar] < m_spectra.shape()[ipar]) {
                break;
            }
            else {
                inx[ipar] = 0;
                ipar--;
            }
        } while (ipar >= 0);

    } // endfor: looped over rows

    // Create binary table
    GFitsBinTable table;

    // Append columns to FITS table
    table.append(col_pars);
    table.append(col_spec);

    // Set table extension name
    table.extname("SPECTRA");

    // Set table keywords
    table.card("HDUCLASS", "OGIP",              "Format conforms to OGIP standard");
    table.card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    table.card("HDUCLAS2", "MODEL SPECTRA",     "Extension containing model spectra");
    table.card("HDUVERS",  "1.0.0",             "Version of format");

    // Return table
    return table;
}


/***********************************************************************//**
 * @brief Load data from PARAMETERS extension
 *
 * @param[in] fits FITS file.
 *
 * @exception GException::invalid_value
 *            Non-positive parameter value encountered for logarithmic
 *            parameters.
 *
 * The method loads data from the PARAMETERS binary table. The format of the
 * table needs to be compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node6.html
 ***************************************************************************/
void GModelSpectralTable::load_par(const GFits& fits)
{
    // Get PARAMETERS extension
    const GFitsTable& table = *fits.table("PARAMETERS");

    // Get number of parameters
    int npars = table.integer("NAXIS2");

    // Loop over all parameters
    for (int i = 0; i < npars; ++i) {

        // Set model parameter
        GModelPar par(table["NAME"]->string(i), table["INITIAL"]->real(i));

        // Apply hard minimum and maximum as parameter range. Note that the
        // minimum and maximum apply to the value factor, hence in case of
        // a negative scale factor the minimum becomes the maximum and vice
        // versa. This is actually a bug in GammaLib, see
        // https://cta-redmine.irap.omp.eu/issues/3072
        double min = table["MINIMUM"]->real(i) / par.scale();
        double max = table["MAXIMUM"]->real(i) / par.scale();
        if (min > max) {
            par.factor_range(max, min);
        }
        else {
            par.factor_range(min, max);
        }

        // Fix or free parameter according to DELTA value
        if (table["DELTA"]->real(i) < 0.0) {
            par.fix();
            par.has_grad(false);
        }
        else {
            par.free();
            par.has_grad(true);
        }

        // Set vector of parameter values
        std::vector<double> values;
        for (int k = 0; k < table["NUMBVALS"]->integer(i); ++k) {

            // Extract value
            double value = table["VALUE"]->real(i,k);

            // If method is logarithmic then store log of the parameter
            // values in the spectral table parameters
            /*
            if (table["METHOD"]->integer(i) == 1) {
                if (value > 0.0) {
                    value = std::log(value);
                }
                else {
                    std::string msg = "Non-positive value "+gammalib::str(value)+
                                      " encountered for logarithmic parameter "
                                      "\""+table["NAME"]->string(i)+"\". Please "
                                      "provide only positive parameter "
                                      "values for logarithmic parameters.";
                    throw GException::invalid_value(G_LOAD_PAR, msg);
                }
            }
            */

            // Append value
            values.push_back(value);

        } // endfor: looped over parameter values

        // Set table model parameter
        GModelSpectralTablePar table_model_par(par, values);

        // Set interpolation method
        table_model_par.method(table["METHOD"]->integer(i));

        // Append table model parameter
        m_table_pars.append(table_model_par);

    } // endfor: looped over all parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data from ENERGIES extension
 *
 * @param[in] fits FITS file.
 *
 * The method loads data from the ENERGIES binary table. The format of the
 * table needs to be compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node7.html
 *
 * If energy units are provided in the ENERGIES extension the method decodes
 * the units and interprets the energy values correspondingly.
 ***************************************************************************/
void GModelSpectralTable::load_eng(const GFits& fits)
{
    // Get ENERGIES extension
    const GFitsTable& table = *fits.table("ENERGIES");

    // Get number of energy bins
    int nebins = table.integer("NAXIS2");

    // Extract energy boundary information from FITS table
    if (nebins > 0) {

        // Get units
        std::string emin_unit = table["ENERG_LO"]->unit();
        std::string emax_unit = table["ENERG_HI"]->unit();
        if (emin_unit.empty()) {
            emin_unit = "keV";
        }
        if (emax_unit.empty()) {
            emax_unit = "keV";
        }

        // Read energy boundaries and append bins
        for (int i = 0; i < nebins; ++i) {
            GEnergy emin(table["ENERG_LO"]->real(i), emin_unit);
            GEnergy emax(table["ENERG_HI"]->real(i), emax_unit);
            m_ebounds.append(emin, emax);
        }

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data from SPECTRA extension
 *
 * @param[in] fits FITS file.
 *
 * The method loads data from the SPECTRA binary table. The format of the
 * table needs to be compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node8.html
 ***************************************************************************/
void GModelSpectralTable::load_spec(const GFits& fits)
{
    // Get number of parameters and energy bins
    int npars  = fits.table("PARAMETERS")->integer("NAXIS2");
    int nebins = fits.table("ENERGIES")->integer("NAXIS2");

    // Setup dimension of spectra array
    std::vector<int> naxis(npars+1, 0);
    for (int i = 0; i < npars; ++i) {
        naxis[i] = (*fits.table("PARAMETERS"))["NUMBVALS"]->integer(i);
    }
    naxis[npars] = nebins;

    // Setup spectra array
    m_spectra = GNdarray(naxis);

    // Get SPECTRA extension
    const GFitsTable& table = *fits.table("SPECTRA");

    // Get number of rows
    int nrows = table.integer("NAXIS2");

    // Loop over rows
    for (int i = 0; i < nrows; ++i) {

        // Initialise spectra array index of first energy bin
        std::vector<int> index(npars+1, 0);

        // Setup spectra array index
        std::vector<double> parval(npars, 0.0);
        for (int k = 0; k < npars; ++k) {

            // Get reference to node array
            const GNodeArray& nodes = m_table_pars[k]->values();

            // Set interpolation value
            nodes.set_value(table["PARAMVAL"]->real(i,k));

            // Get best index
            int inx = (nodes.wgt_left() > nodes.wgt_right()) ? nodes.inx_left()
                                                             : nodes.inx_right();

            // Store index
            index[k] = inx;

        } // endfor: looped over parameters

        // Debug: dump vector array
        #if defined(G_DEBUG_LOAD_SPEC)
        std::cout << i << ": (";
        for (int k = 0; k < index.size(); ++k) {
            if (k > 0) {
                std::cout << ", ";
            }
            std::cout << index[k];
        }
        std::cout << ")" << std::endl;
        #endif

        // Loop over energy bins and store spectrum
        for (int k = 0; k < nebins; ++k, ++index[npars]) {
            m_spectra(index) = table["INTPSPEC"]->real(i,k);
        }

    } // endfor: looped over rows

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return index for parameter name
 *
 * @param[in] name Parameter name.
 * @return Parameter index.
 *
 * @exception GException::invalid_argument
 *            Parameter name not found in spectral table.
 ***************************************************************************/
int GModelSpectralTable::par_index(const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_table_pars[index]->par().name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Parameter name \""+name+"\" not found in spectral "
                          "table. Please specify one of the following parameter "
                          "names:";
        for (int i = 0; i < size(); ++i) {
            if (i > 0) {
                msg += ",";
            }
            msg += " \""+m_table_pars[i]->par().name()+"\"";
        }
        throw GException::invalid_argument(G_PAR_INDEX, msg);
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Update cache for spectral table model computation
 *
 * Update interval vectors for function values and gradients. An update is
 * performed in case that some of the parameter values have changed. The
 * method sets the following cache members:
 *
 *      m_npars (number of model parameters)
 *      m_nebins (number of energy bins in spectra)
 *      m_last_values (last model parameter values)
 *      m_lin_values (function values)
 *
 * The array \f${\tt m\_lin\_values}\f$ holds the vector \f$F(E|p)\f$ as
 * function of energy \f$E\f$, computed for the current set of parameters
 * \f$p\f$. The computation is done using the N-dimensional linear
 * interpolation
 *
 * \f[
 *    F(E|p) = \sum_{k=1}^M \left( \prod_{i=1}^N w_x(p) \right) F_p(E)
 * \f]
 *
 * where
 * - \f$M=2^N\f$ is the number of parameter combinations,
 * - \f$N\f$ is the number of table model parameters,
 * - \f$w_x(p) = \{w_0^l, w_0^r, w_1^l, w_1^r, ...\}\f$
 *   is an array of subsequently arranged left and right weighting factors
 *   for the linear interpolation of parameters, where \f$w_0^l\f$ and
 *   \f$w_0^r\f$ are the left and right weighting factors for the first
 *   parameter, \f$w_1^l\f$ and \f$w_1^r\f$ are the left and right weighting
 *   factors for the second parameter, and so on,
 * - \f$x=2 k + i/2^k \mod 2\f$ is the index in the array of
 *   weighting factors for parameter combination \f$k\f$ and parameter
 *   \f$i\f$, and
 * - \f$F_p(E)\f$ are the table model spectra, computed for a grid of
 *   possible parameter values.
 *
 * For each parameter \f$i\f$, the weighting factors \f$w_i^l(p)\f$
 * and \f$w_i^r(p)\f$ are computed using
 *
 * \f[
 *    w_i^r(p) = \frac{p - p_l}{p_r - p_l}
 * \f]
 *
 * and \f$w_i^l(p) = 1 - w_i^r(p)\f$, where
 * \f$p_l\f$ is the largest parameter value that satisfies \f$p<p_l\f$ and
 * \f$p_r\f$ is the smallest parameter value that satisfies \f$p>p_r\f$.
 *
 * The method also computes the gradients
 *
 * \f[
 *    \frac{\delta F(E|p)}{\delta p} = \sum_{k=1}^M \left(
 *    \frac{\delta w_{x_p}(p)}{\delta p} \prod_{i=1 \land i \neq i_p}^N
 *    w_x(p) \right) F_p(E)
 * \f]
 *
 * where
 * \f$x_p = 2 k + i_p/2^k \mod 2\f$.
 ***************************************************************************/
void GModelSpectralTable::update(void) const
{
    // Get number of parameters and number of energy bins
    m_npars  = m_table_pars.size();
    m_nebins = ebounds().size();

    // Initialise update flag
    bool need_update = false;

    // If dimension of last cached parameter values differ from dimension
    // of spectral table then reallocate cache and request update
    if (m_last_values.size() != m_npars) {
        m_last_values = std::vector<double>(m_npars, 0.0);
        need_update   = true;
    }

    // ... otherwise check if some parameter values have changed
    else {
        for (int i = 0; i < m_npars; ++i) {
            if (m_table_pars[i]->par().value() != m_last_values[i]) {
                need_update = true;
                break;
            }
        }
    }

    // Continue only if update is required
    if (need_update) {

        // Debug option: write header
        #if defined(G_DEBUG_UPDATE)
        std::cout << "GModelSpectralTable::update() required" << std::endl;
        #endif

        // Initialise vectors for weights, weight gradients and indices
        std::vector<double> weights(2*m_npars, 0.0);
        std::vector<double> weight_gradients(2*m_npars, 0.0);
        std::vector<int>    indices(2*m_npars, 0);

        // Loop over all parameters and extract left and right weights,
        // weight gradients and indices and put them into single arrays
        for (int i = 0; i < m_npars; ++i) {

            // Get pointers to node array and parameter
            const GNodeArray* nodes = &(m_table_pars[i]->values());
            const GModelPar*  par   = &(m_table_pars[i]->par());

            // Get parameter value
            double value = par->value();

            // Cache parameter value
            m_last_values[i] = value;

            // Use log of value for logarithmic parameters
            /*
            if (m_table_pars[i]->method() == 1) {
                if (value > 0.0) {
                    value = std::log(value);
                }
                else {
                    std::string msg = "Non-positive value "+gammalib::str(value)+
                                      " encountered for logarithmic parameter "
                                      "\""+m_table_pars[i]->par().name()+"\".";
                    throw GException::invalid_value(G_UPDATE, msg);
                }
            }
            */

            // Set values for node array
            nodes->set_value(value);

            // Compute left and right indices
            int il = 2*i;
            int ir = il + 1;

            // Push back weigths, weight gradients and indices
            weights[il]          = nodes->wgt_left();
            weights[ir]          = nodes->wgt_right();
            weight_gradients[il] = nodes->wgt_grad_left();
            weight_gradients[ir] = nodes->wgt_grad_right();
            indices[il]          = nodes->inx_left();
            indices[ir]          = nodes->inx_right();

            // Debug option: print weights and indices
            #if defined(G_DEBUG_UPDATE)
            std::cout << " wgt_l=" << weights[il];
            std::cout << " wgt_r=" << weights[ir];
            std::cout << " wgt_grad_l=" << weight_gradients[il];
            std::cout << " wgt_grad_r=" << weight_gradients[ir];
            std::cout << " inx_l=" << indices[il];
            std::cout << " inx_r=" << indices[ir] << std::endl;
            #endif

        } // endfor: looped over all parameters

        // Compute number of combinations
        int combinations = 1 << m_npars;

        // Initialise 2d arrays for values and gradients
        m_lin_values = GNdarray(m_nebins, m_npars+1);
        m_log_values = GNdarray(m_nebins, m_npars+1);

        // Debug option: initial sum of weights
        #if defined(G_DEBUG_UPDATE)
        double weight_sum = 0.0;
        #endif

        // Loop over combinations
        for (int i = 0; i < combinations; ++i) {

            // Debug option: start printing combination
            #if defined(G_DEBUG_UPDATE)
            std::cout << " " << i << ": ";
            #endif

            // Initialise weight and gradient weights
            double              weight = 1.0;
            std::vector<double> grad_weight(m_npars, 1.0);

            // Initialise index vector (including the energy dimension)
            std::vector<int> index_shape(m_npars+1,0);

            // Loop over dimensions
            for (int k = 0, div = 1; k < m_npars; ++k, div *= 2) {

                // Compute index for each dimension
                int index = i/div % 2 + k * 2;

                // Update weight
                weight *= weights[index];

                // Add index
                index_shape[k] = indices[index];

                // Update gradient weights
                for (int j = 0; j < m_npars; ++j) {
                    if (k == j) {
                        grad_weight[j] *= weight_gradients[index];
                    }
                    else {
                        grad_weight[j] *= weights[index];
                    }
                } // endfor: update gradient weights

                // Debug option: print index and weight for current dimension
                #if defined(G_DEBUG_UPDATE)
                std::cout << index;
                std::cout << " (" << weights[index];
                std::cout << " @ " << indices[index] << ") ";
                #endif

            } // endfor: looped over dimensions

            // Debug option: print total weight and weight gradient
            #if defined(G_DEBUG_UPDATE)
            std::cout << ": wgt=" << weight;
            std::cout << " [";
            for (int k = 0; k < m_npars; ++k) {
                if (k > 0) {
                    std::cout << ",";
                }
                std::cout << grad_weight[k];
            }
            std::cout << "]" << std::endl;
            weight_sum += weight;
            #endif

            // Compute interpolated value and gradient
            for (int iebin = 0; iebin < m_nebins; ++iebin) {

                // Set energy index
                index_shape[m_npars] = iebin;

                // Get spectral value
                double value = m_spectra(index_shape);

                // Compute contribution and store in value slot
                m_lin_values(iebin,0) += weight * value;

                // Compute gradients and store in gradient slots
                for (int j = 0; j < m_npars; ++j) {
                    m_lin_values(iebin,j+1) += grad_weight[j] * value;
                }

            } // endfor: looped over all energy bins

        } // endfor: looped over combinations

        // Compute log10 values of function values
        /*
        for (int iebin = 0; iebin < m_nebins; ++iebin) {
            if (m_lin_values(iebin,0) > 0.0) {
                m_log_values(iebin,0) = std::log10(m_lin_values(iebin,0));
            }
            else {
                m_log_values(iebin,0) = -1000.0; // Set to a tiny value
            }
        }
        */

        // Debug option: print sum of weights
        #if defined(G_DEBUG_UPDATE)
        std::cout << " sum(wgt)=" << weight_sum << std::endl;
        #endif

    } // endif: updated requested

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update flux cache
 ***************************************************************************/
void GModelSpectralTable::update_flux(void) const
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
        double fmin = m_lin_values(i,0);
        double fmax = m_lin_values(i+1,0);

        // Compute pivot energy (MeV). We use here the geometric mean of
        // the node boundaries.
        double epivot = std::sqrt(emin*emax);

        // Compute spectral index, power law normalisation, photon and
        // energy flux between nodes
        double gamma     = 0.0;
        double prefactor = 0.0;
        double flux      = 0.0;
        double eflux     = 0.0;
        if ((fmin > 0.0) && (fmax > 0.0)) {
            gamma     = std::log(fmin/fmax) / std::log(emin/emax);
            prefactor = fmin / std::pow(emin/epivot, gamma);
            flux      = prefactor *
                        gammalib::plaw_photon_flux(emin, emax, epivot, gamma);
            eflux     = prefactor *
                        gammalib::plaw_energy_flux(emin, emax, epivot, gamma);
        }

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
 * @brief Update MC cache
 *
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * This method sets up an array of indices and the cumulative distribution
 * function needed for MC simulations.
 ***************************************************************************/
void GModelSpectralTable::update_mc(const GEnergy& emin,
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
                       gammalib::plaw_photon_flux(e_min,
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
                       gammalib::plaw_photon_flux(e_min,
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
                       gammalib::plaw_photon_flux(m_lin_nodes[inx_emax],
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
            if (norm > 0.0) {
                for (int i = 0; i < m_mc_cum.size(); ++i) {
                    m_mc_cum[i] /= norm;
                }
            }

            // Set MC values
            for (int i = 0; i < m_mc_cum.size(); ++i) {

                // Compute exponent
                double exponent = m_mc_exp[i] + 1.0;

                // Exponent dependend computation
                if (std::abs(exponent) > 1.0e-11) {

                    // If the exponent is too small then use lower energy
                    // boundary
                    if (exponent < -50.0) {
                        m_mc_exp[i] = 0.0;
                        m_mc_min[i] = std::log(m_mc_min[i]);
                        m_mc_max[i] = m_mc_min[i];
                    }

                    // ... otherwise if exponent is too large then use
                    // upper energy boundary
                    else if (exponent > +50.0) {
                        m_mc_exp[i] = 0.0;
                        m_mc_min[i] = std::log(m_mc_max[i]);
                        m_mc_max[i] = m_mc_min[i];
                    }

                    // ... otherwise use transformation formula
                    else {
                        m_mc_exp[i] = exponent;
                        m_mc_min[i] = std::pow(m_mc_min[i], exponent);
                        m_mc_max[i] = std::pow(m_mc_max[i], exponent);
                    }
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
