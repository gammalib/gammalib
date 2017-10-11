/***************************************************************************
 *          GCTAOnOffObservation.cpp - CTA On/Off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Chia-Chun Lu & Christoph Deil               *
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
 * @file GCTAOnOffObservation.cpp
 * @brief CTA On/Off observation class implementation
 * @author Chia-Chun Lu & Christoph Deil & Pierrick Martin
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GObservationRegistry.hpp"
#include "GTools.hpp"
#include "GModels.hpp"
#include "GModelSky.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GSkyRegionMap.hpp"
#include "GOptimizerPars.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTACubeBackground.hpp"
#include "GCTAModelIrfBackground.hpp"
#include "GCTAOnOffObservation.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Globals ____________________________________________________________ */
const GCTAOnOffObservation g_onoff_obs_cta_seed;
const GObservationRegistry g_onoff_obs_cta_registry(&g_onoff_obs_cta_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE_SET           "GCTAOnOffObservation::response(GResponse&)"
#define G_RESPONSE_GET                     "GCTAOnOffObservation::response()"
#define G_WRITE                   "GCTAOnOffObservation::write(GXmlElement&)"
#define G_READ                     "GCTAOnOffObservation::read(GXmlElement&)"
#define G_LIKELIHOOD            "GCTAOnOffObservation::likelihood(GModels&, "\
                                          "GOptimizerPars&, GMatrixSparse&, "\
                                                "GVector&, double&, double&)"
#define G_SET                   "GCTAOnOffObservation::set(GCTAObservation&)"
#define G_COMPUTE_ARF   "GCTAOnOffObservation::compute_arf(GCTAObservation&)"
#define G_COMPUTE_BGD   "GCTAOnOffObservation::compute_bgd(GCTAObservation&)"
#define G_COMPUTE_ALPHA                "GCTAOnOffObservation::compute_alpha("\
                                                          "GCTAObservation&)"
#define G_COMPUTE_RMF   "GCTAOnOffObservation::compute_rmf(GCTAObservation&)"

/* __ Constants __________________________________________________________ */
const double minmod = 1.0e-100;                      //!< Minimum model value
const double minerr = 1.0e-100;                //!< Minimum statistical error

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_LIKELIHOOD_DEBUG                //!< Debug likelihood computation


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty On/Off observation.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(void) : GObservation()
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs On/Off observation.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GCTAOnOffObservation& obs) :
                      GObservation(obs)
{ 
    // Initialise private
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief CTA observation constructor
 *
 * @param[in] obs CTA observation.
 * @param[in] etrue True energy boundaries.
 * @param[in] ereco Reconstructed energy boundaries.
 * @param[in] on On regions.
 * @param[in] off Off regions.
 *
 * Constructs On/Off observation a CTA observation by filling the On and Off
 * spectra and computing the Auxiliary Response File (ARF) and Redistribution
 * Matrix File (RMF). The method requires the specification of the true and
 * reconstructed energy boundaries, as well as the definition of On and Off
 * regions.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GCTAObservation& obs,
                                           const GEbounds&        etrue,
                                           const GEbounds&        ereco,
                                           const GSkyRegionMap&   on,
                                           const GSkyRegionMap&   off,
                                           const GSkyDir          src_dir)
{
    // Initialise private
    init_members();

    // Initialise spectra
    m_on_spec  = GPha(ereco);
    m_off_spec = GPha(ereco);

    // Initialise response information
    m_arf = GArf(etrue);
    m_rmf = GRmf(etrue, ereco);

    // Store regions
    m_on_regions  = on;
    m_off_regions = off;
    
    // Set direction to source
    m_src_dir=src_dir;

    // Set On/Off observation from CTA observation
    set(obs);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAOnOffObservation::~GCTAOnOffObservation(void)
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
 * @param[in] obs On/Off observation.
 * @return On/Off observation.
 *
 * Assigns one On/Off observation to another On/Off observation object.
 ***************************************************************************/
GCTAOnOffObservation& GCTAOnOffObservation::operator=(const GCTAOnOffObservation& obs)
{ 
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * Clears the On/Off observation. All class members will be set to the
 * initial state. Any information that was present in the object before will
 * be lost.
 ***************************************************************************/
void GCTAOnOffObservation::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of On/Off observation.
 *
 * Returns a pointer to a deep copy of an On/Off observation.
 **************************************************************************/
GCTAOnOffObservation* GCTAOnOffObservation::clone(void) const
{
    return new GCTAOnOffObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GException::invalid_argument
 *            Invalid response class specified.
 *
 * Sets the response function for the On/Off observation.
 ***************************************************************************/
void GCTAOnOffObservation::response(const GResponse& rsp)
{
    // Cast response dynamically
    const GCTAResponse* ptr = dynamic_cast<const GCTAResponse*>(&rsp);

    // Throw exception if response is not of correct type
    if (ptr == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Invalid response type \""+cls+"\" provided on "
                          "input. Please specify a \"GCTAResponse\" "
                          "as argument.";
        throw GException::invalid_argument(G_RESPONSE_SET, msg);
    }

    // Free response
    if (m_response != NULL) delete m_response;

    // Clone response function
    m_response = ptr->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer to CTA response function
 *
 * @return Pointer to CTA response function.
 *
 * @exception GException::invalid_value
 *            No valid response found in CTA observation.
 *
 * Returns a pointer to the CTA response function. An exception is thrown if
 * the pointer is not valid, hence the user does not need to verify the
 * validity of the pointer.
 ***************************************************************************/
const GCTAResponse* GCTAOnOffObservation::response(void) const
{
    // Throw an exception if the response pointer is not valid
    if (m_response == NULL) {
        std::string msg = "No valid response function found in CTA On/Off "
                          "observation.\n";
        throw GException::invalid_value(G_RESPONSE_GET, msg);
    }

    // Return pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Read On/Off observation from an XML element
 *
 * @param[in] xml XML element.
 *
 * Reads information for an On/Off observation from an XML element. The
 * expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="Pha_on"      file="..."/>
 *       <parameter name="Pha_off"     file="..."/>
 *       <parameter name="Regions_on"  file="..."/>
 *       <parameter name="Regions_off" file="..."/>
 *       <parameter name="Arf"         file="..."/>
 *       <parameter name="Rmf"         file="..."/>
 *     </observation>
 *
 ***************************************************************************/
void GCTAOnOffObservation::read(const GXmlElement& xml)
{
    // clean object
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Get file names
    std::string pha_on  = gammalib::xml_get_attr(G_READ, xml, "Pha_on",      "file");
    std::string pha_off = gammalib::xml_get_attr(G_READ, xml, "Pha_off",     "file");
    std::string reg_on  = gammalib::xml_get_attr(G_READ, xml, "Regions_on",  "file");
    std::string reg_off = gammalib::xml_get_attr(G_READ, xml, "Regions_off", "file");
    std::string arf     = gammalib::xml_get_attr(G_READ, xml, "Arf",         "file");
    std::string rmf     = gammalib::xml_get_attr(G_READ, xml, "Rmf",         "file");

    // Expand file names
    pha_on  = gammalib::xml_file_expand(xml, pha_on);
    pha_off = gammalib::xml_file_expand(xml, pha_off);
    reg_on  = gammalib::xml_file_expand(xml, reg_on);
    reg_off = gammalib::xml_file_expand(xml, reg_off);
    arf     = gammalib::xml_file_expand(xml, arf);
    rmf     = gammalib::xml_file_expand(xml, rmf);

    // Load files
    m_on_spec.load(pha_on);
    m_off_spec.load(pha_off);
    m_on_regions.load(reg_on);
    m_off_regions.load(reg_off);
    m_arf.load(arf);
    m_rmf.load(rmf);

	// Return
	return;
}


/***********************************************************************//**
 * @brief write observation to an xml element
 *
 * @param[in] xml XML element.
 *
 * Writes information for an On/Off observation into an XML element. The
 * expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="Pha_on"      file="..."/>
 *       <parameter name="Pha_off"     file="..."/>
 *       <parameter name="Regions_on"  file="..."/>
 *       <parameter name="Regions_off" file="..."/>
 *       <parameter name="Arf"         file="..."/>
 *       <parameter name="Rmf"         file="..."/>
 *     </observation>
 *
 * The actual files described in the XML elements are not written.
 ***************************************************************************/
void GCTAOnOffObservation::write(GXmlElement& xml) const
{
    // Allocate XML element pointer
    GXmlElement* par;

    // Set Pha_on parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Pha_on");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_on_spec.filename()));

    // Set Pha_off parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Pha_off");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_off_spec.filename()));

    // Set Regions_on parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Regions_on");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_on_regions.filename()));

    // Set Regions_off parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Regions_off");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_off_regions.filename()));

    // Set Arf parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Arf");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_arf.filename()));

    // Set Rmf parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Rmf");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_rmf.filename()));

	// Return
	return;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for On/Off analysis
 *
 * @param[in] models Models.
 * @param[in,out] gradient Pointer to gradients.
 * @param[in,out] curvature Pointer to curvature matrix.
 * @param[in,out] npred Pointer to Npred value.
 * @return Log-likelihood value.
 *
 * @exception GException::invalid_value
 *            There are no model parameters.
 *
 * Computes the log-likelihood value for the On/Off observation. The
 * method loops over the energy bins to update the function value, its
 * derivatives and the curvature matrix. The number of On counts
 * \f$N_{\rm on}\f$ and Off counts \f$N_{\rm off}\f$ are taken from the
 * On and Off spectra, the expected number of gamma-ray events
 * \f$N_{\gamma}\f$ and background events \f$N_{\rm bgd}\f$ are
 * computed from the spectral models of the relevant components in the
 * model container (spatial and temporal components are ignored so far).
 * See the N_gamma() and N_bgd() methods for details about the model
 * computations.
 *
 * The log-likelihood is given by
 *
 * \f[
 *    \ln L = \sum_i \ln L_i
 * \f]
 *
 * where the sum is taken over all energy bins \f$i\f$ and
 *
 * \f[
 *    \ln L_i = - N_{\rm on}  \ln N_{\rm pred} + N_{\rm pred}
 *              - N_{\rm off} \ln N_{\rm bgd}  + N_{\rm bgd}
 * \f]
 *
 * with
 *
 * \f[
 *    N_{\rm pred} = N_{\gamma} + \alpha N_{\rm bgd}
 * \f]
 *
 * being the total number of predicted events for an energy bin in the On
 * region,
 * \f$N_{\rm on}\f$ is the total number of observed events for an energy
 * bin in the On region,
 * \f$N_{\rm off}\f$ is the total number of observed events for an energy
 * bin in the Off region, and
 * \f$N_{\rm bgd}\f$ is the predicted number of background events for an
 * energy bin in the Off region.
 *
 * The log-likelihood gradient with respect to sky model parameters
 * \f$p_{\rm sky}\f$ is given by
 *
 * \f[
 *    \left( \frac{\partial \ln L_i}{\partial p_{\rm sky}} \right) =
 *    \left( 1 - \frac{N_{\rm on}}{N_{\rm pred}} \right)
 *    \left( \frac{\partial N_{\gamma}}{\partial p_{\rm sky}} \right)
 * \f]
 *
 * and with respect to background model parameters \f$p_{\rm bgd}\f$ is
 * given by
 *
 * \f[
 *    \left( \frac{\partial \ln L_i}{\partial p_{\rm bgd}} \right) =
 *    \left( 1 + \alpha - \frac{N_{\rm off}}{N_{\rm bgd}} -
 *           \frac{\alpha N_{\rm on}}{N_{\rm pred}} \right)
 *    \left( \frac{\partial N_{\rm bgd}}{\partial p_{\rm bgd}} \right)
 * \f]
 *
 * The curvature matrix elements are given by
 *
 * \f[
 *    \left( \frac{\partial^2 \ln L_i}{\partial^2 p_{\rm sky}} \right) =
 *    \left( \frac{N_{\rm on}}{N_{\rm pred}^2} \right)
 *    \left( \frac{\partial N_{\gamma}}{\partial p_{\rm sky}} \right)^2
 * \f]
 *
 * \f[
 *    \left( \frac{\partial^2 \ln L_i}{\partial p_{\rm sky}
 *                                     \partial p_{\rm bgd}} \right) =
 *    \left( \frac{\alpha N_{\rm on}}{N_{\rm pred}^2} \right)
 *    \left( \frac{\partial N_{\gamma}}{\partial p_{\rm sky}} \right)
 *    \left( \frac{\partial N_{\rm bgd}}{\partial p_{\rm bgd}} \right)
 * \f]
 *
 * \f[
 *    \left( \frac{\partial^2 \ln L_i}{\partial p_{\rm bgd}
 *                                     \partial p_{\rm sky}} \right) =
 *    \left( \frac{\alpha N_{\rm on}}{N_{\rm pred}^2} \right)
 *    \left( \frac{\partial N_{\rm bgd}}{\partial p_{\rm bgd}} \right)
 *    \left( \frac{\partial N_{\gamma}}{\partial p_{\rm sky}} \right)
 * \f]
 *
 * \f[
 *    \left( \frac{\partial^2 \ln L_i}{\partial^2 p_{\rm bgd}} \right) =
 *    \left( \frac{N_{\rm off}}{N_{\rm bgd}^2} +
 *           \frac{\alpha^2 N_{\rm on}}{N_{\rm pred}^2} \right)
 *    \left( \frac{\partial N_{\rm bgd}}{\partial p_{\rm bgd}} \right)^2
 * \f]
 ***********************************************************************/
double GCTAOnOffObservation::likelihood(const GModels& models,
                                        GVector*       gradient,
                                        GMatrixSparse* curvature,
                                        double*        npred) const
{
    // Timing measurement
    #if defined(G_LIKELIHOOD_DEBUG)
    #ifdef _OPENMP
    double t_start = omp_get_wtime();
    #else
    clock_t t_start = clock();
    #endif
    #endif
	
    // Initialise statistics
    #if defined(G_LIKELIHOOD_DEBUG)
    int    n_bins        = m_on_spec.size();
    int    n_used        = 0;
    int    n_small_model = 0;
    int    n_zero_data   = 0;
    double sum_data      = 0.0;
    double sum_model     = 0.0;
    double init_npred    = *npred;
    #endif
	
    // Initialise log-likelihood value
    double value = 0.0;

    // Get number of model parameters in model container
    int npars = models.npars();
    
    // Create model gradient vectors for sky and background parameters
    GVector sky_grad(npars);
    GVector bgd_grad(npars);

    // Allocate working array
    GVector colvar(npars);

    // Get reference to alpha parameters
    const std::vector<double>& alpha = m_arf["ALPHA"];

    // Check that there is at least one parameter
    if (npars > 0) {

        // Loop over all energy bins
        for (int i = 0; i < m_on_spec.size(); ++i) {

            // Reinitialize working arrays
            for (int j = 0; j < npars; ++j) {
                sky_grad[j] = 0.0;
                bgd_grad[j] = 0.0;
            }

            // Get number of On and Off counts
            double non  = m_on_spec[i];
            double noff = m_off_spec[i];

            // Get number of gamma and background events (and corresponding
            // spectral model gradients)
            double ngam = N_gamma(models, i, &sky_grad);
            double nbgd = N_bgd(models, i, &bgd_grad);

            // Skip bin if model is too small (avoids -Inf or NaN gradients)
            double nonpred = ngam + alpha[i] * nbgd;
            if ((nbgd <= minmod) || (nonpred <= minmod)) {
                #if defined(G_LIKELIHOOD_DEBUG)
                n_small_model++;
                #endif
                continue;
            }

            // Now we have all predicted gamma and background events for
            // current energy bin. Update the log(likelihood) and predicted
            // number of events
            value  += -non * log(nonpred) + nonpred - noff * log(nbgd) + nbgd;
            *npred += nonpred;

            // Update statistics
            #if defined(G_LIKELIHOOD_DEBUG)
            n_used++;
            sum_data  += non;
            sum_model += nonpred;
            #endif

            // Fill derivatives
            double fa         = non/nonpred;
            double fb         = fa/nonpred;
            double fc         = alpha[i] * fb;
            double fd         = fc * alpha[i] + noff/(nbgd*nbgd);
            double sky_factor = 1.0 - fa;
            double bgd_factor = 1.0 + alpha[i] - alpha[i] * fa - noff/nbgd;

            // Loop over all parameters
            for (int j = 0; j < npars; ++j) {

                // If spectral model for sky component is non-zero and
                // non-infinite then handle sky component gradients and
                // second derivatives including at least a sky component ...
                if (sky_grad[j] != 0.0  && !gammalib::is_infinite(sky_grad[j])) {

                    // Gradient
                    (*gradient)[j] += sky_factor * sky_grad[j];

                    // Hessian (from first-order derivatives only)
                    for (int k = 0; k < npars; ++k) {

                        // If spectral model for sky component is non-zero and
                        // non-infinite then we have the curvature element
                        // of a sky component
                        if (sky_grad[k] != 0.0  &&
                            !gammalib::is_infinite(sky_grad[k])) {
                            colvar[k] = sky_grad[j] * sky_grad[k] * fb;
                        }

                        // ... else if spectral model for background component
                        // is non-zero and non-infinite then we have the mixed
                        // curvature element between a sky and a background
                        // component
                        else if (bgd_grad[k] != 0.0  &&
                                 !gammalib::is_infinite(bgd_grad[k])) {
                            colvar[k] = sky_grad[j] * bgd_grad[k] * fc;
                        }

                        // ...else neither sky nor background
                        else {
                            colvar[k] = 0.0;
                        }

                    } // endfor: Hessian computation

                    // Update matrix
                    curvature->add_to_column(j, colvar);

                } // endif: spectral model is non-zero and non-infinite

                // ... otherwise if spectral model for background component is
                // non-zero and non-infinite then handle background component
                // gradients and second derivatives including at least a
                // background component
                else if (bgd_grad[j] != 0.0  &&
                         !gammalib::is_infinite(bgd_grad[j])) {

                    // Gradient
                    (*gradient)[j] += bgd_factor * bgd_grad[j];

                    // Hessian (from first-order derivatives only)
                    for (int k = 0; k < npars; ++k) {

                        // If spectral model for sky component is non-zero and
                        // non-infinite then we have the mixed curvature element
                        // between a sky and a background component
                        if (sky_grad[k] != 0.0  &&
                            !gammalib::is_infinite(sky_grad[k])) {
                            colvar[k] = bgd_grad[j] * sky_grad[k] * fc;
                        }

                        // ... else if spectral model for background component
                        // is non-zero and non-infinite then we have the
                        // curvature element of a background component
                        else if (bgd_grad[k] != 0.0  &&
                                 !gammalib::is_infinite(bgd_grad[k])) {
                            colvar[k] = bgd_grad[j] * bgd_grad[k] * fd;
                        }

                        // ... else neither sky nor background
                        else {
                            colvar[k] = 0.0;
                        }

                    } // endfor: Hessian computation

                    // Update matrix
                    curvature->add_to_column(j, colvar);

                } // endif: spectral model for background component valid

            } // endfor: looped over all parameters for derivatives computation

        } // endfor: looped over energy bins

    } // endif: number of parameters was positive
				
	// ... else there are no parameters, so throw an exception
	else {
		std::string msg ="No model parameter for the computation of the "
						 "likelihood in observation \""+this->name()+
		                 "\" (ID "+this->id()+").\n";
		throw GException::invalid_value(G_LIKELIHOOD, msg);
	}
		
    // Dump statistics
    #if defined(G_LIKELIHOOD_DEBUG)
    std::cout << "Number of parameters: " << npars << std::endl;
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
    std::cout << "Number of bins excluded due to small model: ";
    std::cout << n_small_model << std::endl;
    std::cout << "Number of bins with zero data: " << n_zero_data << std::endl;
    std::cout << "Sum of data (On): " << sum_data << std::endl;
    std::cout << "Sum of model (On): " << sum_model << std::endl;
    std::cout << "Statistics: " << value << std::endl;
    #endif
	
    // Optionally dump gradient and curvature matrix
    #if defined(G_LIKELIHOOD_DEBUG)
    std::cout << *gradient << std::endl;
    std::cout << *curvature << std::endl;
    #endif
	
    // Timing measurement
    #if defined(G_LIKELIHOOD_DEBUG)
    #ifdef _OPENMP
    double t_elapse = omp_get_wtime()-t_start;
    #else
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    #endif
    std::cout << "GCTAOnOffObservation::optimizer::likelihood: CPU usage = "
	          << t_elapse << " sec" << std::endl;
    #endif
	
    // Return
    return value;
}


/***********************************************************************//**
 * @brief Print On/Off observation information
 *
 * @return String containing On/Off observation information.
 ***************************************************************************/
std::string GCTAOnOffObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {
    
        // Append header
        result.append("=== GCTAOnOffObservation ===");

        // Append parameters
        result.append("\n"+gammalib::parformat("Name")+m_name);
        result.append("\n"+gammalib::parformat("Identifier")+m_id);

        // Append spectra, ARF and RMF
        result.append("\n"+m_on_spec.print(gammalib::reduce(chatter)));
        result.append("\n"+m_off_spec.print(gammalib::reduce(chatter)));
        result.append("\n"+m_arf.print(gammalib::reduce(chatter)));
        result.append("\n"+m_rmf.print(gammalib::reduce(chatter)));

        // Append regions
        result.append("\n"+m_on_regions.print(gammalib::reduce(chatter)));
        result.append("\n"+m_off_regions.print(gammalib::reduce(chatter)));
    }

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
void GCTAOnOffObservation::init_members(void)
{
    // Initialise members
    m_instrument = "CTAOnOff";
    m_response   = NULL;
    m_ontime     = 0.0;
    m_livetime   = 0.0;
    m_deadc      = 1.0;
    m_on_spec.clear();
    m_off_spec.clear();
    m_arf.clear();
    m_rmf.clear();
    m_on_regions.clear();
    m_off_regions.clear();
    m_src_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs CTA On/Off observation.
 ***************************************************************************/
void GCTAOnOffObservation::copy_members(const GCTAOnOffObservation& obs)
{
    // Copy attributes
    m_instrument  = obs.m_instrument;
    m_ontime      = obs.m_ontime;
    m_livetime    = obs.m_livetime;
    m_deadc       = obs.m_deadc;
    m_on_spec     = obs.m_on_spec;
    m_off_spec    = obs.m_off_spec;
    m_arf         = obs.m_arf;
    m_rmf         = obs.m_rmf;
    m_on_regions  = obs.m_on_regions;
    m_off_regions = obs.m_off_regions;
    m_src_dir     = obs.m_src_dir;

    // Clone members
    m_response = (obs.m_response != NULL) ? obs.m_response->clone() : NULL;

    // Return
    return;
}



/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAOnOffObservation::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set On/Off observation from a CTA observation
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No CTA event list found in CTA observation.
 *
 * Sets an On/Off observation from a CTA observation by filling the events
 * that fall in the On and Off regions into the PHA spectra and by computing
 * the corresponding ARF and RMF response functions.
 ***************************************************************************/
void GCTAOnOffObservation::set(const GCTAObservation& obs)
{
    // Get CTA event list pointer
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
    if (events == NULL) {
        std::string msg = "No event list found in CTA observation \""+
                          obs.name()+"\" (ID="+obs.id()+"). ON/OFF observation "
                          "can only be filled from event list.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Loop over all events
    for (int i = 0; i < events->size(); ++i) {

        // Get measured event direction
        const GCTAEventAtom* atom = (*events)[i];
        GSkyDir              dir  = atom->dir().dir();

        // Fill in spectrum according to region containment
        if (m_on_regions.contains(dir)) {
            m_on_spec.fill(atom->energy());
        }
        if (m_off_regions.contains(dir)) {
            m_off_spec.fill(atom->energy());
        }
        
    } // endfor: looped over all events

    // Store the livetime as exposures of the spectra
    m_on_spec.exposure(obs.livetime());
    m_off_spec.exposure(obs.livetime());

    // Store the ontime, livetime and deadtime correction in the observation
    m_ontime   = obs.ontime();
    m_livetime = obs.livetime();
    m_deadc    = obs.deadc();

    // Compute response components
    compute_arf(obs);
    compute_bgd(obs);
    compute_alpha(obs);
    compute_rmf(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute ARF of On/Off observation
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No CTA response found in CTA observation.
 *
 * @todo Implement GCTAResponse::npred usage instead of computing Aeff
 *       (required the integration of PSF over On region).
 ***************************************************************************/
void GCTAOnOffObservation::compute_arf(const GCTAObservation& obs)
{
    // Get reconstructed energy boundaries from on ARF
    GEbounds etrue = m_arf.ebounds();
    int      ntrue = etrue.size();

    // Continue only if there are ARF bins
    if (ntrue > 0) {
    
        // Get CTA response pointer. Throw an exception if no response is found
        const GCTAResponseIrf* response =
              dynamic_cast<const GCTAResponseIrf*>(obs.response());
        if (response == NULL) {
            std::string msg = "Response in CTA observation \""+obs.name()+"\" "
                              "(ID="+obs.id()+") is not of the GCTAResponseIrf "
                              "type.";
            throw GException::invalid_value(G_COMPUTE_ARF, msg);
        }

        // Get CTA observation pointing direction, zenith, and azimuth
        GCTAPointing obspnt  = obs.pointing();
        GSkyDir      obsdir  = obspnt.dir();
        double       zenith  = obspnt.zenith();
        double       azimuth = obspnt.azimuth();

        // Retrieve ON map pixel indices and declare working variables
        std::vector<int> pixidx = m_on_regions.nonzeroindices();
        int pixnum = pixidx.size();
        
        // Loop over true energies
        for (int i = 0; i < ntrue; ++i) {
        
            // Get mean energy of bin
            double logEtrue = etrue.elogmean(i).log10TeV();

            // Initialize effective area for this bin
            m_arf[i] = 0.0;

            // Initialize totals
            double totsolid = 0.0;
            double totpsf = 0.0;

            // Loop over pixels in ON region map and integrate effective area 
            for (int j = 0; j < pixnum; j++ ) {
                
                // Get direction to pixel center
                GSkyDir pixdir=m_on_regions.map().inx2dir(pixidx[j]);
                // Get solid angle subtended by this pixel
                double pixsolid=m_on_regions.map().solidangle(pixidx[j]);
                
                // Compute position of pixel centre in instrument coordinates
                double theta = obsdir.dist(pixdir);
                double phi   = obsdir.posang(pixdir);
                
                // Add up effective area
                m_arf[i] += response->aeff(theta,
                                           phi,
                                           zenith,
                                           azimuth,
                                           logEtrue)*pixsolid;
                                           
                // Add pixel solid angle to total for averaging later
                totsolid += pixsolid;
                
                // Integrate PSF
                double delta = m_src_dir.dist_deg(pixdir);
                totpsf += response->psf(delta,
                                        theta,
                                        phi,
                                        zenith,
                                        azimuth,
                                        logEtrue)*pixsolid;

            } // Looped over all pixels in map
            
            // Average effective area over solid angle
            if (totsolid > 0.0) {
                m_arf[i] /= totsolid;
            }
            
            // Correct effective area by containment fraction
            if (totpsf >= 0.0) {
                m_arf[i] *= totpsf;
            }
       
        } // endfor: looped over true energies
        
	} // endif: there were energy bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute background rate in Off region
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_argument
 *            Observation does not contain relevant response or background
 *            information
 *
 * Compute the background rate in units of events/s/MeV in all Off regions
 * and stores the background rate as additional column with name
 * `BACKGROUND` in the Auxiliary Response File (ARF).
 *
 ***************************************************************************/
void GCTAOnOffObservation::compute_bgd(const GCTAObservation& obs)
{	
    // Get reconstructed energy boundaries from on ARF
	GEbounds ereco = m_arf.ebounds();
    int      nreco = ereco.size();
    
    // Continue only if there are ARF bins
    if (nreco > 0) {

		// Initialise background rates to zero
        std::vector<double> background(nreco, 0.0);
    		
		// Get CTA observation pointing direction
		GCTAPointing obspnt = obs.pointing();
	
        // Get pointer on CTA IRF response
        const GCTAResponseIrf* rsp =
              dynamic_cast<const GCTAResponseIrf*>(obs.response());
        if (rsp == NULL) {
            std::string msg = "Specified observation does not contain an "
                              "IRF response.\n" + obs.print();
            throw GException::invalid_argument(G_COMPUTE_BGD, msg);
        }

        // Get pointer to CTA background
        const GCTABackground* bgd = rsp->background();
        if (bgd == NULL) {
            std::string msg = "Specified observation contains no "
                              "background information.\n" + obs.print();
            throw GException::invalid_argument(G_COMPUTE_BGD, msg);
        }
        
        // Loop over pixels in OFF region map and integrate background rate 
        std::vector<int> pixidx= m_off_regions.nonzeroindices();
        int pixnum=pixidx.size();
        for (int j = 0; j < pixnum; j++ ) {
            
            // Get direction to pixel center
            GSkyDir pixdir=m_off_regions.map().inx2dir(pixidx[j]);
            // Translate sky direction into instrument direction
            GCTAInstDir pixinstdir = obspnt.instdir(pixdir);
            // Get solid angle subtended by this pixel
            double pixsolid = m_off_regions.map().solidangle(pixidx[j]);
            
            // Loop over energy bins
            for (int i = 0; i < nreco; ++i) {

                // Get log10(E/TeV) of mean reconstructed bin energy
                double logEreco = m_on_spec.ebounds().elogmean(i).log10TeV();

                // Get background rate in events/s/MeV
                background[i] += (*bgd)(logEreco,
                                        pixinstdir.detx(),
                                        pixinstdir.dety()) * pixsolid;
            } // Looped over energy bins
        
        } // Looped over all pixels in map
        
        // Append background vector to ARF
        m_arf.append("BACKGROUND", background);

    } // endif: there were spectral bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute vector of alpha parameters
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No CTA response found in CTA observation.
 *
 * Compute the alpha parameters for all energy bins. The alpha parameter
 * gives the ratio between the On and Off region background acceptance
 * multiplied by the ratio between On and Off region solid angles.
 *
 ***************************************************************************/
void GCTAOnOffObservation::compute_alpha(const GCTAObservation& obs)
{
    // Get reconstructed energy boundaries from on ARF
    GEbounds ereco = m_rmf.emeasured();
    int      nreco = ereco.size();
    
    // Continue only if there are ARF bins
    if (nreco > 0) {

        // Initialise On/Off exposure ratios
        std::vector<double> alpha(nreco, 0.0);
    
        // Get CTA response pointer. Throw an exception if no response is found
        const GCTAResponseIrf* response =
              dynamic_cast<const GCTAResponseIrf*>(obs.response());
        if (response == NULL) {
            std::string msg = "Response in CTA observation \""+obs.name()+"\" "
                              "(ID="+obs.id()+") is not of the GCTAResponseIrf "
                              "type.";
            throw GException::invalid_value(G_COMPUTE_ALPHA, msg);
        }

        // Get CTA observation pointing direction, zenith, and azimuth
        GCTAPointing obspnt  = obs.pointing();
        GSkyDir      obsdir  = obspnt.dir();

        // Declare variables for map handling
        std::vector<int> pixidx;
        int pixnum = 0;

        // Loop over reconstructed energies 
        for (int i = 0; i < nreco; ++i) {
        
            // Get mean log10 energy in TeV of bin
            double logEreco = ereco.elogmean(i).log10TeV();

            // Initialise background rate totals
            double aon(0.0);
            double aoff(0.0);
            
            // Loop over pixels in ON region map and integrate acceptance 
            pixidx = m_on_regions.nonzeroindices();
            pixnum = pixidx.size();
            for (int j = 0; j < pixnum; j++ ) {
                
                // Get direction to pixel center
                GSkyDir pixdir = m_on_regions.map().inx2dir(pixidx[j]);
                // Translate sky direction into instrument direction
                GCTAInstDir pixinstdir = obspnt.instdir(pixdir);
                // Get solid angle subtended by this pixel
                double pixsolid = m_on_regions.map().solidangle(pixidx[j]);
             
                // Add up acceptance
                aon += (*response->background())(logEreco,
                                                 pixinstdir.detx(),
                                                 pixinstdir.dety()) * pixsolid;
            } // Looped over all pixels in map

            // Loop over pixels in OFF region map and integrate acceptance 
            pixidx = m_off_regions.nonzeroindices();
            pixnum = pixidx.size();
            for (int j = 0; j < pixnum; j++ ) {
            
                // Get direction to pixel center
                GSkyDir pixdir = m_off_regions.map().inx2dir(pixidx[j]);
                // Translate sky direction into instrument direction
                GCTAInstDir pixinstdir = obspnt.instdir(pixdir);
                // Get solid angle subtended by this pixel
                double pixsolid = m_off_regions.map().solidangle(pixidx[j]);
            
                // Add up acceptance
                aoff += (*response->background())(logEreco,
                                                 pixinstdir.detx(),
                                                 pixinstdir.dety()) * pixsolid;
            } // Looped over all pixels in map
            
            // Compute alpha for this energy bin
            if (aoff > 0.0) {
                alpha[i] = aon/aoff;
            }

        } // endfor: looped over reconstructed energies

        // Append alpha vector to ARF
        m_arf.append("ALPHA", alpha);

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute RMF of On/Off observation
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            Observation does not contain IRF response
 *
 * Compute the energy redistribution matrix for an On/Off observation. The
 * method requires that the RMF energy axes have been defined before.
 ***************************************************************************/
void GCTAOnOffObservation::compute_rmf(const GCTAObservation& obs)
{
    // Get true and reconstructed energy boundaries from RMF
    GEbounds etrue = m_rmf.etrue();
    GEbounds ereco = m_rmf.emeasured();
    int      ntrue = etrue.size();
    int      nreco = ereco.size();

    // Continue only if there are RMF bins
    if (ntrue > 0 && nreco > 0) {
    
        // Get CTA response pointer
        const GCTAResponseIrf* response =
              dynamic_cast<const GCTAResponseIrf*>(obs.response());
        if (response == NULL) {
            std::string msg = "Response in CTA observation \""+obs.name()+"\" "
                              "(ID="+obs.id()+") is not of the GCTAResponseIrf "
                              "type.";
            throw GException::invalid_value(G_COMPUTE_RMF, msg);
        }

        // Get CTA observation pointing direction, zenith, and azimuth
        GCTAPointing obspnt  = obs.pointing();
        GSkyDir      obsdir  = obspnt.dir();
        double       zenith  = obspnt.zenith();
        double       azimuth = obspnt.azimuth();

        // Retrieve map information
        std::vector<int> pixidx= m_on_regions.nonzeroindices();
        int pixnum=pixidx.size();
        
        // Loop over reconstructed energy
        for (int ireco = 0; ireco < nreco; ++ireco) {

            // Get log10(E)
            double logEreco = ereco.elogmean(ireco).log10TeV();

            // Loop over true energy
            for (int itrue = 0; itrue < ntrue; ++itrue) {

                // Initialise RMF element
                m_rmf(itrue, ireco) = 0.0;
                
                // Compute true energy and initialise weight for this bin
                double logEtrue = etrue.elogmean(itrue).log10TeV();
                double weight   = 0.0;

                // Loop over pixels in ON region map and integrate acceptance
                for (int j = 0; j < pixnum; j++ ) {
                
                    // Get direction to pixel center
                    GSkyDir pixdir = m_on_regions.map().inx2dir(pixidx[j]);
                    // Translate sky direction into instrument direction
                    GCTAInstDir pixinstdir = obspnt.instdir(pixdir);
                    // Get solid angle subtended by this pixel
                    double pixsolid=m_on_regions.map().solidangle(pixidx[j]);
             
                    // Compute position of pixel centre in instrument coordinates
                    double theta = obsdir.dist(pixdir);
                    double phi   = obsdir.posang(pixdir);

                    // Get effective area
                    double aeff = response->aeff(theta,
                                                 phi,
                                                 zenith,
                                                 azimuth,
                                                 logEreco);

                    // Add up energy dispersion weighted by effective area
                    weight              += aeff;
                    m_rmf(itrue, ireco) += response->edisp(ereco.elogmean(ireco),
                                                           theta,
                                                           phi,
                                                           zenith,
                                                           azimuth,
                                                           logEtrue) * aeff;
                } // Looped over all pixels in map
                
                // Complete weighing by dividing by total effective area
                if (weight > 0.0) {
                    m_rmf(itrue, ireco) /= weight;
                }	

            } // endfor: looped over true energy

        } // endfor: looped over reconstructed energy

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************
 * @brief Compute \f$N_{\gamma}\f$ value and model parameter gradients
 *
 * @param[in] models Model container.
 * @param[in] ibin Energy bin number.
 * @param[in,out] grad Model gradient vector.
 *
 * Returns the predicted number of source events \f$N_{\gamma}\f$
 * in the On regions for a given energy bin. The method computes also
 *
 * \f[
 *    \left( \frac{\partial N_{\gamma}}{\partial p_{\rm sky}} \right)
 * \f]
 *
 * which are the gradients in the predicted number of source events with
 * respect to all model parameters.
 *
 * The method assumes that parameters are stored in the order
 * spatial-spectral-temporal.
 *
 * @todo I think this method only works for point sources. What happens
 *       for an extended source?
 ***********************************************************************/
double GCTAOnOffObservation::N_gamma(const GModels& models,
                                     const int&     ibin,
                                     GVector*       grad) const
{
    // Get total number of model parameters
    int npars = models.npars();

    // Initialize results
    double value = 0.0;
    for (int i = 0; i < npars; ++i) {
        (*grad)[i] = 0.0;
    }
    
    // Continue only if bin number is in range and there are model parameters
    if ((ibin >= 0) && (ibin < m_on_spec.size()) && (npars > 0)) {

        // Initialise parameter index
        int ipar = 0;
        
        // Get true energy binning (loop over bins further down)
        GEbounds etrue = m_arf.ebounds();
        int      ntrue = etrue.size();
        
        /* Previous version without RMF

        // Get reconstructed energy bin properties
        const GEnergy erecomin   = m_on_spec.ebounds().emin(ibin);
        const GEnergy erecomax   = m_on_spec.ebounds().emax(ibin);
        const GEnergy erecomean  = m_on_spec.ebounds().elogmean(ibin);
        const double  erecowidth = m_on_spec.ebounds().ewidth(ibin).MeV();
        
        // Compute normalisation factors
        double exposure  = m_on_spec.exposure();
        double norm_flux = m_arf[ibin] * exposure; // cm2 s
        double norm_grad = norm_flux   * ewidth;   // cm2 s MeV
        
        */

        // Loop over models
        for (int j = 0; j < models.size(); ++j) {

            // Get model pointer. Fall through if pointer is not valid
            const GModel* mptr = models[j];
            if (mptr == NULL) {
                continue;
            }

            // Fall through if model does not apply to specific instrument
            // and observation identifier
            if (!mptr->is_valid(instrument(), id())) {
                ipar += mptr->size();
                continue;
            }

            // Fall through if this model component is a not sky component
            const GModelSky* sky = dynamic_cast<const GModelSky*>(mptr);
            if (sky == NULL) {
                ipar += mptr->size();
                continue;
            }

            // Increase parameter counter for spatial parameter
            GModelSpatial* spatial = sky->spatial();
            if (spatial != NULL)  {
                ipar += spatial->size();
            }

            // Spectral component (the useful one)
            GModelSpectral* spectral = sky->spectral();
            if (spectral != NULL)  {
                  
                    // Loop over true energy bins
                    for (int itrue = 0; itrue < ntrue; ++itrue) {
                        
                        // True energy bin properties
                        const GEnergy etruemin   = etrue.emin(itrue);
                        const GEnergy etruemax   = etrue.emax(itrue);
                        const GEnergy etruemean  = etrue.elogmean(itrue);
                        const double  etruewidth = etrue.ewidth(itrue).MeV();
                        
                        // Determine number of gamma-ray events in model by
                        // computing the flux over the true  energy bin 
                        // in ph/cm2/s and multiplying this by effective area (cm2)
                        // and livetime (s) and redistribution probability
                        double norm_flux = m_arf[itrue] * m_on_spec.exposure() * m_rmf(itrue, ibin);
                        value += spectral->flux(etruemin, etruemax) * norm_flux;
                        
                        // Determine the model gradients at the current true energy. The
                        // eval() method needs a time in case that the spectral model
                        // has a time dependence. We simply use a dummy time here.
                        double norm_grad = norm_flux * etruewidth;
                        spectral->eval(etruemean, GTime(), true);

                        // Loop over spectral model parameters
                        for (int k = 0; k < spectral->size(); ++k, ++ipar)  {
                            GModelPar& par = (*spectral)[k];
                            if (par.is_free() && ipar < npars)  {
                                (*grad)[ipar] += par.factor_gradient() * norm_grad;
                            }
                        } // Looped over model parameters
                        
                    } // Looped over true energy bins
                
                
                /* Previous version without RMF
                
                // Determine the number of gamma-ray events in model by
                // computing the flux over the energy bin in ph/cm2/s
                // and multiplying this flux by the effective area (cm2)
                // and the livetime (s)
                value += spectral->flux(emin, emax) * norm_flux;

                // Determine the model gradients at the current energy. The
                // eval() method needs a time in case that the spectral model
                // has a time dependence. We simply use a dummy time here.
                spectral->eval(emean, GTime(), true);

                // Loop over spectral model parameters
                for (int k = 0; k < spectral->size(); ++k, ++ipar)  {
                    GModelPar& par = (*spectral)[k];
                    if (par.is_free() && ipar < npars)  {
                        (*grad)[ipar] = par.factor_gradient() * norm_grad;
                    }
                }
                
                */

            } // endif: spectral component

            // Increase parameter counter for temporal parameter
            GModelTemporal* temporal = sky->temporal();
            if (temporal != NULL)  {
                ipar += temporal->size();
            }

        } // endfor: looped over model components

    } // endif: bin number is in the range and model container is not empty

    // Return number of gamma-ray events
    return value;
}


/***********************************************************************
 * @brief Compute \f$N_{\rm bgd}\f$ value and model parameter gradients
 *
 * @param[in] models Model container.
 * @param[in] ibin Energy bin index.
 * @param[in,out] grad Model gradient vector.
 * @return Predicted number of background events in Off regions.
 *
 * Returns the predicted number of background events \f$N_{\rm bgd}\f$
 * in the Off regions for a given energy bin. The method computes also
 *
 * \f[
 *    \left( \frac{\partial N_{\rm bgd}}{\partial p_{\rm bgd}} \right)
 * \f]
 *
 * which are the gradients in the predicted number of background events
 * with respect to all model parameters.
 *
 * The method assumes that the model parameters are stored in the order
 * spectral-temporal.
 *
 * @todo So far only handles the GCTAModelIrfBackground background model;
 *       Should be more generic, but there is maybe for the moment no
 *       other choice than implementing something for all possible
 *       background models; that's not very satisfactory ...
 ***********************************************************************/
double GCTAOnOffObservation::N_bgd(const GModels& models,
                                   const int&     ibin,
                                   GVector*       grad) const
{
    // Get total number of model parameters
    int npars = models.npars();
    
    // Initialize results
    double value = 0.0;
    for (int i = 0; i < npars; ++i) {
        (*grad)[i] = 0.0;
    }

    // Continue only if bin number is valid and if there are model parameters
    if ((ibin >= 0) && (ibin < m_off_spec.size()) && (npars > 0))  {

        // Initialise parameter index
        int ipar = 0;

        // Get reference to background rates (events/MeV/s)
        const std::vector<double>& background = m_arf["BACKGROUND"];

        // Get reconstructed energy bin mean and width
        const GEnergy emean  = m_off_spec.ebounds().elogmean(ibin);
        const double  ewidth = m_off_spec.ebounds().ewidth(ibin).MeV();

        // Compute normalisation factor (events)
        double exposure = m_off_spec.exposure();
        double norm     = background[ibin] * exposure * ewidth;

        // Loop over models
        for (int j = 0; j < models.size(); ++j) {

            // Get model pointer. Fall through if pointer is not valid
            const GModel* mptr = models[j];
            if (mptr == NULL) {
                continue;
            }

            // Fall through if model does not apply to specific instrument
            // and observation identifier.
            if (!mptr->is_valid(this->instrument(), this->id())) {
                ipar += mptr->size();
                continue;
            }

            // Fall through if model is not an IRF background component
            const GCTAModelIrfBackground* bgd =
                  dynamic_cast<const GCTAModelIrfBackground*>(mptr);
            if (bgd == NULL) {
                ipar += mptr->size();
                continue;
            }

            // Get spectral component
            GModelSpectral* spectral = bgd->spectral();
            if (spectral != NULL)  {

                // Determine the number of background events in model by
                // computing the model normalization at the mean bin energy bin
                // and multiplying the normalisation with the number of
                // background events. The eval() method needs a time in case
                // that the spectral model has a time dependence. We simply
                // use a dummy time here.
                value += spectral->eval(emean, GTime(), true) * norm;

                // Compute the parameter gradients for all spectral model
                // parameters
                for (int k = 0; k < spectral->size(); ++k, ++ipar)  {
                    GModelPar& par = (*spectral)[k];
                    if (par.is_free() && ipar < npars)  {
                        (*grad)[ipar] += par.factor_gradient() * norm;
                    }
                }

            } // endif: pointer to spectral component was not NULL

            // Increase parameter counter for temporal parameter
            GModelTemporal* temporal = bgd->temporal();
            if (temporal != NULL)  {
                ipar += temporal->size();
            }

        } // endfor: looped over model components

    } // endif: bin number is in the range and model container is not empty

    // Return
    return value;
}
