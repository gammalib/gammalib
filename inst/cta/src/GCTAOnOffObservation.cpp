/***************************************************************************
 *          GCTAOnOffObservation.cpp - CTA on-off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Chia-Chun Lu & Christoph Deil               *
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
 * @brief CTA on-off observation class implementation
 * @author Chia-Chun Lu & Christoph Deil
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
#include "GSkyRegionCircle.hpp"
#include "GOptimizerPars.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTACubeBackground.hpp"
#include "GCTAModelIrfBackground.hpp"
#include "GCTAOnOffObservation.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAOnOffObservation g_onoff_obs_cta_seed;
const GObservationRegistry g_onoff_obs_cta_registry(&g_onoff_obs_cta_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE_SET           "GCTAOnOffObservation::response(GResponse&)"
#define G_RESPONSE_GET                     "GCTAOnOffObservation::response()"
#define G_WRITE                   "GCTAOnOffObservation::write(GXmlElement&)"
#define G_READ                     "GCTAOnOffObservation::read(GXmlElement&)"
#define G_FILL                 "GCTAOnOffObservation::fill(GCTAObservation&)"
#define G_COMPUTE_ALPHA                "GCTAOnOffObservation::compute_alpha("\
                                                          "GCTAObservation&)"
#define G_COMPUTE_ARF   "GCTAOnOffObservation::compute_arf(GCTAObservation&)"
#define G_COMPUTE_BGD   "GCTAOnOffObservation::compute_bgd(GCTAObservation&,"\
                                                                  "GModels&)"
#define G_COMPUTE_RMF   "GCTAOnOffObservation::compute_rmf(GCTAObservation&,"\
                                                                " GEbounds&)"
#define G_POISSON_ONOFF     "GCTAOnOffObservation::likelihood_poisson_onoff("\
                                 "GModels&, GOptimizerPars&, GMatrixSparse&,"\
                                               " GVector&, double&, double&)"

/* __ Constants __________________________________________________________ */
const double minmod = 1.0e-100;                      //!< Minimum model value
const double minerr = 1.0e-100;                //!< Minimum statistical error

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs CTA ON/OFF observation.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GCTAOnOffObservation& obs)
{ 
    // Initialise private
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] ereco Reconstructed energy bins.
 * @param[in] on ON regions.
 * @param[in] off OFF regions.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GEbounds&    ereco,
                                           const GSkyRegions& on,
                                           const GSkyRegions& off)
{
    // Initialise private
    init_members();

    // Initialise spectra
    m_on_spec  = GPha(ereco);
    m_off_spec = GPha(ereco);

    // Store regions
    m_on_regions  = on;
    m_off_regions = off;

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
 * @param[in] obs CTA ON/OFF observation.
 * @return CTA ON/OFF observation.
 *
 * Assigns one CTA ON/OFF observation to another ON/OFF observation object.
 ***************************************************************************/
GCTAOnOffObservation& GCTAOnOffObservation::operator=(const GCTAOnOffObservation& obs)
{ 
    // Execute only if object is not identical
    if (this != &obs) {

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
 * Clears the ON/OFF observation. All class members will be set to the
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
 * @return Pointer to deep copy of ON/OFF observation.
 *
 * Returns a pointer to a deep copy of an ON/OFF observation.
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
 * Sets the response function for the observation.
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
 * @brief Read ON/OFFobservation from an xml element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::xml_invalid_parnum
 *            Invalid number of parameters found in XML element.
 * @exception GException::xml_invalid_parnames
 *            Invalid parameter names found in XML element.
 *
 * Reads information for a CTA ON/OFF observation from an XML element.
 * The expected format of the XML element is
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

	// Determine number of parameter nodes in XML element
	int npars = xml.elements("parameter");

	// Verify that XML element has exactly 6 parameters
	if (xml.elements() != 6 || npars != 6) {
		throw GException::xml_invalid_parnum(G_READ, xml,
			  "CTA observation requires exactly 6 parameters.");
	}

	// Extract parameters
	int npar[] = {0, 0, 0, 0, 0, 0};
	for (int i = 0; i < npars; ++i) {

		// Get parameter element
		const GXmlElement* par = xml.element("parameter", i);

		// Handle Pha_on
		if (par->attribute("name") == "Pha_on") {

			// Read pha_on file name
			std::string filename = par->attribute("file");

			// Load pha_on
			m_on_spec.load(filename);

			// Increment parameter counter
			npar[0]++;
		}

		// Handle Pha_off
		else if (par->attribute("name") == "Pha_off") {

			// Read countsmap file name
			std::string filename = par->attribute("file");

			// Load pha_off
			m_off_spec.load(filename);

			// Increment parameter counter
			npar[1]++;
		}

		// Handle on regions
		else if(par->attribute("name") == "Regions_on") {

			// Get filename
			std::string filename = par->attribute("file");

			// load on regions
			m_on_regions.load(filename);

			// Increase number of parameters
			npar[2]++;
		}

		// Handle off regions
		else if(par->attribute("name") == "Regions_off") {

			// Get filename
			std::string filename = par->attribute("file");

			// load off regions
			m_on_regions.load(filename);

			// Increase number of parameters
			npar[3]++;
		}

		// Handle Arf
		else if(par->attribute("name") == "Arf") {

			// Get filename
			std::string filename = par->attribute("file");

			// load arf
			m_arf.load(filename);

			// Increase number of parameters
			npar[4]++;
		}

		// Handle Rmf
		else if(par->attribute("name") == "Rmf") {

			// Get filename
			std::string filename = par->attribute("file");

			// load Rmf
			m_rmf.load(filename);

			// Increase number of parameters
			npar[5]++;
		}

	} // endfor: looped over all parameters

	// Verify that all parameters were found
	if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 ||
        npar[3] != 1 || npar[4] != 1 || npar[5] != 1) {
		throw GException::xml_invalid_parnames(G_READ, xml,
			  "Require \"Pha_on\" or \"Pha_off\" and \"Regions_on\""
			  ", \"Regions_off\",\"Arf\" and \"Rmf\" parameters.");
	}

	// Return
	return;
}


/***********************************************************************//**
 * @brief write observation to an xml element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::xml_invalid_parnum
 *            Invalid number of parameters found in XML element.
 * @exception GException::xml_invalid_parnames
 *            Invalid parameter names found in XML element.
 *
 * Writes information for a CTA ON/OFF observation into an XML element. The
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
void GCTAOnOffObservation::write(GXmlElement& xml) const
{
	// If XML element has 0 nodes then append 6 parameter nodes
	if (xml.elements() == 0) {
		xml.append(GXmlElement("parameter name=\"Pha_on\""));
		xml.append(GXmlElement("parameter name=\"Pha_off\""));
		xml.append(GXmlElement("parameter name=\"Regions_on\""));
		xml.append(GXmlElement("parameter name=\"Regions_off\""));
		xml.append(GXmlElement("parameter name=\"Arf\""));
		xml.append(GXmlElement("parameter name=\"Rmf\""));
	}

	// Verify that XML element has exactly 4 parameters
	if (xml.elements() != 6 || xml.elements("parameter") != 6) {
		throw GException::xml_invalid_parnum(G_WRITE, xml,
			  "CTAOnOffObservation requires exactly 6 parameters.");
	}

	// Set or update parameter attributes
	int npar[] = {0, 0, 0, 0, 0, 0};
	for (int i = 0; i < 6; ++i) {

		// Get parameter element
		GXmlElement* par = xml.element("parameter", i);

		// Handle on counts
		if (par->attribute("name") == "Pha_on") {
			par->attribute("file", m_on_spec.filename());
			npar[0]++;
		}

		// Handle off counts
		else if (par->attribute("name") == "Pha_off") {
			par->attribute("file", m_off_spec.filename());
			npar[1]++;
		}

		// Handle on regions
		else if (par->attribute("name") == "Regions_on") {
			par->attribute("file", m_on_regions.filename());
			npar[2]++;
		}

		// Handle off regions
		else if (par->attribute("name") == "Regions_off") {
			par->attribute("file", m_off_regions.filename());
			npar[3]++;
		}

		// Handle effective area
		else if (par->attribute("name") == "Arf") {
			par->attribute("file", m_arf.filename());
			npar[4]++;
		}

		// Handle energy resolution
		else if (par->attribute("name") == "Rmf") {
			par->attribute("file", m_rmf.filename());
			npar[5]++;
		}
	}

	// Verify that all required parameters are present
	if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 ||
        npar[3] != 1 || npar[4] != 1 || npar[5] != 1) {
		throw GException::xml_invalid_parnames(G_WRITE, xml,
			  "Require \"Pha_on\" or \"Pha_off\" and \"Regions_on\","
              " \"Regions_off\", \"Arf\", and \"Rmf\" parameters.");
	}

	// Return
	return;
}


/***********************************************************************//**
 * @brief Fill events in the On and Off spectra
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No CTA event list found in CTA observation.
 ***************************************************************************/
void GCTAOnOffObservation::fill(const GCTAObservation& obs)
{
    // Get CTA event list pointer
	const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
	if (events == NULL) {
        std::string msg = "No event list found in CTA observation \""+
                          obs.name()+"\" (ID="+obs.id()+"). ON/OFF observation "
                          "can only be filled from event list.";
        throw GException::invalid_value(G_FILL, msg);
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
	
	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute response for On/Off observation
 *
 * @param[in] obs CTA observation.
 * @param[in] etrue True energy boundaries.
 ***************************************************************************/
void GCTAOnOffObservation::compute_response(const GCTAObservation& obs,
                                            const GModels&         models,
                                            const GEbounds&        etrue)
{
	// Compute response components
	compute_arf(obs);
	compute_bgd(obs, models);
	compute_rmf(obs, etrue);
	compute_alpha(obs);

	// Return
	return;
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


/***********************************************************************
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
 * derivatives and the curvature matrix. The number of On and Off counts
 * are taken from the On and Off spectra, the expected number of
 * gamma-ray and background events are computed from the spectral models
 * of the relevant components in the model container (spatial and
 * temporal components are ignored so far). See the model_on() and
 * model_off() for details about the model computations.
 *
 * @todo Add formula for curvature matrix computation
 ***********************************************************************/
double GCTAOnOffObservation::likelihood(const GModels& models,
                                        GVector*       gradient,
                                        GMatrixSparse* curvature,
                                        double*        npred) const
{
    // Timing measurement
    #if G_EVAL_TIMING
    clock_t t_start = clock();
    #endif
	
    // Initialise statistics
    #if G_OPT_DEBUG
	int    n_bins        = m_on_spec.size();
    int    n_used        = 0;
    int    n_small_model = 0;
    int    n_zero_data   = 0;
    double sum_data      = 0.0;
    double sum_model     = 0.0;
	double init_npred    = *npred;
    #endif
	
	// Initialise likelihood value
    double value = 0.0;
	
	// Get number of parameters
    int npars = models.size();
	
	// Create model gradients array pointers (one for sky parameters, one for
    // background parameters)
	GVector sky_grad(npars);
	GVector bgd_grad(npars);

	// Working arrays
	GVector colvar(npars);
	
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
            double ngam = model_on(models, i, &sky_grad);
            double nbgd = model_off(models, i, &bgd_grad);

            // Skip bin if model is too small (avoids -Inf or NaN gradients)
            double nonpred = ngam + m_alpha[i] * nbgd;
            if ((nbgd <= minmod) || (nonpred <= minmod)) {
                #if G_OPT_DEBUG
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
            #if G_OPT_DEBUG
            n_used++;
            sum_data  += non;
            sum_model += nonpred;
            #endif
		  
            // Fill derivatives
            double fa         = non/nonpred;
            double fb         = fa/nonpred;
            double fc         = m_alpha[i] * fb;
            double fd         = fc * m_alpha[i] + noff/(nbgd*nbgd);
            double sky_factor = 1.0 - fa;
            double bgd_factor = 1.0 + m_alpha[i] - m_alpha[i] * fa - noff/nbgd;

            // Loop over all parameters
            for (int j = 0; j < npars; ++j) {
			
                // If spectral model for sky component is non-zero and non-infinite
                if (sky_grad[j] != 0.0  && !gammalib::is_infinite(sky_grad[j])) {
				
                    // Gradient
                    (*gradient)[j] += sky_factor * sky_grad[j];
				
                    // Hessian (from first-order derivatives only)
                    for (int k = 0; k < npars; ++k) {
				
                        // If spectral model for sky component is non-zero and
                        // non-infinite ...
                        if (sky_grad[k] != 0.0  && !gammalib::is_infinite(sky_grad[k])) {
                            colvar[k] = sky_grad[j] * sky_grad[k] * fb;
                        }
					
                        // ... else if spectral model for sky component is
                        // non-zero and non-infinite
                        // ???????????????????? sky_grad[j]=0 ?????????????????
                        else if (bgd_grad[k] != 0.0  && !gammalib::is_infinite(bgd_grad[k])) {
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
			
                // ... else if spectral model for background component is
                // non-zero and non-infinite
                else if (bgd_grad[j] != 0.0  && !gammalib::is_infinite(bgd_grad[j])) {
				
                    // Gradient
                    (*gradient)[j] += bgd_factor * bgd_grad[j];
				
                    // Hessian (from first-order derivatives only)
                    for (int k = 0; k < npars; ++k) {
					
                        // If spectral model for sky component is non-zero and non-infinite
                        // ???????????????????? sky_grad[j] ?????????????????
                        if (sky_grad[k] != 0.0  && !gammalib::is_infinite(sky_grad[k])) {
                            colvar[k] = bgd_grad[j] * bgd_grad[k] * fc;
                        }
						
						// ... else if spectral model for sky component is
                        // non-zero and non-infinite
                        else if (bgd_grad[k] != 0.0  && !gammalib::is_infinite(bgd_grad[k])) {
                            colvar[k] = bgd_grad[j] * bgd_grad[k] * fd;
						
						// ... else neither sky nor background
                        } else {
                            colvar[k] = 0.0;
                        }

                    } // endfor: Hessian computation
				
                    // Update matrix
                    curvature->add_to_column(j, colvar);
			
                } // endif: spectral model for background component is non-
                  //        zero and non-infinite
						
            } // endfor: looped over all parameters for derivatives computation
			
        } // endfor: looped over energy bins

    } // endif: number of parameters was positive
				
	// ... else there are no parameters, so throw an exception
	else {
		std::string msg ="No model parameter for the computation of the "
						 "likelihood in observation \""+this->name()+
		                 "\" (ID "+this->id()+").\n";
		throw GException::invalid_value(G_POISSON_ONOFF,msg);
	}
		
    // Dump statistics
    #if G_OPT_DEBUG
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
    std::cout << "Number of bins excluded due to small model: " << n_small_model << std::endl;
    std::cout << "Number of bins with zero data: " << n_zero_data << std::endl;
    std::cout << "Sum of data (ON): " << sum_data << std::endl;
    std::cout << "Sum of model (ON): " << sum_model << std::endl;
    std::cout << "Statistics: " << value << std::endl;
    #endif
	
    // Optionally dump gradient and curvature matrix
    #if G_EVAL_DEBUG
    std::cout << *gradient << std::endl;
    std::cout << *curvature << std::endl;
    #endif
	
    // Timing measurement
    #if G_EVAL_TIMING
    #ifdef _OPENMP
    double t_elapse = omp_get_wtime()-t_start;
    #else
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    #endif
    std::cout << "GCTAOnOffObservation::optimizer::likelihood_poisson_onoff: CPU usage = "
	          << t_elapse << " sec" << std::endl;
    #endif
	
    // Return
    return value;
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
    m_instrument = "CTA";
    m_response   = NULL;
    m_ontime     = 0.0;
    m_livetime   = 0.0;
    m_deadc      = 1.0;
    m_on_spec.clear();
    m_off_spec.clear();
    m_arf.clear();
	m_bgd.clear();
    m_rmf.clear();
    m_on_regions.clear();
    m_off_regions.clear();
    m_alpha.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs CTA on-off observation.
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
	m_bgd         = obs.m_bgd;
    m_rmf         = obs.m_rmf;
    m_on_regions  = obs.m_on_regions;
    m_off_regions = obs.m_off_regions;
    m_alpha       = obs.m_alpha;

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
 * @brief Compute alpha parameter (ratio of On to Off effective areas)
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No CTA response found in CTA observation.
 ***************************************************************************/
void GCTAOnOffObservation::compute_alpha(const GCTAObservation& obs)
{
    // Get energy boundaries from On spectrum
	GEbounds ereco = m_on_spec.ebounds();
    
    // Continue only if there are spectral bins
    int nreco = ereco.size();
    if (nreco > 0) {
		
		// Initialise On/Off exposure ratios
        m_alpha.assign(nreco, 0.0);
    
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
		double       zenith  = obspnt.zenith();
		double       azimuth = obspnt.azimuth();

        // Loop over reconstructed energies 
        for (int i = 0; i < nreco; ++i) {
        
            // Get mean log10 energy in TeV of bin
            double logEreco = ereco.elogmean(i).log10TeV();
			
			// Initialise effective area totals
			double aon  = 0.0;
			double aoff = 0.0;
			
			// Loop over ON sky regions
            for (int m = 0; m < m_on_regions.size(); ++m)  {
				
				// If region is of type GSkyRegionCircle
				const GSkyRegionCircle* skyreg =
                      dynamic_cast<const GSkyRegionCircle*>(m_on_regions[m]);
				if (skyreg != NULL) {
					
					// Get centre of circular region
					GSkyDir centreg = skyreg->centre();
					
					// Compute position of region centre in instrument
                    // coordinates (should be replaced by proper method to
                    // yield theta and phi)
				    double theta = obsdir.dist(centreg);
				    double phi   = obsdir.posang(centreg);
					
					// Add up effective area
					aon += response->aeff(theta,
                                          phi,
                                          zenith,
                                          azimuth,
                                          logEreco);
					
				} // endif: sky region is of type GSkyRegionCircle
				
			} // Looped over ON regions
						
			// Loop over OFF sky regions
            for (int m = 0; m < m_off_regions.size(); ++m)  {
				
				// If region is of type GSkyRegionCircle
				const GSkyRegionCircle* skyreg =
                      dynamic_cast<const GSkyRegionCircle*>(m_off_regions[m]);
				if (skyreg != NULL) {
					
					// Get centre of circular region
					GSkyDir centreg = skyreg->centre();
					
					// Compute position of region centre in instrument
                    // coordinates (should be replaced by proper method to
                    // yield theta and phi)
				    double theta = obsdir.dist(centreg);
				    double phi   = obsdir.posang(centreg);
					
					// Add up effective area
					aoff += response->aeff(theta,
                                           phi,
                                           zenith,
                                           azimuth,
                                           logEreco);
					
                } // endif: sky region is of type GSkyRegionCircle
				
            } // endfor: looped over OFF regions
			
			// Compute alpha for this energy bin
			if (aoff > 0.0) {
				m_alpha[i] = aon/aoff;
			}
												 
        } // endfor: looped over reconstructed energies
	
    } // endif: there were energy bins

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
 * @todo Implement GCTAResponse::npred usage.
 ***************************************************************************/
void GCTAOnOffObservation::compute_arf(const GCTAObservation& obs)
{
    // Get energy boundaries from on spectrum
	GEbounds ereco = m_on_spec.ebounds();
    
    // Continue only if there are spectral bins
    int nreco = ereco.size();
    if (nreco > 0) {
    
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

        // Initialize ARF
        m_arf = GArf(ereco);

        // Loop over reconstructed energies
        for (int i = 0; i < nreco; ++i) {
        
            // Get mean energy of bin
            double logEreco = ereco.elogmean(i).log10TeV();
			
			// Initialize effective area for this bin
			m_arf[i] = 0.0;
			
			// Loop over On sky regions
            for (int m = 0; m < m_on_regions.size(); ++m)  {
				
				// If region is of type GSkyRegionCircle
				const GSkyRegionCircle* skyreg =
                      dynamic_cast<const GSkyRegionCircle*>(m_on_regions[m]);
				if (skyreg != NULL) {
					
					// Get centre of circular region
					GSkyDir centreg(skyreg->centre());
					
					// Compute position of region centre in instrument coordinates
				    double theta = obsdir.dist(centreg);
				    double phi   = obsdir.posang(centreg);
					
		            // Add up effective area
		            m_arf[i] += response->aeff(theta,
                                               phi,
                                               zenith,
                                               azimuth,
                                               logEreco);
									
				} // endif: sky region is of type GSkyRegionCircle
				
			} // endfor: looped over On regions
       
        } // endfor: looped over reconstructed energies
        
	} // endif: there were energy bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute background rate in Off observations
 *        (taken from IRF, does not include rescaling by spectral function)
 *
 * @param[in] obs CTA observation.
 * @param[in] models Model container.
 *
 * @exception GException::invalid_argument
 *            Specified observation does not contain relevant response or
 *            background information
 ***************************************************************************/
void GCTAOnOffObservation::compute_bgd(const GCTAObservation& obs,
                                       const GModels&         models)
{	
    // Get energy boundaries from on spectrum
	GEbounds ereco = m_on_spec.ebounds();
    
    // Continue only if there are spectral bins
    int nreco = ereco.size();
    if (nreco > 0) {
    		
		// Get CTA observation pointing direction, zenith, and azimuth
		GCTAPointing obspnt(obs.pointing());
		double       zenith  = obspnt.zenith();
		double       azimuth = obspnt.azimuth();
	
		// Initialize background rate array
    	m_bgd = GArf(ereco);
	
		// Loop over models to find background
    	for (int j = 0; j < models.size(); ++j) {
			
	    	// Get model pointer. Fall through if pointer is not valid
	    	const GModel* mptr = models[j];
	    	if (mptr == NULL) {
                continue;
            }
				
            // Fall through if model does not apply to specific instrument and
            // observation identifier
            if (!mptr->is_valid(this->instrument(), this->id())) {
                continue;
            }
				
            // Fall through if this model component is not a background
            // component of type GCTAModelIrfBackground
            if (dynamic_cast<const GCTAModelIrfBackground*>(mptr) == NULL) {
                continue;
            }
					
            // Get pointer on CTA IRF response
            const GCTAResponseIrf* rsp =
                  dynamic_cast<const GCTAResponseIrf*>(obs.response());
            if (rsp == NULL) {
                std::string msg = "Specified observation does not contain an "
                                  "IRF response.\n" + obs.print();
                throw GException::invalid_argument(G_COMPUTE_BGD, msg);
            }

            // Retrieve pointer to CTA background
            const GCTABackground* bgd = rsp->background();
            if (bgd == NULL) {
                std::string msg = "Specified observation contains no "
                                  "background information.\n" + obs.print();
                throw GException::invalid_argument(G_COMPUTE_BGD, msg);
            }
					
            // Loop over Off regions
            for (int m = 0; m < m_off_regions.size(); ++m)  {
						
                // Get region centre direction, fall through if region is not
                // of type GSkyRegionCircle
                const GSkyRegionCircle* skyreg =
                      dynamic_cast<const GSkyRegionCircle*>(m_off_regions[m]);
				if (skyreg == NULL) {
                    continue;
                }
								
                // Get centre of circular region
                GSkyDir centreg = skyreg->centre();
								
                // Set direction to centre of the region in instrument
                // coordinates
                GCTAInstDir offdir = obspnt.instdir(centreg);
							
                // Loop over energy bins
                for (int i = 0; i < nreco; ++i) {

                    // Get log energy bin mean
                    double logEreco = m_on_spec.ebounds().elogmean(i).log10TeV();

                    // Get rate in events/s/MeV
                    m_bgd[i] += (*bgd)(logEreco, offdir.detx(), offdir.dety()) *
                                skyreg->solidangle();

                } // endfor: looped over energy bins
				
            } // endfor: looped over sky regions
			
		} // endfor: looped over model component
		
	} // endif: there were spectral bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute RMF of On/Off observation
 *
 * @param[in] obs CTA observation.
 * @param[in] etrue True energy boundaries.
 *
 * @exception GException::invalid_value
 *            Observation does not contain IRF response
 *
 * @todo Check if we really should use reconstructed energy in the Aeff
 *       access.
 ***************************************************************************/
void GCTAOnOffObservation::compute_rmf(const GCTAObservation& obs,
                                       const GEbounds&        etrue)
{
    // Get reconstructed energy boundaries from on spectrum
	GEbounds ereco = m_on_spec.ebounds();
    
    // Continue only if there are spectral bins
    int ntrue = etrue.size();
    int nreco = ereco.size();
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
		
        // Initialize RMF
        m_rmf = GRmf(etrue, ereco);

        // Loop over reconstructed energy
        for (int ireco = 0; ireco < nreco; ++ireco) {

            // Get log10(E)
            double logEreco = ereco.elogmean(ireco).log10TeV();

            // Loop over true energy
            for (int itrue = 0; itrue < ntrue; ++itrue) {
                
                // Compute true energy and initialise weight for this bin
                double logEtrue = etrue.elogmean(itrue).log10TeV();
				double weight   = 0.0;
				
				// Loop over ON sky regions
	            for (int m = 0; m < m_on_regions.size(); ++m)  {
				
					// Fall through if region is not of type GSkyRegionCircle
					const GSkyRegionCircle* skyreg =
                          dynamic_cast<const GSkyRegionCircle*>(m_on_regions[m]);
					if (skyreg == NULL) {
                        continue;
                    }
					
                    // Get centre of circular region
                    GSkyDir centreg(skyreg->centre());
					
                    // Compute position of region centre in instrument coordinates
                    double theta = obsdir.dist(centreg);
                    double phi   = obsdir.posang(centreg);

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
									
				} // endfor: looped over ON regions
				
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
 * @brief Evaluate model for sky contribution and fill model gradients
 *
 * @param[in] models Model container.
 * @param[in,out] ibin Energy bin number.
 * @param[in,out] mod_grad Model gradient array.
 *
 * Computes the number of expected gamma-ray events in an On region for a
 * given energy bin, and returns it. Also computes the gradients of the
 * sky model spectral components. The method assumes that 
 * parameters are stored in the order spatial-spectral-temporal.
 *
 * @todo I think this method only works for point sources. What happens
 *       for an extended source?
 ***********************************************************************/
double GCTAOnOffObservation::model_on(const GModels& models,
                                      int            ibin,
                                      GVector*       mod_grad) const
{
	// Get total number of model parameters
	int npars = models.npars();
	
	// Initialize variables (vector has 0.0 values)
	double value = 0.0;
	int    ipar  = 0;
	
	// Continue only if bin number is in range and there are model parameters
	if ((ibin < m_on_spec.size()) && (npars > 0)) {
		
        // Get energy bin bounds
        const GEnergy emin   = m_on_spec.ebounds().emin(ibin);
        const GEnergy emax   = m_on_spec.ebounds().emax(ibin);
        const GEnergy emean  = m_on_spec.ebounds().elogmean(ibin);
        const double  ewidth = m_on_spec.ebounds().ewidth(ibin).MeV();

        // Compute normalisation factors
        double exposure  = m_on_spec.exposure();
        double norm_flux = m_arf[ibin] * exposure;
        double norm_grad = norm_flux   * ewidth;
				
        // Loop over models
        for (int j = 0; j < models.size(); ++j) {
			
            // Get model pointer. Fall through if pointer is not valid
            const GModel* mptr = models[j];
            if (mptr == NULL) {
                continue;
            }
				
            // Fall through if model does not apply to specific instrument
            // and observation identifier
            if (!mptr->is_valid(this->instrument(), this->id())) {
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
								
                // Determine the number of gamma-ray events in model by
                // computing the flux over the energy bin in ph/cm2/s
                // and multiplying this flux by the effective area (cm2)
                // and the livetime (s)
                value += spectral->flux(emin,emax) * norm_flux;
								
                // Determine the model gradients at the current energy
                spectral->eval_gradients(emean);

                // Loop over spectral model parameters
                for (int k = 0; k < spectral->size(); ++k, ++ipar)  {
                    GModelPar& par = (*spectral)[k];
                    if (par.is_free() && ipar < npars)  {
                        (*mod_grad)[ipar] = par.factor_gradient() * norm_grad;
                    }
                }
								
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
 * @brief Evaluate model for bgd contribution and fill model gradients
 *
 * @param[in] models Model container.
 * @param[in,out] ibin Energy bin number.
 * @param[in,out] mod_grad Model gradient array.
 * @return Number of expected background events.
 *
 * Returns the number of expected background events in the Off spectrum
 * for a given energy bin. The method computes also the gradients of the
 * background model spectral components.
 *
 * The method assumes that the model parameters are stored in the order
 * spatial-spectral-temporal.
 ***********************************************************************/
double GCTAOnOffObservation::model_off(const GModels& models,
							           int            ibin,
						               GVector*       mod_grad) const
{
	// Get total number of model parameters
	int npars = models.npars();
	
	// Initialize variables
	double value  = 0.0;
	int    ipar = 0;
	
	// Continue only if bin number is in range and if there are parameters
	if ((ibin < m_off_spec.size()) && (npars > 0))  {
				
        // Get energy bin mean and width
        const GEnergy emean  = m_off_spec.ebounds().elogmean(ibin);
        const double  ewidth = m_off_spec.ebounds().ewidth(ibin).MeV();

        // Compute normalisation factor
        double exposure = m_off_spec.exposure();
        double norm     = m_bgd[ibin] * exposure * ewidth;
			
        // Loop over models
        for (int j = 0; j < models.size(); ++j) {
				
            // Get model pointer. Fall through if pointer is not valid
            const GModel* mptr = models[j];
            if (mptr == NULL) {
                continue;
            }
					
            // Fall through if model does not apply to specific instrument
            // and observation identifier
            if (!mptr->is_valid(this->instrument(), this->id())) {
                ipar += mptr->size();
                continue;
            }
						
            // Fall through if model is not a IRF background component
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
                // computing the rate at the energy bin in events/s/MeV and
                // by multiplying the rate with the off livetime (s), the
                // energy bin width (MeV) and the spectral normalisation
                value += spectral->eval_gradients(emean) * norm;

                // Loop over spectral model parameters
                for (int k = 0; k < spectral->size(); ++k, ++ipar)  {
                    GModelPar& par=(*spectral)[k];
                    if (par.is_free() && ipar < npars)  {
                        (*mod_grad)[ipar] = par.factor_gradient() * norm;
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
