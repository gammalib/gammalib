/***************************************************************************
 *          GCTAOnOffObservation.cpp - CTA on-off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Chia-Chun Lu & Christoph Deil                    *
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
#include "GCTAOnOffObservation.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_WRITE                   "GCTAOnOffObservation::write(GXmlElement&)"
#define G_READ                     "GCTAOnOffObservation::read(GXmlElement&)"
#define G_FILL                 "GCTAOnOffObservation::fill(GCTAObservation&)"
#define G_COMPUTE_ARF   "GCTAOnOffObservation::compute_arf(GCTAObservation&)"
#define G_COMPUTE_RMF   "GCTAOnOffObservation::compute_rmf(GCTAObservation&,"\
                                                                " GEbounds&)"
#define G_POISSON_ONOFF  "GCTAOnOffObservation::likelihood_poisson_onoff(GModels&,"\
															            "GOptimizerPars&,"\
															            "GMatrixSparse&,"\
															            "GVector&,"\
																		"double&,"\
															            "double&)"
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
 * @brief Fill events in the ON and OFF spectra
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
        std::string msg =
            "No event list found in CTA observation \""+obs.name()+"\" "
            "(ID="+obs.id()+").\n"
            "ON/OFF observation can only be filled from event list.";
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
	
	// Set exposure times and compute alpha parameter 
	// (default just based on numbers of regions, eventually should include ratio of effective areas)
	m_ontime=obs.livetime();
	m_offtime=obs.livetime()*m_off_regions.size();
	m_alpha=m_ontime/m_offtime;
	
	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute response of ON/OFF observation
 *
 * @param[in] obs CTA observation.
 * @param[in] etrue True energy boundaries.
 ***************************************************************************/
void GCTAOnOffObservation::compute_response(const GCTAObservation& obs,
                                            const GEbounds&        etrue)
{
	// Compute response components
	compute_arf(obs);
	compute_rmf(obs, etrue);

	// Return
	return;
}

/***********************************************************************//**
 * @brief Compute ARF of ON/OFF observation
 *
 * @param[in] obs CTA observation.
 *
 * @todo Implement GCTAResponse::npred usage.
 ***************************************************************************/
void GCTAOnOffObservation::compute_arf(const GCTAObservation& obs)
{
    // Set constant response parameters
    const double theta   = 0.0;
    const double phi     = 0.0;
    const double zenith  = 0.0;
    const double azimuth = 0.0;

    // Get energy boundaries from on spectrum
	GEbounds ereco = m_on_spec.ebounds();
    
    // Continue only if there are spectral bins
    int nreco = ereco.size();
    if (nreco > 0) {
    
        // Get CTA response pointer. Throw an exception if no response is
        // found
        const GCTAResponse& response = obs.response();

        // Initialize ARF
        m_arf = GArf(ereco);

        // Loop over reconstructed energies
        for (int i = 0; i < nreco; ++i) {
        
            // Get mean energy of bin
            double logenergy = ereco.elogmean(i).log10TeV();

            // Set specresp value
            m_arf[i] = response.aeff(theta,
                                     phi,
                                     zenith,
                                     azimuth,
                                     logenergy);
        
        } // endfor: looped over reconstructed energies
        
	} // endif: there were energy bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute RMF of ON/OFF observation
 *
 * @param[in] obs CTA observation.
 * @param[in] etrue True energy boundaries.
 ***************************************************************************/
void GCTAOnOffObservation::compute_rmf(const GCTAObservation& obs,
                                       const GEbounds&        etrue)
{
    // Set constant response parameters
    const double theta   = 0.0;
    const double phi     = 0.0;
    const double zenith  = 0.0;
    const double azimuth = 0.0;
    
    // Get reconstructed energy boundaries from on spectrum
	GEbounds ereco = m_on_spec.ebounds();
    
    // Continue only if there are spectral bins
    int ntrue = etrue.size();
    int nreco = ereco.size();
    if (ntrue > 0 && nreco > 0) {
    
        // Get CTA response pointer
        const GCTAResponse& response = obs.response();

        // Initialize RMF
        m_rmf = GRmf(etrue, ereco);

        // Loop over reconstructed energy
        for (int ireco = 0; ireco < nreco; ++ireco) {

            // Compute reconstructed energy
            double eng_reco = ereco.elogmean(ireco).log10TeV();

            // Loop over true energy
            for (int itrue = 0; itrue < ntrue; ++itrue) {
                
                // Compute true energy
                double eng_true = etrue.elogmean(itrue).log10TeV();

                // Set RMF value
                m_rmf(itrue, ireco) = response.edisp(eng_reco,
                                                     theta,
                                                     phi,
                                                     zenith,
                                                     azimuth,
                                                     eng_true);
            } // endfor: looped over true energy
        } // endfor: looped over reconstructed energy
    } // endif: there were energy bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Print ON/OFF observation information
 *
 * @return String containing ON/OFF observation information.
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
 * @brief Evaluate model for sky contribution and fill model gradients
 *
 * @param[in] pars Model parameters.
 * @param[in,out] ibin Energy bin number.
 * @param[in,out] mod_grad Model gradient array.
 *
 * Computes the number of expected gamma events in an ON region for a
 * given energy bin, and returns it. Also computes the gradients of the
 * sky model spectral components. The method assumes that 
 * parameters are stored in the order spatial-spectral-temporal.
 *
 ***********************************************************************/
double GCTAOnOffObservation::model_on(const GModels&            models,
									  const GOptimizerPars&     pars,
										    int                 ibin,
										    GVector*            mod_grad)
{
	// Get number of parameters
	int npars = pars.size();
	
	// Initialize variables (vector has 0.0 values)
	double ngam=0.0;
	int i_par=0;
	if (mod_grad != NULL) delete mod_grad;
	mod_grad = new GVector(npars);
	
	// If bin number is in range
	if (ibin < m_on_spec.size())  {
		
	    // Check that there are parameters
	    if (npars > 0) {
		
		    // Create time object (empty, just needed in some calls)
			GTime time;
			
			// Get energy bin bounds
		    const GEnergy emin=m_on_spec.ebounds().emin(ibin);
		    const GEnergy emax=m_on_spec.ebounds().emax(ibin);
		    const GEnergy emean=m_on_spec.ebounds().elogmean(ibin);
				
		    // Loop over models 
		    for (int j = 0; j < models.size(); ++j) {
			
			    // Pointers to model components
			    GModelSpatial*  spaskyptr=NULL;
				GModelSpectral* speskyptr=NULL;
				GModelTemporal* temskyptr=NULL;
			
			    // Get model pointer. Continue only if pointer is valid
			    const GModel* mptr = models[j];
			    if (mptr != NULL) {
				
				    // Continue only if model applies to specific instrument and observation identifier
				    if (mptr->is_valid(this->instrument(), this->id())) {
					
					    // If this model component is a sky component
					    if (dynamic_cast<const GModelSky*>(mptr) != NULL) {
													
						    // Set pointer to sky
							const GModelSky* skyptr=dynamic_cast<const GModelSky*>(mptr);
							
							// Increase parameter counter for spatial parameter
							spaskyptr=skyptr->spatial();
							if (spaskyptr != NULL)  {
								i_par += spaskyptr->size();
							}
						
						    // Spectral component (the useful one)
						    speskyptr=skyptr->spectral();
						    
							// Number of gamma events in model
						    // (Get flux over energy bin in ph/cm2/s and multiply by effective area and time)
						    if (speskyptr != NULL)  {
							    ngam += speskyptr->flux(emin,emax)*m_arf[ibin]*m_ontime;
						    }
							
							// Gradients
							// Evaluate model at current energy (also fills gradients)
							double v=speskyptr->eval_gradients(emean,time);
							// Loop over spectral model parameters
							for (int k = 0; k < speskyptr->size(); ++k)  {  
								const GModelPar sppar=(*speskyptr)[k];
								if (sppar.is_free())  {
									mod_grad[i_par]=sppar.gradient();
									i_par++;
								}
								
							} // Looped over parameters of the spectral component of the current model
							
							// Increase parameter counter for temporal parameter
							temskyptr=skyptr->temporal();
							if (temskyptr != NULL)  {
								i_par += temskyptr->size();
							}
						
					    // ... else model is not of sky type, increase parameter counter and move on
					    } else {
						    i_par += mptr->size();
						    continue;
					    }
															
				    // ... else model does not apply to this obs, increase parameter counter
				    } else {
					    i_par += mptr->size();
				    }
				
			    } // Pointer to model component is not NULL
			
	        } // Looped over model components		
		
	    } // Model container is not empty	
		
	} // Bin number is in the range
	
	// Exit
	return ngam;

}


/***********************************************************************
 * @brief Evaluate model for bgd contribution and fill model gradients
 *
 * @param[in] pars Model parameters.
 * @param[in,out] ibin Energy bin number.
 * @param[in,out] mod_grad Model gradient array.
 *
 * Computes the number of expected background events in an OFF region for a
 * given energy bin, and returns it. Also computes the gradients of the
 * background model spectral components.  The method assumes that 
 * parameters are stored in the order spatial-spectral-temporal.
 *
 ***********************************************************************/
double GCTAOnOffObservation::model_off(const GModels&            models,
									   const GOptimizerPars&     pars,
									         int                 ibin,
									         GVector*            mod_grad)
{
	// Get number of parameters
	int npars = pars.size();
	
	// Initialize variables (vector has 0.0 values)
	double nbgd=0.0;
	int i_par=0;
	if (mod_grad != NULL) delete mod_grad;
	mod_grad = new GVector(npars);
	
	// If bin number is in range
	if (ibin < m_off_spec.size())  {
				
	    // Check that there are parameters
	    if (npars > 0) {
			
		    // Create time object (empty, just needed in some calls)
			GTime time;
			
			// Get energy bin bounds
		    const GEnergy emin=m_on_spec.ebounds().emin(ibin);
		    const GEnergy emax=m_on_spec.ebounds().emax(ibin);
		    const GEnergy emean=m_on_spec.ebounds().elogmean(ibin);
			
		    // Loop over models 
		    for (int j = 0; j < models.size(); ++j) {
				
			    // Pointers to model components
			    GModelSpatial*  spabgdptr=NULL;
				GModelSpectral* spebgdptr=NULL;
				GModelTemporal* tembgdptr=NULL;
				
			    // Get model pointer. Continue only if pointer is valid
			    const GModel* mptr = models[j];
			    if (mptr != NULL) {
					
				    // Continue only if model applies to specific instrument and observation identifier
				    if (mptr->is_valid(this->instrument(), this->id())) {
						
					    // If this model component is a background component
					    if (dynamic_cast<const GCTAModelBackground*>(mptr) != NULL) {
							
							// Set pointer to background
							const GCTAModelBackground* bgdptr=dynamic_cast<const GCTAModelBackground*>(mptr);
							
						    // Increase parameter counter for spatial parameter
							spabgdptr=bgdptr->spatial();
							if (spabgdptr != NULL)  {
								i_par += spabgdptr->size();
							}
							
						    // Spectral component (the useful one)
						    spebgdptr=bgdptr->spectral();
						    
							// Number of gamma events in model
						    // (Get flux over energy bin in ph/cm2/s and multiply by effective area and time)
						    if (spebgdptr != NULL)  {
							    nbgd += spebgdptr->flux(emin,emax)*m_offtime;
						    }
							
							// Gradients
							// Evaluate model at current energy (also fills gradients)
							double v=spebgdptr->eval_gradients(emean,time);
							// Loop over spectral model parameters
							for (int k = 0; k < spebgdptr->size(); ++k)  {  
								const GModelPar sppar=(*spebgdptr)[k];
								if (sppar.is_free())  {
									mod_grad[i_par]=sppar.gradient();
									i_par++;
								}
								
							} // Looped over parameters of the spectral component of the current model
							
							// Increase parameter counter for temporal parameter
							tembgdptr=bgdptr->temporal();
							if (tembgdptr != NULL)  {
								i_par += tembgdptr->size();
							}
							
							// ... else model is not of sky type, increase parameter counter and move on
					    } else {
						    i_par += mptr->size();
						    continue;
					    }
						
						// ... else model does not apply to this obs, increase parameter counter
				    } else {
					    i_par += mptr->size();
				    }
					
			    } // Pointer to model component is not NULL
				
	        } // Looped over model components		
			
	    } // Model container is not empty	
		
	} // Bin number is in the range
	
	// Exit
	return nbgd;
	
}


/***********************************************************************
 * @brief Evaluate log-likelihood function for ON-OFF analysis
 *
 * @param[in] obs Observation.
 * @param[in] pars Optimizer parameters.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] gradient Gradient.
 * @param[in,out] value Likelihood value.
 * @param[in,out] npred Number of predicted events.
 *
 * This method computes the log(likelihood) for one observation in the 
 * case of an ON-OFF analysis. The method loops over energy bins to 
 * update the function value and its first and second derivatives.
 * The number of ON and OFF counts Non et Noff are taken from the class 
 * members. The number of expected gamma and background events are computed 
 * from the spectral models of the relevant components in the model 
 * container (spatial and temporal models ignored so far). This is done
 * in methods model_on() and model_off().
 *
 * The curvature matrix includes only terms containing first derivatives.
 *
 ***********************************************************************/
double GCTAOnOffObservation::likelihood_poisson_onoff(const GModels&            models,
													  const GOptimizerPars&     pars,
											                GMatrixSparse*      curvature,
											                GVector*            gradient,
											                double&             value,
											                double&             npred)
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
    double init_value    = value;
	double init_npred    = npred;
    #endif
	
	// Get number of parameters
	int npars = pars.size();
	
	// Create model gradients array pointers (one for sky parameters, one for background parameters)
	GVector* sky_grad=NULL;
	GVector* bgd_grad=NULL;
	// Working arrays
	double* values = new double[npars];
	int*    inx    = new int[npars];
	
	// Create time object (empty, just needed in some calls)
	GTime time;
		
	// Check that there is at least one parameter
	if (npars > 0) {
		
		// Loop over all energy bins
	  for (int i = 0; i < m_on_spec.size(); ++i) {
				
	    // Get energy bin bounds
		const GEnergy emin=m_on_spec.ebounds().emin(i);
		const GEnergy emax=m_on_spec.ebounds().emax(i);
		const GEnergy emean=m_on_spec.ebounds().elogmean(i);
		  
		// Number of ON and OFF counts
		double non=m_on_spec[i];
		double noff=m_off_spec[i];
		double ngam=0.0;
		double nbgd=0.0;
		double nonpred=0.0;
			
		// Get number of gamma and background events (and corresponding spectral model gradients)
		ngam=model_on(models,pars,i,sky_grad);
		nbgd=model_off(models,pars,i,bgd_grad);
		  
		// Skip bin if model is too small (avoids -Inf or NaN gradients)
		nonpred= ngam+m_alpha*nbgd;
		if ((nbgd <= minmod) || (nonpred <= minmod)) {
		  #if G_OPT_DEBUG
		  n_small_model++;
		  #endif
		  continue;
		}
		  
		// Now we have all predicted gamma and background events for current energy bin
		// Update the log(likelihood) and predicted number of events
		value += -non*log(nonpred)+nonpred-noff*log(nbgd)+nbgd;
		npred += nonpred;
		  
		// Update statistics
		#if G_OPT_DEBUG
		n_used++;
		sum_data  += non;
		sum_model += nonpred;
		#endif
		
	    // Create index array of non-zero derivatives and initialise working array
		// (just needed in call to update hessian matrix)
		int ndev = 0;
		for (int j = 0; j < npars; ++j) {
		     if (((*sky_grad)[j] != 0.0  && !gammalib::is_infinite((*sky_grad)[j])) || 
				 ((*bgd_grad)[j] != 0.0  && !gammalib::is_infinite((*bgd_grad)[j])))  {
				  inx[ndev] = j;
				  ndev++;
			 }
		}
		  
		// Fill derivatives
		double fa=non/nonpred;
		double fb=fa/nonpred;
		double fc=m_alpha*fb;
		double fd=fc*m_alpha+noff/nbgd/nbgd;
		double sky_factor=1.0-fa;
		double bgd_factor=1.0+m_alpha-m_alpha*fa-noff/nbgd;
		// Loop over all parameters
		for (int j = 0; j < npars; ++j) {
			
			// If spectral model for sky component is non-zero and non-infinite
			if ((*sky_grad)[j] != 0.0  && !gammalib::is_infinite((*sky_grad)[j])) {
				
				// Gradient
				(*gradient)[j] += sky_factor*(*sky_grad)[j];
				
				// Hessian (from first-order derivatives only)
				for (int k = 0; k < npars; ++k) {
				
					// If spectral model for sky component is non-zero and non-infinite
					if ((*sky_grad)[k] != 0.0  && !gammalib::is_infinite((*sky_grad)[k])) {	
					    values[k]=(*sky_grad)[j]*(*sky_grad)[k]*fb;
					
					// If spectral model for sky component is non-zero and non-infinite
					} else if ((*bgd_grad)[k] != 0.0  && !gammalib::is_infinite((*bgd_grad)[k])) {
						values[k]=(*sky_grad)[j]*(*bgd_grad)[k]*fc;
						
					// ...else neither sky nor background
					} else {
					    values[k]=0.0;
					}
					
					// Update matrix
					curvature->add_to_column(j, values, inx, ndev);
					
				}
			
			// If spectral model for sky component is non-zero and non-infinite
			} else if ((*bgd_grad)[j] != 0.0  && !gammalib::is_infinite((*bgd_grad)[j])) {
				
				// Gradient
				(*gradient)[j] += bgd_factor*(*bgd_grad)[j];
				
				// Hessian (from first-order derivatives only)
				for (int k = 0; k < npars; ++k) {
					
					// If spectral model for sky component is non-zero and non-infinite
					if ((*sky_grad)[k] != 0.0  && !gammalib::is_infinite((*sky_grad)[k])) {	
					    values[k]=(*bgd_grad)[j]*(*sky_grad)[k]*fc;
						
						// If spectral model for sky component is non-zero and non-infinite
					} else if ((*bgd_grad)[k] != 0.0  && !gammalib::is_infinite((*bgd_grad)[k])) {
						values[k]=(*bgd_grad)[j]*(*bgd_grad)[k]*fd;
						
						// ...else neither sky nor background
					} else {
					    values[k]=0.0;
					}
					
					// Update matrix
					curvature->add_to_column(j, values, inx, ndev);
					
				}
			
			} 
						
		} // Looped over all parameters for derivatives computation
			
	  } // Looped over energy bins
		
	// ... else parameter array is not of type GModel so raise exception	
	} else {
		std::string msg ="No model parameter for the computation of the "
						 "likelihood in observation "+this->name()+
		                 " (ID "+this->id()+").\n";
		throw GException::invalid_value(G_POISSON_ONOFF,msg);
	}
	
    // Free working array
    if (values != NULL) delete [] values;
	if (inx    != NULL) delete [] inx;
	
    // Dump statistics
    #if G_OPT_DEBUG
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
    std::cout << "Number of bins excluded due to small model: " << n_small_model << std::endl;
    std::cout << "Number of bins with zero data: " << n_zero_data << std::endl;
    std::cout << "Sum of data (ON): " << sum_data << std::endl;
    std::cout << "Sum of model (ON): " << sum_model << std::endl;
    std::cout << "Initial statistics: " << init_value << std::endl;
    std::cout << "Delta statistics: " << value-init_value << std::endl;
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
    std::cout << "GCTAOnOffObservation::optimizer::poisson_onoff: CPU usage = "
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
    m_name.clear();
    m_id.clear();
    m_instrument = "CTAOnOff";
    m_on_spec.clear();
    m_off_spec.clear();
    m_arf.clear();
    m_rmf.clear();
    m_on_regions.clear();
    m_off_regions.clear();

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
    m_name        = obs.m_name;
    m_id          = obs.m_id;
    m_instrument  = obs.m_instrument;
    m_on_spec     = obs.m_on_spec;
    m_off_spec    = obs.m_off_spec;
    m_arf         = obs.m_arf;
    m_rmf         = obs.m_rmf;
    m_on_regions  = obs.m_on_regions;
    m_off_regions = obs.m_off_regions;

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
