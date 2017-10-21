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
#include "GIntegral.hpp"
#include "GMatrixSparse.hpp"
#include "GModels.hpp"
#include "GModelSky.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GSkyRegions.hpp"
#include "GSkyRegionMap.hpp"
#include "GOptimizerPars.hpp"
#include "GObservations.hpp"
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
#define G_CONSTRUCTOR1   "GCTAOnOffObservation::GCTAOnOffObservation(GPha&, "\
                                                       "GPha&, GArf&, GRmf&)"
#define G_CONSTRUCTOR2             "GCTAOnOffObservation(GObservations& obs)"
#define G_RESPONSE_SET           "GCTAOnOffObservation::response(GResponse&)"
#define G_RESPONSE_GET                     "GCTAOnOffObservation::response()"
#define G_WRITE                   "GCTAOnOffObservation::write(GXmlElement&)"
#define G_READ                     "GCTAOnOffObservation::read(GXmlElement&)"
#define G_LIKELIHOOD            "GCTAOnOffObservation::likelihood(GModels&, "\
                                          "GOptimizerPars&, GMatrixSparse&, "\
                                                "GVector&, double&, double&)"
#define G_LIKELIHOOD_CSTAT          "GCTAOnOffObservation::likelihood_cstat("\
                                "GModels&, GOptimizerPars&, GMatrixSparse&, "\
                                                "GVector&, double&, double&)"
#define G_LIKELIHOOD_WSTAT          "GCTAOnOffObservation::likelihood_wstat("\
                                "GModels&, GOptimizerPars&, GMatrixSparse&, "\
                                                "GVector&, double&, double&)"
#define G_SET        "GCTAOnOffObservation::set(GCTAObservation&, GSkyDir&, "\
                                            "GSkyRegionMap&, GSkyRegionMap&)"
#define G_COMPUTE_ARF  "GCTAOnOffObservation::compute_arf(GCTAObservation&, "\
                                                  "GSkyDir&, GSkyRegionMap&)"
#define G_COMPUTE_BGD  "GCTAOnOffObservation::compute_bgd(GCTAObservation&, "\
                                                            "GSkyRegionMap&)"
#define G_COMPUTE_ALPHA                "GCTAOnOffObservation::compute_alpha("\
                          "GCTAObservation&, GSkyRegionMap&, GSkyRegionMap&)"
#define G_COMPUTE_RMF  "GCTAOnOffObservation::compute_rmf(GCTAObservation&, "\
                                                            "GSkyRegionMap&)"

/* __ Constants __________________________________________________________ */
const double minmod = 1.0e-100;                      //!< Minimum model value
const double minerr = 1.0e-100;                //!< Minimum statistical error

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_LIKELIHOOD_DEBUG                //!< Debug likelihood computation
//#define G_N_GAMMA_DEBUG                   //!< Debug N_gamma computation


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
 * @param[in] pha_on On spectrum.
 * @param[in] pha_off Off spectrum.
 * @param[in] arf Auxiliary Response File.
 * @param[in] rmf Redistribution Matrix File.
 *
 * Constructs On/Off observation from On and Off spectra, an Auxiliary
 * Response File and a Redistribution Matrix File.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GPha& pha_on,
                                           const GPha& pha_off,
                                           const GArf& arf,
                                           const GRmf& rmf)
{
    // Initialise private
    init_members();

    // Set data members
    m_on_spec  = pha_on;
	m_off_spec = pha_off;
    m_arf      = arf;
    m_rmf      = rmf;

    // Set log true energy node array
    set_logetrue();

    // Check consistency of On/Off observation
    check_consistency(G_CONSTRUCTOR1);

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
 * @param[in] srcdir Point source location.
 *
 * Constructs On/Off observation by filling the On and Off spectra and
 * computing the Auxiliary Response File (ARF) and Redistribution Matrix
 * File (RMF). The method requires the specification of the true and
 * reconstructed energy boundaries, as well as the definition of On and Off
 * regions.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GCTAObservation& obs,
                                           const GSkyDir&         srcdir,
                                           const GEbounds&        etrue,
                                           const GEbounds&        ereco,
                                           const GSkyRegions&     on,
                                           const GSkyRegions&     off)
{
    // Initialise private
    init_members();

    // Initialise spectra
    m_on_spec  = GPha(ereco);
    m_off_spec = GPha(ereco);

    // Initialise response information
    m_arf = GArf(etrue);
    m_rmf = GRmf(etrue, ereco);

    // Set On/Off observation from CTA observation
    set(obs, srcdir, on, off);

    // Return
    return;
}


/***********************************************************************//**
 * @brief CTA On/Off observation stacking constructor
 *
 * @param[in] obs Observation container.
 *
 * GException::invalid_value
 *             Incompatible On/Off observation definition.
 *
 * Constructs On/Off observation by stacking all On/Off observation in the
 * observation container into a single On/Off observation.
 *
 * The constructor uses the following formulae:
 *
 * \f[
 *    N^{\rm on}(E_{\rm reco}) = \sum_i N^{\rm on}_i(E_{\rm reco})
 * \f]
 *
 * \f[
 *    N^{\rm off}(E_{\rm reco}) = \sum_i N^{\rm off}_i(E_{\rm reco})
 * \f]
 *
 * \f[
 *    \alpha(E_{\rm reco}) = \frac{\sum_i \alpha_i(E_{\rm reco})
 *                                        N^{\rm off}_i(E_{\rm reco})}
 *                                {\sum_i N^{\rm off}_i(E_{\rm reco})}
 * \f]
 *
 * \f[
 *    ARF(E_{\rm true}) = \frac{\sum_i ARF_i(E_{\rm true}) \tau_i}
 *                             {\sum_i \tau_i}
 * \f]
 *
 * \f[
 *    RMF(E_{\rm true}, E_{\rm reco}) =
 *    \frac{\sum_i RMF_i(E_{\rm true}, E_{\rm reco}) ARF_i(E_{\rm true}) \tau_i}
 *         {\sum_i ARF_i(E_{\rm true}) \tau_i}
 * \f]
 *
 * where
 * \f$N^{\rm on}_i(E_{\rm reco})\f$ is the On spectrum of observation \f$i\f$,
 * \f$N^{\rm off}_i(E_{\rm reco})\f$ is the Off spectrum of observation
 * \f$i\f$,
 * \f$\alpha_i(E_{\rm reco})\f$ is the background scaling of observation
 * \f$i\f$,
 * \f$ARF_i(E_{\rm true})\f$ is the Auxiliary Response File of observation
 * \f$i\f$,
 * \f$RMF_i(E_{\rm true}, E_{\rm reco})\f$ is the Redistribution Matrix File
 * of observation \f$i\f$, and
 * \f$\tau_i\f$ is the livetime of observation \f$i\f$.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GObservations& obs) :
                      GObservation()
{
    // Initialise private
    init_members();

    // Signal first On/Off observation
    bool first = true;

    // Initialise exposure
    double exposure = 0.0;

    // Loop over all observation in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get pointer to On/Off observation
        const GCTAOnOffObservation* onoff =
              dynamic_cast<const GCTAOnOffObservation*>(obs[i]);

        // Skip observation if it is not a On/Off observation
        if (onoff == NULL) {
            continue;
        }

        // Check consistency of On/Off observation
        onoff->check_consistency(G_CONSTRUCTOR2);

        // If this is the first On/Off observation then store the data to
        // initialise the data definition
        if (first) {

            // Store PHA, ARF and RMF
            m_on_spec  = onoff->on_spec();
            m_off_spec = onoff->off_spec();
            m_arf      = onoff->arf() * onoff->on_spec().exposure();
            m_rmf      = onoff->rmf();
            m_ontime   = onoff->ontime();
            m_livetime = onoff->livetime();
            exposure   = onoff->on_spec().exposure();

            // Compute RMF contribution
            for (int itrue = 0; itrue < m_rmf.ntrue(); ++itrue) {
                double arf = m_arf[itrue];
                for (int imeasured = 0; imeasured < m_rmf.nmeasured(); ++imeasured) {
                    m_rmf(itrue,imeasured) *= arf;
                }
            }

            // Signal that the On/Off definition has been set
            first = false;

        }

        // ... otherwise stack data
        else {

            // Check consistency of On spectrum
            if (m_on_spec.ebounds() != onoff->on_spec().ebounds()) {
                std::string msg = "Incompatible energy binning of On spectrum.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }

            // Check consistency of Off spectrum
            if (m_off_spec.ebounds() != onoff->off_spec().ebounds()) {
                std::string msg = "Incompatible energy binning of Off spectrum.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }

            // Check consistency of Arf
            if (m_arf.ebounds() != onoff->arf().ebounds()) {
                std::string msg = "Incompatible energy binning of ARF.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }

            // Check consistency of Rmf
            if (m_rmf.etrue() != onoff->rmf().etrue()) {
                std::string msg = "Incompatible true energy binning of RMF.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }
            if (m_rmf.emeasured() != onoff->rmf().emeasured()) {
                std::string msg = "Incompatible measured energy binning of RMF.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }

            // Compute background scaling factor
            for (int i = 0; i < m_on_spec.size(); ++i) {

                // Compute background scaling factor for spectral bin
                double n1 = m_off_spec[i];                // Counts from Off
                double n2 = onoff->off_spec()[i];         // Counts from Off
                double a1 = m_on_spec.backscal(i);        // Alpha from On
                double a2 = onoff->on_spec().backscal(i); // Alpha from On
                double n  = n1 + n2;
                double a  = (n > 0.0) ? (a1 * n1 + a2 * n2)/n : 1.0;

                // Set background scaling factor of On Spectrum
                m_on_spec.backscal(i,a);

            } // endfor: computed background scaling factors

            // Add On and Off spectrum
            m_on_spec  += onoff->on_spec();
            m_off_spec += onoff->off_spec();

            // Add ARF
            m_arf += onoff->arf() * onoff->on_spec().exposure();

            // Add RMF
            for (int itrue = 0; itrue < m_arf.size(); ++itrue) {
                double arf = onoff->arf()[itrue] * onoff->on_spec().exposure();
                for (int imeasured = 0; imeasured < m_rmf.nmeasured(); ++imeasured) {
                    m_rmf(itrue,imeasured) += onoff->rmf()(itrue,imeasured) * arf;
                }
            }

            // Add exposure
            exposure += onoff->on_spec().exposure();

            // Add ontime and livetime
            m_ontime   += onoff->ontime();
            m_livetime += onoff->livetime();

        } // endelse: stacked data

    } // endfor: looped over observations

    // Compute RMF
    for (int itrue = 0; itrue < m_arf.size(); ++itrue) {
        double arf = m_arf[itrue];
        for (int imeasured = 0; imeasured < m_rmf.nmeasured(); ++imeasured) {
            m_rmf(itrue,imeasured) /= arf;
        }
    }

    // Compute ARF
    if (exposure > 0.0) {
        m_arf /= exposure;
    }

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
 * @exception GException::invalid_value
 *            Invalid statistic attribute encountered
 *
 * Reads information for an On/Off observation from an XML element. The
 * expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="Pha_on"  file="..."/>
 *       <parameter name="Pha_off" file="..."/>
 *       <parameter name="Arf"     file="..."/>
 *       <parameter name="Rmf"     file="..."/>
 *     </observation>
 *
 * Optionally, the statistic used for maximum likelihood fitting can be
 * specified:
 *
 *     <observation name="..." id="..." instrument="..." statistic="...">
 *
 ***************************************************************************/
void GCTAOnOffObservation::read(const GXmlElement& xml)
{
    // clean object
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Read in user defined statistic for this observation
    if (xml.attribute("statistic") != "") {

        // Extract statistic value
        std::string statistic = gammalib::toupper(xml.attribute("statistic"));

        // If statistic is not POISSON, CSTAT or WSTAT than throw an exception
        if ((statistic != "POISSON") &&
            (statistic != "CSTAT")   &&
            (statistic != "WSTAT")) {
            std::string msg = "Invalid statistic \""+statistic+"\" encountered "
                              "in observation definition XML file for "
                              "\""+m_instrument+"\" observation with identifier "
                              "\""+xml.attribute("id")+"\". Only \"POISSON\" "
                              ", \"CSTAT\" or \"WSTAT\" are supported.";
            throw GException::invalid_value(G_READ, msg);
        }

        // Save statistic value
        this->statistic(xml.attribute("statistic"));

    }

    // Get file names
    std::string pha_on  = gammalib::xml_get_attr(G_READ, xml, "Pha_on",  "file");
    std::string pha_off = gammalib::xml_get_attr(G_READ, xml, "Pha_off", "file");
    std::string arf     = gammalib::xml_get_attr(G_READ, xml, "Arf",     "file");
    std::string rmf     = gammalib::xml_get_attr(G_READ, xml, "Rmf",     "file");

    // Expand file names
    pha_on  = gammalib::xml_file_expand(xml, pha_on);
    pha_off = gammalib::xml_file_expand(xml, pha_off);
    arf     = gammalib::xml_file_expand(xml, arf);
    rmf     = gammalib::xml_file_expand(xml, rmf);

    // Load files
    m_on_spec.load(pha_on);
	m_off_spec.load(pha_off);
    m_arf.load(arf);
    m_rmf.load(rmf);

    // Set log true energy node array
    set_logetrue();

    // Check consistency of On/Off observation
    check_consistency(G_READ);

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
 *     <observation name="..." id="..." instrument="..." statistic="...">
 *       <parameter name="Pha_on"  file="..."/>
 *       <parameter name="Pha_off" file="..."/>
 *       <parameter name="Arf"     file="..."/>
 *       <parameter name="Rmf"     file="..."/>
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

    // Set Arf parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Arf");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_arf.filename()));

    // Set Rmf parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Rmf");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_rmf.filename()));

    // Add user defined statistic attributes
    if (statistic() != "") {
        xml.attribute("statistic", statistic());
    }

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
 *            Invalid statistic encountered.
 ***************************************************************************/
double GCTAOnOffObservation::likelihood(const GModels& models,
                                        GVector*       gradient,
                                        GMatrixSparse* curvature,
                                        double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Extract statistic for this observation
    std::string statistic = gammalib::toupper(this->statistic());

    // Poisson statistic with modeled background
    if ((statistic == "POISSON") || (statistic == "CSTAT")) {
        value = likelihood_cstat(models, gradient, curvature, npred);
    }

    // ... or Poisson statistic with measured background
    else if (statistic == "WSTAT") {
        value = likelihood_wstat(models, gradient, curvature, npred);
    }

    // ... or unsupported
    else {
        std::string msg = "Invalid statistic \""+statistic+"\" encountered. "
                          "Either specify \"POISSON\", \"CSTAT\" or "
                          "\"WSTAT\".";
        throw GException::invalid_value(G_LIKELIHOOD, msg);
    }

    // Return likelihood
    return value;

}


/***********************************************************************//**
 * @brief Print On/Off observation information
 *
 * @param[in] chatter Chattiness.
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
    m_logetrue.clear();

    // Overwrite base class statistic
    m_statistic = "Poisson";

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
    m_instrument = obs.m_instrument;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;
    m_on_spec    = obs.m_on_spec;
    m_off_spec   = obs.m_off_spec;
    m_arf        = obs.m_arf;
    m_rmf        = obs.m_rmf;
    m_logetrue   = obs.m_logetrue;

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
 * @brief Set true energy node array
 ***************************************************************************/
void GCTAOnOffObservation::set_logetrue(void)
{
    // Clear node array
    m_logetrue.clear();

    // Continue only if there are true energies in Arf
    int netrue = m_arf.size();
    if (netrue > 0) {

        // Reserve space in node array
        m_logetrue.reserve(netrue);

        // Append all log mean energies to node array
        for (int i = 0; i < netrue; ++i) {

            // Get log mean of true energy in TeV
            double logE = m_arf.ebounds().elogmean(i).log10TeV();

            // Append energy to node array
            m_logetrue.append(logE);

        } // endfor: appended log mean energies

    } // endif: there were true energies in Arf

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check consistency of data members
 *
 * @param[in] method Calling method.
 *
 * Checks the consistency of data members and throw exceptions in case that
 * inconsistencies are found.
 ***************************************************************************/
void GCTAOnOffObservation::check_consistency(const std::string& method) const
{
    // Check that On spectrum is consistent with Off spectrum
    if (m_on_spec.ebounds() != m_off_spec.ebounds()) {
        std::string msg = "On and Off spectra are incompatible.";
        throw GException::invalid_value(method, msg);
    }

    // Check that On spectrum is consistent with RMF
    if (m_on_spec.ebounds() != m_rmf.emeasured()) {
        std::string msg = "Redistribution Matrix File is incompatible with "
                          "spectra.";
        throw GException::invalid_value(method, msg);
    }

    // Check that ARF is consistent with RMF
    if (m_arf.ebounds() != m_rmf.etrue()) {
        std::string msg = "Redistribution Matrix File is incompatible with "
                          "Auxiliary Response File.";
        throw GException::invalid_value(method, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set On/Off observation from a CTA observation
 *
 * @param[in] obs CTA observation.
 * @param[in] srcdir Point source location.
 * @param[in] on On regions.
 * @param[in] off Off regions.
 *
 * @exception GException::invalid_value
 *            No CTA event list found in CTA observation.
 *
 * Sets an On/Off observation from a CTA observation by filling the events
 * that fall in the On and Off regions into the PHA spectra and by computing
 * the corresponding ARF and RMF response functions.
 ***************************************************************************/
void GCTAOnOffObservation::set(const GCTAObservation& obs,
                               const GSkyDir&         srcdir,
                               const GSkyRegions&     on,
                               const GSkyRegions&     off)
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
		if (on.contains(dir)) {
			m_on_spec.fill(atom->energy());
		}
		if (off.contains(dir)) {
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

    // Convert all regions into region maps
    GSkyRegions reg_on;
    GSkyRegions reg_off;
    for (int i = 0; i < on.size(); ++i) {
        reg_on.append(GSkyRegionMap(on[i]));
    }
    for (int i = 0; i < off.size(); ++i) {
        reg_off.append(GSkyRegionMap(off[i]));
    }

	// Compute response components
	compute_arf(obs, srcdir, reg_on);
	compute_bgd(obs, reg_off);
	compute_alpha(obs, reg_on, reg_off);
	compute_rmf(obs, reg_on);

    // Set log true energy node array
    set_logetrue();

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute ARF of On/Off observation
 *
 * @param[in] obs CTA observation.
 * @param[in] srcdir Point source location.
 * @param[in] on On regions.
 *
 * @exception GException::invalid_value
 *            No CTA response found in CTA observation.
 ***************************************************************************/
void GCTAOnOffObservation::compute_arf(const GCTAObservation& obs,
                                       const GSkyDir&         srcdir,
                                       const GSkyRegions&     on)
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

        // Loop over true energies
        for (int i = 0; i < ntrue; ++i) {

            // Get mean energy of bin
            double logEtrue = etrue.elogmean(i).log10TeV();

            // Initialize effective area for this bin
            m_arf[i] = 0.0;

            // Initialize totals
            double totsolid = 0.0;
            double totpsf   = 0.0;

            // Loop over regions
            for (int k = 0; k < on.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* on_map = static_cast<const GSkyRegionMap*>(on[k]);

                // Loop over pixels in On region map and integrate effective
                // area
                for (int j = 0; j < on_map->nonzero_indices().size(); ++j) {

                    // Get pixel index
                    int pixidx = on_map->nonzero_indices()[j];

                    // Get direction to pixel center
                    GSkyDir pixdir = on_map->map().inx2dir(pixidx);

                    // Get solid angle subtended by this pixel
                    double pixsolid = on_map->map().solidangle(pixidx);

                    // Compute position of pixel centre in instrument coordinates
                    double theta = obsdir.dist(pixdir);
                    double phi   = obsdir.posang(pixdir);

                    // Add up effective area
                    m_arf[i] += response->aeff(theta,
                                               phi,
                                               zenith,
                                               azimuth,
                                               logEtrue) * pixsolid;

                    // Add pixel solid angle to total for averaging later
                    totsolid += pixsolid;

                    // Compute offset angle to source
                    double delta = srcdir.dist(pixdir);

                    // Integrate PSF
                    totpsf += response->psf(delta,
                                            theta,
                                            phi,
                                            zenith,
                                            azimuth,
                                            logEtrue) * pixsolid;

                } // endfor: looped over all pixels in region map

            } // endfor: looped over all regions

            // Average effective area over solid angle
            if (totsolid > 0.0) {
                m_arf[i] /= totsolid;
            }

            // Correct effective area by containment fraction
            if (totpsf >= 0.0 && totpsf <= 1.0) {
                m_arf[i] *= totpsf;
            }

        } // endfor: looped over true energies
        
	} // endif: there were energy bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute background rate in Off regions
 *
 * @param[in] obs CTA observation.
 * @param[in] off Off regions.
 *
 * @exception GException::invalid_argument
 *            Observation does not contain relevant response or background
 *            information
 *
 * Compute the background rate in units of events/s/MeV in the Off region
 * map and stores the result as additional column with name `BACKRESP` in
 * the Auxiliary Response File (ARF).
 ***************************************************************************/
void GCTAOnOffObservation::compute_bgd(const GCTAObservation& obs,
                                       const GSkyRegions&     off)
{
    // Get true energy boundaries from on Arf and interpret them as
    // reconstructed energies for the background rates
	GEbounds ereco = m_arf.ebounds();
    int      nreco = ereco.size();

    // Continue only if there are Arf bins
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

        // Loop over regions
        for (int k = 0; k < off.size(); ++k) {

            // Get pointer on sky region map
            const GSkyRegionMap* off_map = static_cast<const GSkyRegionMap*>(off[k]);

            // Loop over pixels in Off region map and integrate background
            // rate
            for (int j = 0; j < off_map->nonzero_indices().size(); ++j) {

                // Get pixel index
                int pixidx = off_map->nonzero_indices()[j];

                // Get direction to pixel center
                GSkyDir pixdir = off_map->map().inx2dir(pixidx);

                // Translate sky direction into instrument direction
                GCTAInstDir pixinstdir = obspnt.instdir(pixdir);

                // Get solid angle subtended by this pixel
                double pixsolid = off_map->map().solidangle(pixidx);

                // Loop over energy bins
                for (int i = 0; i < nreco; ++i) {

                    // Get log10(E/TeV) of mean reconstructed bin energy
                    double logEreco = ereco.elogmean(i).log10TeV();

                    // Get background rate in events/s/MeV
                    background[i] += (*bgd)(logEreco,
                                            pixinstdir.detx(),
                                            pixinstdir.dety()) * pixsolid;

                } // endfor: looped over energy bins

            } // endfor: looped over all pixels in map

        } // endfor: looped over all regions

        // Append background vector to ARF
        m_arf.append("BACKRESP", background);

    } // endif: there were spectral bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute vector of alpha parameters
 *
 * @param[in] obs CTA observation.
 * @param[in] on On regions.
 * @param[in] off Off regions.
 *
 * @exception GException::invalid_value
 *            No CTA response found in CTA observation.
 *
 * Compute the alpha parameters for all energy bins. The alpha parameter
 * gives the ratio between the On and Off region background acceptance
 * multiplied by the ratio between On and Off region solid angles.
 *
 ***************************************************************************/
void GCTAOnOffObservation::compute_alpha(const GCTAObservation& obs,
                                         const GSkyRegions&     on,
                                         const GSkyRegions&     off)
{
    // Get reconstructed energy boundaries from RMF
	GEbounds ereco = m_rmf.emeasured();
    int      nreco = ereco.size();

    // Continue only if there are reconstructed energy bins
    if (nreco > 0) {

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

        // Loop over reconstructed energies
        for (int i = 0; i < nreco; ++i) {

            // Get mean log10 energy in TeV of bin
            double logEreco = ereco.elogmean(i).log10TeV();

            // Initialise background rate totals
            double aon  = 0.0;
            double aoff = 0.0;

            // Loop over On regions
            for (int k = 0; k < on.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* on_map = static_cast<const GSkyRegionMap*>(on[k]);

                // Loop over pixels in On region map and integrate acceptance
                for (int j = 0; j < on_map->nonzero_indices().size(); ++j) {

                    // Get pixel index
                    int pixidx = on_map->nonzero_indices()[j];

                    // Get direction to pixel center
                    GSkyDir pixdir = on_map->map().inx2dir(pixidx);

                    // Translate sky direction into instrument direction
                    GCTAInstDir pixinstdir = obspnt.instdir(pixdir);

                    // Get solid angle subtended by this pixel
                    double pixsolid = on_map->map().solidangle(pixidx);

                    // Add up acceptance
                    aon += (*response->background())(logEreco,
                                                     pixinstdir.detx(),
                                                     pixinstdir.dety()) *
                                                     pixsolid;

                } // endfor: looped over all pixels in map

            } // endfor: looped over regions

            // Loop over Off regions
            for (int k = 0; k < off.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* off_map = static_cast<const GSkyRegionMap*>(off[k]);

                // Loop over pixels in Off region map and integrate acceptance
                for (int j = 0; j < off_map->nonzero_indices().size(); ++j) {

                    // Get pixel index
                    int pixidx = off_map->nonzero_indices()[j];

                    // Get direction to pixel center
                    GSkyDir pixdir = off_map->map().inx2dir(pixidx);

                    // Translate sky direction into instrument direction
                    GCTAInstDir pixinstdir = obspnt.instdir(pixdir);

                    // Get solid angle subtended by this pixel
                    double pixsolid = off_map->map().solidangle(pixidx);

                    // Add up acceptance
                    aoff += (*response->background())(logEreco,
                                                      pixinstdir.detx(),
                                                      pixinstdir.dety()) *
                                                      pixsolid;

                } // endfor: looped over all pixels in map

            } // endfor: looped over all regions

			// Compute alpha for this energy bin
            double alpha = (aoff > 0.0) ? aon/aoff : 1.0;

            // Set background scaling in On spectra
            m_on_spec.backscal(i, alpha);

        } // endfor: looped over reconstructed energies

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute RMF of On/Off observation
 *
 * @param[in] obs CTA observation.
 * @param[in] on On regions.
 *
 * @exception GException::invalid_value
 *            Observation does not contain IRF response
 *
 * Compute the energy redistribution matrix for an On/Off observation. The
 * method requires that the RMF energy axes have been defined before.
 ***************************************************************************/
void GCTAOnOffObservation::compute_rmf(const GCTAObservation& obs,
                                       const GSkyRegions&     on)
{
    // Get true and reconstructed energy boundaries from Rmf
    GEbounds etrue = m_rmf.etrue();
    GEbounds ereco = m_rmf.emeasured();
    int      ntrue = etrue.size();
    int      nreco = ereco.size();

    // Continue only if there are Rmf bins
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

        // Initialise Rmf matrix
        for (int itrue = 0; itrue < ntrue; ++itrue) {
            for (int ireco = 0; ireco < nreco; ++ireco) {
                m_rmf(itrue, ireco) = 0.0;
            }
        }

        // Initialise weight matrix
        GMatrixSparse weight(ntrue, nreco);

        // Loop over On regions
        for (int k = 0; k < on.size(); ++k) {

            // Get pointer on sky region map
            const GSkyRegionMap* on_map = static_cast<const GSkyRegionMap*>(on[k]);

            // Loop over pixels in On region map and integrate acceptance
            for (int j = 0; j < on_map->nonzero_indices().size(); ++j) {

                // Get pixel index
                int pixidx = on_map->nonzero_indices()[j];

                // Get direction to pixel center
                GSkyDir pixdir = on_map->map().inx2dir(pixidx);

                // Compute position of pixel centre in instrument coordinates
                double theta = obsdir.dist(pixdir);
                double phi   = obsdir.posang(pixdir);

                // Loop over true energy
                for (int itrue = 0; itrue < ntrue; ++itrue) {

                    // Compute log10 of true energy in TeV
                    double logEtrue = etrue.elogmean(itrue).log10TeV();

                    // Get effective area for weighting
                    double aeff = response->aeff(theta, phi,
                                                 zenith, azimuth,
                                                 logEtrue);

                    // Setup energy dispersion integral
                    GCTAOnOffObservation::edisp_kern integrand(response,
                                                               theta,
                                                               phi,
                                                               zenith,
                                                               azimuth,
                                                               logEtrue);
                    GIntegral integral(&integrand);
                    integral.eps(1.0e-4);

                    // Loop over reconstructed energy
                    for (int ireco = 0; ireco < nreco; ++ireco) {

                        // Get log of reconstructed energy boundaries
                        double ereco_min = std::log(ereco.emin(ireco).MeV());
                        double ereco_max = std::log(ereco.emax(ireco).MeV());

                        // Do Romberg integration
                        double value = integral.romberg(ereco_min, ereco_max);

                        // Update Rmf value and weight
                        m_rmf(itrue, ireco)  += value * aeff;
                        weight(itrue, ireco) += aeff;

                    } // endfor: looped over reconstructed energy

                } // endfor: looped over true energy

            } // endfor: looped over all pixels in map

        } // endfor: looped over all regions

        // Normalise Rmf matrix
        for (int itrue = 0; itrue < ntrue; ++itrue) {
            for (int ireco = 0; ireco < nreco; ++ireco) {
                if (weight(itrue, ireco) > 0.0) {
                    m_rmf(itrue, ireco) /= weight(itrue, ireco);
                }
            }
        }

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

                // Debug code
                #if defined(G_N_GAMMA_DEBUG)
                double rmf_sum = 0.0;
                #endif

                // Loop over true energy bins
                for (int itrue = 0; itrue < m_arf.size(); ++itrue) {

                    // Get RMF value. Continue only if it is positive
                    double rmf = m_rmf(itrue, ibin);
                    if (rmf <= 0.0) {
                        continue;
                    }

                    // Debug code
                    #if defined(G_N_GAMMA_DEBUG)
                    rmf_sum += rmf;
                    #endif

                    // Get true energy bin properties
                    GEnergy etruemin   = m_arf.ebounds().emin(itrue);
                    GEnergy etruemax   = m_arf.ebounds().emax(itrue);
                    GEnergy etruemean  = m_arf.ebounds().elogmean(itrue);
                    double  etruewidth = m_arf.ebounds().ewidth(itrue).MeV();

                    // Compute normalisation factors
                    double exposure  = m_on_spec.exposure();
                    double norm_flux = m_arf[itrue] * exposure * rmf;
                    double norm_grad = norm_flux * etruewidth;

                    // Determine number of gamma-ray events in model by
                    // computing the flux over the true energy bin in
                    // ph/cm2/s and multiplying this by effective area (cm2),
                    // livetime (s) and redistribution probability
                    value += spectral->flux(etruemin, etruemax) * norm_flux;

                    // Determine the model gradients at the current true
                    // energy. The eval() method needs a time in case that the
                    // spectral model has a time dependence. We simply use a
                    // dummy time here.
                    spectral->eval(etruemean, GTime(), true);

                    // Loop over spectral model parameters
                    for (int k = 0; k < spectral->size(); ++k) {
                        GModelPar& par = (*spectral)[k];
                        if (par.is_free() && (k+ipar) < npars)  {
                            (*grad)[k+ipar] += par.factor_gradient() * norm_grad;
                        }
                    } // endfor: looped over model parameters

                } // endfor: looped over true energy bins

                // Debug code
                #if defined(G_N_GAMMA_DEBUG)
                std::cout << "sum(Rmf) = " << rmf_sum << std::endl;
                #endif

                // Increment parameter counter for spectral parameter
                ipar += spectral->size();

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
    if ((ibin >= 0) && (ibin < m_on_spec.size()) && (npars > 0))  {

        // Initialise parameter index
        int ipar = 0;

        // Get reference to background response (events/MeV/s)
        const std::vector<double>& backresp = m_arf["BACKRESP"];

        // Get reconstructed energy bin mean and width
        GEnergy emean  = m_on_spec.ebounds().elogmean(ibin);
        double  ewidth = m_on_spec.ebounds().ewidth(ibin).MeV();

        // Perform log-log interpolation of background rate at reconstructed
        // energy
        m_logetrue.set_value(emean.log10TeV());
        double wgt_left   = m_logetrue.wgt_left();
        double wgt_right  = m_logetrue.wgt_right();
        double bkg_left   = backresp[m_logetrue.inx_left()];
        double bkg_right  = backresp[m_logetrue.inx_right()];
        double background = 0.0;
        if (bkg_left > 0.0 && bkg_right > 0.0) {
            background = std::exp(wgt_left  * std::log(bkg_left) +
                                  wgt_right * std::log(bkg_right));
        }
        if (background < 0.0) {
            background = 0.0;
        }

        // Compute normalisation factor (events)
        double exposure = m_on_spec.exposure();
        double norm     = background * exposure * ewidth;

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
                // computing the model normalization at the mean value of the
                // energy bin and multiplying the normalisation with the number
                // of background events. The eval() method needs a time in case
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


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for On/Off analysis in the
 * case of Poisson signal with modeled Poisson background
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
double GCTAOnOffObservation::likelihood_cstat(const GModels& models,
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

            // Get background scaling
            double alpha = m_on_spec.backscal(i);

            // Get number of gamma and background events (and corresponding
            // spectral model gradients)
            double ngam = N_gamma(models, i, &sky_grad);
            double nbgd = N_bgd(models, i, &bgd_grad);

            // Skip bin if model is too small (avoids -Inf or NaN gradients)
            double nonpred = ngam + alpha * nbgd;
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
            double fc         = alpha * fb;
            double fd         = fc * alpha + noff/(nbgd*nbgd);
            double sky_factor = 1.0 - fa;
            double bgd_factor = 1.0 + alpha - alpha * fa - noff/nbgd;

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
		throw GException::invalid_value(G_LIKELIHOOD_CSTAT, msg);
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
    std::cout << "Statistic: " << value << std::endl;
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
    std::cout << "GCTAOnOffObservation::optimizer::likelihood_cstat: CPU usage = "
	          << t_elapse << " sec" << std::endl;
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for On/Off analysis in the
 * case of Poisson signal with measured Poisson background
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
 * \f$N_{\gamma}\f$ is computed from the spectral models of the relevant
 * components in the model container (spatial and temporal components
 * are ignored so far). See the N_gamma() method for details about the
 * model computations. The number of background counts is derived based
 * on the measurement in the Off region by analytical minimization of
 * the Poisson likelihood, i.e., it is treated as a nuisance parameter.
 * See Appendix B of the XSpec user manual, section Poisson data with
 * Poisson background (cstat), for more details.
 ***********************************************************************/
double GCTAOnOffObservation::likelihood_wstat(const GModels& models,
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

    // Create model gradient vectors for sky parameters
    GVector sky_grad(npars);

    // Allocate working array
    GVector colvar(npars);

    // Check that there is at least one parameter
    if (npars > 0) {

        // Loop over all energy bins
        for (int i = 0; i < m_on_spec.size(); ++i) {

            // Reinitialize working arrays
            for (int j = 0; j < npars; ++j) {
                sky_grad[j] = 0.0;
            }

            // Get number of On and Off counts
            double non  = m_on_spec[i];
            double noff = m_off_spec[i];

            // Get background scaling
            double alpha = m_on_spec.backscal(i);

            // Get number of gamma and background events (and corresponding
            // spectral model gradients)
            double ngam = N_gamma(models, i, &sky_grad);

            // Precompute some values
            double alphap1  = alpha + 1.0;
            double alpharat = alpha / alphap1;

            // Calculate number of background events, profile likelihood value
            // and likelihood derivatives
            double nbgd;
            double nonpred;
            double dlogLdsky;
            double d2logLdsky2;

            // Special case noff = 0
            if (noff == 0.0) {
                if (ngam < non * alpharat) {
                    nbgd        = non / alphap1 - ngam / alpha;
                    nonpred     = ngam + alpha * nbgd;
                    value      += -ngam / alpha - non * std::log(alpharat);
                    dlogLdsky   = -1.0/alpha;
                    d2logLdsky2 = 0.0;
                }
                else {
                    nbgd        = 0.0;
                    nonpred     = ngam;
                    value      += ngam + non * (std::log(non) - std::log(ngam) - 1.0);
                    dlogLdsky   = 1.0 - non / ngam;
                    d2logLdsky2 = non / (ngam * ngam);
                }
            } // endif: noff = 0
            
            // Special case non = 0
            else if (non == 0.0) {
                nbgd        = noff / alphap1;
                nonpred     = ngam + alpha * nbgd;
                value      += ngam + noff * std::log(alphap1);
                dlogLdsky   = 1.0;
                d2logLdsky2 = 0.0;
            } // endif: non = 0

            // General case
            else {
                double alphat2  = 2.0 * alpha;
                double n1       = alpha * (non + noff) - alphap1 * ngam;
                double n2       = std::sqrt(n1 * n1 +
                                  4.0 * alpha * alphap1 * noff * ngam);
                nbgd            = (n1 + n2) / (alphat2 * alphap1);
                nonpred         = ngam + alpha * nbgd;
                value          += ngam + alphap1 * nbgd - non - noff;
                value          += -non * (std::log(nonpred) - std::log(non));
                value          += -noff * (std::log(nbgd) - std::log(noff));
                double dbgddgam = ((alphat2 * noff) / n2 - n1 / n2 -1.0) /
                                  alphat2;
                dlogLdsky       = 1.0 - non / nonpred +
                                  (1.0 - noff / nbgd) * dbgddgam;
                d2logLdsky2     = alphap1 / (alphat2 * n2);
                d2logLdsky2    -= (alphat2 * alphap1 * noff / n2 -
                                  alphap1 * n1 / n2) / (n2*n2);
            } // endelse: general case

            // Update Npred
            *npred += nonpred;

            // Update statistics
            #if defined(G_LIKELIHOOD_DEBUG)
            n_used++;
            sum_data  += non;
            sum_model += nonpred;
            #endif

            // Loop over all parameters
            for (int j = 0; j < npars; ++j) {

                // If spectral model for sky component is non-zero and
                // non-infinite then handle sky component gradients and
                // second derivatives including at least a sky component ...
                if (sky_grad[j] != 0.0  && !gammalib::is_infinite(sky_grad[j])) {

                    // Gradient
                    (*gradient)[j] +=  dlogLdsky * sky_grad[j];

                    // Hessian (from first-order derivatives only)
                    for (int k = 0; k < npars; ++k) {

                        // If spectral model for sky component is non-zero and
                        // non-infinite then we have the curvature element
                        // of a sky component
                        if (sky_grad[k] != 0.0  &&
                            !gammalib::is_infinite(sky_grad[k])) {
                            colvar[k] = sky_grad[j] * sky_grad[k] * d2logLdsky2;
                        }
                        
                        // ...... else if spectral model for background component
                        // or neither sky nor background
                        else {
                            colvar[k] = 0.0;
                        }

                    } // endfor: Hessian computation

                    // Update matrix
                    curvature->add_to_column(j, colvar);

                } // endif: spectral model is non-zero and non-infinite

            } // endfor: looped over all parameters for derivatives computation

        } // endfor: looped over energy bins

    } // endif: number of parameters was positive

    // ... else there are no parameters, so throw an exception
    else {
        std::string msg ="No model parameter for the computation of the "
                         "likelihood in observation \""+this->name()+
                         "\" (ID "+this->id()+").\n";
        throw GException::invalid_value(G_LIKELIHOOD_WSTAT, msg);
    }

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
    std::cout << "GCTAOnOffObservation::optimizer::likelihood_wstat: CPU usage = "
              << t_elapse << " sec" << std::endl;
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Energy dispersion integration kernel evaluation
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel for the energy dispersion.
 ***************************************************************************/
double GCTAOnOffObservation::edisp_kern::eval(const double& x)
{
    // Set energy
    GEnergy ereco;
    double expx = std::exp(x);
    ereco.MeV(expx);

    // Get function value
    double value = m_irf->edisp(ereco, m_theta, m_phi, m_zenith, m_azimuth,
                                m_logEtrue);

    // Correct for variable substitution
    value *= expx;

    // Return value
    return value;
}
