/***************************************************************************
 *          GCTAOnOffObservation.cpp - CTA On/Off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2019 by Chia-Chun Lu & Christoph Deil               *
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
#include "GModelData.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GSource.hpp"
#include "GSkyRegions.hpp"
#include "GSkyRegionMap.hpp"
#include "GSkyRegionCircle.hpp"
#include "GOptimizerPars.hpp"
#include "GObservations.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAAeff2D.hpp"
#include "GCTACubeBackground.hpp"
#include "GCTAOnOffObservation.hpp"
#include "GCTASupport.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Globals ____________________________________________________________ */
const GCTAOnOffObservation g_onoff_obs_cta_seed(true, "CTAOnOff");
const GCTAOnOffObservation g_onoff_obs_hess_seed(true, "HESSOnOff");
const GCTAOnOffObservation g_onoff_obs_magic_seed(true, "MAGICOnOff");
const GCTAOnOffObservation g_onoff_obs_veritas_seed(true, "VERITASOnOff");
const GObservationRegistry g_onoff_obs_cta_registry(&g_onoff_obs_cta_seed);
const GObservationRegistry g_onoff_obs_hess_registry(&g_onoff_obs_hess_seed);
const GObservationRegistry g_onoff_obs_magic_registry(&g_onoff_obs_magic_seed);
const GObservationRegistry g_onoff_obs_veritas_registry(&g_onoff_obs_veritas_seed);

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
#define G_SET  "GCTAOnOffObservation::set(GCTAObservation&, GModelSpatial&, "\
                                            "GSkyRegionMap&, GSkyRegionMap&)"
#define G_COMPUTE_ARF  "GCTAOnOffObservation::compute_arf(GCTAObservation&, "\
                                            "GModelSpatial&, GSkyRegionMap&)"
#define G_COMPUTE_ARF_CUT            "GCTAOnOffObservation::compute_arf_cut("\
                          "GCTAObservation&, GModelSpatial&, GSkyRegionMap&)"
#define G_COMPUTE_BGD  "GCTAOnOffObservation::compute_bgd(GCTAObservation&, "\
                                                            "GSkyRegionMap&)"
#define G_COMPUTE_ALPHA                "GCTAOnOffObservation::compute_alpha("\
                          "GCTAObservation&, GSkyRegionMap&, GSkyRegionMap&)"
#define G_COMPUTE_RMF  "GCTAOnOffObservation::compute_rmf(GCTAObservation&, "\
                                                            "GSkyRegionMap&)"
#define G_MODEL_BACKGROUND "GCTAOnOffObservation::model_background(GModels&)"
#define G_ARF_RADIUS_CUT              "GCTAOnOffObservation::arf_radius_cut("\
                                            "GCTAObservation&, GSkyRegions&)"

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
 * @brief Instrument constructor
 *
 * @param[in] dummy Dummy flag.
 * @param[in] instrument Instrument name.
 *
 * Constructs an empty CTA On/Off observation for a given instrument.
 *
 * This method is specifically used allocation of global class instances for
 * observation registry. By specifying explicit instrument names it is
 * possible to use the "CTA" module are for other Imaging Air Cherenkov
 * Telescopes. So far, the following instrument codes are supported:
 * "CTAOnOff", "HESSOnOff", "VERITASOnOff", "MAGICOnOff".
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const bool&        dummy,
                                           const std::string& instrument) :
                      GObservation()
{
    // Initialise members
    init_members();

    // Set instrument
    m_instrument = instrument;

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

    // Set the ontime, livetime and deadtime correction
    set_exposure();

    // Check consistency of On/Off observation
    check_consistency(G_CONSTRUCTOR1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief CTA observation constructor
 *
 * @param[in] obs CTA observation.
 * @param[in] models Models including source and background model.
 * @param[in] srcname Name of source in models.
 * @param[in] etrue True energy boundaries.
 * @param[in] ereco Reconstructed energy boundaries.
 * @param[in] on On regions.
 * @param[in] off Off regions.
 * @param[in] use_model_bkg Use model background.
 *
 * Constructs On/Off observation by filling the On and Off spectra and
 * computing the Auxiliary Response File (ARF) and Redistribution Matrix
 * File (RMF). The method requires the specification of the true and
 * reconstructed energy boundaries, as well as the definition of On and Off
 * regions.
 *
 * The @p use_model_bkg flag controls whether a background model should be
 * used for the computations or not. This impacts the computations of the
 * @c BACKSCAL column in the On spectrum and the @c BACKRESP column in
 * the Off spectrum. See the compute_alpha() and compute_bgd() methods for
 * more information on the applied formulae.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GCTAObservation& obs,
                                           const GModels&         models,
                                           const std::string&     srcname,
                                           const GEbounds&        etrue,
                                           const GEbounds&        ereco,
                                           const GSkyRegions&     on,
                                           const GSkyRegions&     off,
                                           const bool&            use_model_bkg)
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
    set(obs, models, srcname, on, off, use_model_bkg);

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
 *                                        B_i(E_{\rm reco}) \tau_i}
 *                                {\sum_i B_i(E_{\rm reco}) \tau_i}
 * \f]
 *
 * \f[
 *    B(E_{\rm reco}) = \frac{\sum_i B_i(E_{\rm reco}) \tau_i}{\sum_i \tau_i}
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
 * \f$B_i(E_{\rm reco})\f$ is the background rate for observation \f$i\f$,
 * \f$ARF_i(E_{\rm true})\f$ is the Auxiliary Response File of observation
 * \f$i\f$,
 * \f$RMF_i(E_{\rm true}, E_{\rm reco})\f$ is the Redistribution Matrix File
 * of observation \f$i\f$, and
 * \f$\tau_i\f$ is the livetime of observation \f$i\f$.
 *
 * The method extracts the instrument name from the first On/Off observation
 * in the observation container and only stacks subsequent On/Off
 * observations with the same instrument name.
 ***************************************************************************/
GCTAOnOffObservation::GCTAOnOffObservation(const GObservations& obs) :
                      GObservation()
{
    // Initialise private
    init_members();

    // Signal first On/Off observation
    bool first = true;

    // Initialise total exposure
    double total_exposure = 0.0;

    // Allocate observation energy range
    GEnergy emin_obs;
    GEnergy emax_obs;

    // Allocate instrument
    std::string instrument;

    // Loop over all observation in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get pointer to On/Off observation
        const GCTAOnOffObservation* onoff =
              dynamic_cast<const GCTAOnOffObservation*>(obs[i]);

        // Skip observation if it is not a On/Off observation
        if (onoff == NULL) {
            continue;
        }

        // If the instrument name is empty then set the instrument now.
        // Otherwise, skip the observation if it does not correspond
        // to the same instrument.
        if (instrument.empty()) {
            instrument = onoff->instrument();
        }
        else if (instrument != onoff->instrument()) {
            continue;
        }

        // Check consistency of On/Off observation
        onoff->check_consistency(G_CONSTRUCTOR2);

        // Get energy boundaries of observation
        GEnergy emin = onoff->on_spec().emin_obs();
        GEnergy emax = onoff->on_spec().emax_obs();

        // Get stacked ARF and RMF
        GArf arf_stacked = onoff->arf();
        GRmf rmf_stacked = onoff->rmf();

        // Get exposure for observation
        double exposure = onoff->on_spec().exposure();

        // If this is the first On/Off observation then store the data to
        // initialise the data definition
        if (first) {

            // Store PHA, ARF and RMF, set exposure, ontime and livetime
            m_on_spec      = onoff->on_spec();
            m_off_spec     = onoff->off_spec();
            m_arf          = arf_stacked * exposure;
            m_rmf          = rmf_stacked;
            total_exposure = exposure;
            m_ontime       = onoff->ontime();
            m_livetime     = onoff->livetime();

            // Compute number of background events/MeV from BACKRESP column
            // and store result intermediately into BACKRESP column of off
            // spectum
            std::vector<double>& backresp = m_off_spec["BACKRESP"];
            for (int k = 0; k < backresp.size(); ++k) {
                backresp[k] *= exposure;
            }

            // Compute background scaling factor contribution
            for (int k = 0; k < m_on_spec.size(); ++k) {
                double  background = onoff->off_spec()["BACKRESP"][k];
                double  alpha      = onoff->on_spec().backscal(k);
                double  scale      = alpha * background * exposure;
                m_on_spec.backscal(k,scale);
            }

            // Compute RMF contribution
            for (int itrue = 0; itrue < m_rmf.ntrue(); ++itrue) {
                double arf = m_arf[itrue];
                for (int ireco = 0; ireco < m_rmf.nmeasured(); ++ireco) {
                    m_rmf(itrue,ireco) *= arf;
                }
            }

            // Set observation energy range
            emin_obs = emin;
            emax_obs = emax;

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

            // Check consistency of ARF
            if (m_arf.ebounds() != onoff->arf().ebounds()) {
                std::string msg = "Incompatible energy binning of ARF.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }

            // Check consistency of RMF
            if (m_rmf.etrue() != onoff->rmf().etrue()) {
                std::string msg = "Incompatible true energy binning of RMF.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }
            if (m_rmf.emeasured() != onoff->rmf().emeasured()) {
                std::string msg = "Incompatible measured energy binning of RMF.";
                throw GException::invalid_value(G_CONSTRUCTOR2, msg);
            }

            // Add On and Off spectrum
            m_on_spec  += onoff->on_spec();
            m_off_spec += onoff->off_spec();

            // Compute background scaling factor contribution
            for (int k = 0; k < m_on_spec.size(); ++k) {
                double  background = onoff->off_spec()["BACKRESP"][k];
                double  alpha      = onoff->on_spec().backscal(k);
                double  scale      = alpha * background * exposure;
                m_on_spec.backscal(k, m_on_spec.backscal(k) + scale);
            }

            // Add ARF
            m_arf += arf_stacked * exposure;

            // Add number of background events/MeV from BACKRESP column
            const std::vector<double>& src = onoff->off_spec()["BACKRESP"];
            std::vector<double>&       dst = m_off_spec["BACKRESP"];
            for (int k = 0; k < dst.size(); ++k) {
                dst[k] += src[k] * exposure;
            }

            // Add RMF
            for (int itrue = 0; itrue < m_rmf.ntrue(); ++itrue) {
                double arf = arf_stacked[itrue] * exposure;
                for (int ireco = 0; ireco < m_rmf.nmeasured(); ++ireco) {
                    m_rmf(itrue,ireco) += rmf_stacked(itrue,ireco) * arf;
                }
            }

            // Add exposure, ontime and livetime
            total_exposure += exposure;
            m_ontime       += onoff->ontime();
            m_livetime     += onoff->livetime();

            // Update observation energy range
            if (emin < emin_obs) {
                emin_obs = emin;
            }
            if (emax > emax_obs) {
                emax_obs = emax;
            }

        } // endelse: stacked data

    } // endfor: looped over observations

    // Compute stacked background scaling factor
    for (int k = 0; k < m_on_spec.size(); ++k) {
        double norm = m_off_spec["BACKRESP"][k];
        if (norm > 0.0) {
            double scale = m_on_spec.backscal(k) / norm;
            m_on_spec.backscal(k,scale);
        }
        else {
            m_on_spec.backscal(k,0.0);
        }
    }

    // Compute RMF
    for (int itrue = 0; itrue < m_rmf.ntrue(); ++itrue) {
        double arf = m_arf[itrue];
        if (arf > 0.0) {
            for (int ireco = 0; ireco < m_rmf.nmeasured(); ++ireco) {
                m_rmf(itrue,ireco) /= arf;
            }
        }
    }

    // Compute ARF
    if (total_exposure > 0.0) {
        m_arf /= total_exposure;
    }

    // Compute background events/MeV/sec and store them in BACKRESP column
    if (total_exposure > 0.0) {
        std::vector<double>& backresp = m_off_spec["BACKRESP"];
        for (int k = 0; k < backresp.size(); ++k) {
            backresp[k] /= total_exposure;
        }
    }

    // Set energy boundaries of observation
    m_on_spec.emin_obs(emin_obs);
    m_on_spec.emax_obs(emax_obs);
    m_off_spec.emin_obs(emin_obs);
    m_off_spec.emax_obs(emax_obs);

    // Set instrument name
    m_instrument = instrument;

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
    // Retrieve CTA response pointer
    const GCTAResponse* ptr = gammalib::cta_rsp(G_RESPONSE_SET, rsp);

    // Free existing response only if it differs from current response. This
    // prevents unintential deallocation of the argument
    if ((m_response != NULL) && (m_response != &rsp)) {
        delete m_response;
    }

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

    // Set the ontime, livetime and deadtime correction
    set_exposure();

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


/***********************************************************************
 * @brief Compute predicted gamma-ray counts for given model
 *
 * @param[in] models Model container.
 * @return Model PHA with predicted gamma-ray counts.
 *
 * Returns spectrum of predicted number of source events in the On
 * regions.
 ***********************************************************************/
GPha GCTAOnOffObservation::model_gamma(const GModels& models) const
{
    // Get number of model parameters in model container
    int npars = models.npars();

    // Create dummy gradient vectors to provide required argument to N_gamma
    GVector dummy_grad(npars);

    // Initialise output spectrum
    GEbounds ereco  = m_on_spec.ebounds();
    GPha     gammas = GPha(ereco);

    // Loop over all energy bins
    for (int i = 0; i < m_on_spec.size(); ++i) {
        gammas.fill(ereco.emean(i), N_gamma(models, i, &dummy_grad));
    }

    // Return predicted gamma-ray counts spectrum
    return gammas;
}


/***********************************************************************
 * @brief Compute predicted background counts PHA for given model
 *
 * @param[in] models Model container.
 * @return Model background PHA.
 *
 * Returns spectrum of predicted number of background events in the Off
 * regions. The computation method changed depending on the statistic used
 * to match the spectrum used for likelihood fitting. 
 * To be noted: for the WSTAT statistic the model background depends on
 * on the model adopted for gamma rays, that therefore need to be passed
 * to this function.
 ***********************************************************************/
GPha GCTAOnOffObservation::model_background(const GModels& models) const
{
    // Get number of model parameters in model container
    int npars = models.npars();

    // Create dummy gradient vectors to provide required argument to N_gamma
    GVector dummy_grad(npars);

    // Initialise output spectrum
    GEbounds ereco = m_on_spec.ebounds();
    GPha     bgds  = GPha(ereco);

    // Extract statistic for this observation
    std::string statistic = gammalib::toupper(this->statistic());

    // Loop over all energy bins
    for (int i = 0; i < m_on_spec.size(); ++i) {

        // Initialise variable to store number of background counts
        double nbgd = 0.0;

        // Choose background evaluation method based on fit statistics

        // CSTAT
        if ((statistic == "POISSON") || (statistic == "CSTAT")) {
            nbgd = N_bgd(models, i, &dummy_grad);
        }

        // ... or Poisson statistic with measured background
        else if (statistic == "WSTAT") {

            // Get number of On and Off counts
            double non  = m_on_spec[i];
            double noff = m_off_spec[i];

            // Get background scaling
            double alpha = m_on_spec.backscal(i);

            // Get number of gamma and background events (and corresponding
            // spectral model gradients)
            double ngam = N_gamma(models, i,  &dummy_grad);

            // Initialise variables for likelihood computation
            double nonpred     = 0.0;
            double dlogLdsky   = 0.0;
            double d2logLdsky2 = 0.0;

            // Perform likelihood profiling and derive likelihood value
            // The number of background counts is updated to the profile value
            double logL = wstat_value(non, noff, alpha,  ngam, nonpred, nbgd,
                                      dlogLdsky,d2logLdsky2);

        } // endelse: WSTAT

        // ... or unsupported
        else {
            std::string msg = "Invalid statistic \""+statistic+"\" encountered. "
                              "Either specify \"POISSON\", \"CSTAT\" or "
                              "\"WSTAT\".";
            throw GException::invalid_value(G_MODEL_BACKGROUND, msg);
        }

        // Fill background spectrum
        bgds.fill(ereco.emean(i), nbgd);

    } // endfor: looped over energy bins

    // Return model count spectrum
    return bgds;
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
        result.append("\n"+gammalib::parformat("Instrument")+instrument());
        result.append("\n"+gammalib::parformat("Statistic")+statistic());
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(ontime())+" s");
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(livetime())+" s");
        result.append("\n"+gammalib::parformat("Deadtime correction"));
        result.append(gammalib::str(m_deadc));

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
 * @brief Set ontime, livetime and deadtime correction factor
 ***************************************************************************/
void GCTAOnOffObservation::set_exposure(void)
{
    // Set the ontime, livetime and deadtime correction
    m_ontime   = m_on_spec.exposure();
    m_livetime = m_on_spec.exposure();
    m_deadc    = 1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check consistency of data members
 *
 * @param[in] method Calling method.
 *
 * @exception GException::invalid_value
 *            Inconsistency found.
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
 * @brief Return ARF for stacking
 *
 * @param[in] arf Auxiliary Response File.
 * @param[in] emin Minimum observation energy.
 * @param[in] emax Maximum observation energy.
 * @return Auxiliary Response File for stacking.
 *
 * Returns an Auxiliary Response File for stacking that has all elements with
 * true energies outside the interval [@p emin, @p emax] clipped to zero.
 * This prevents spill over from one observation into another.
 ***************************************************************************/
GArf GCTAOnOffObservation::arf_stacked(const GArf&    arf,
                                       const GEnergy& emin,
                                       const GEnergy& emax) const
{
    // Copy ARF
    GArf arf_stacked = arf;

    // Get true energy boundaries of ARF
    const GEbounds& etrue = arf_stacked.ebounds();

    // Set all ARF bins outside energy interval to zero
    for (int itrue = 0; itrue < etrue.size(); ++itrue) {
        if ((etrue.emax(itrue) < emin) || (etrue.emin(itrue) > emax)) {
            arf_stacked[itrue] = 0.0;
        }
    }

    // Return ARF
    return arf_stacked;
}


/***********************************************************************//**
 * @brief Return RMF for stacking
 *
 * @param[in] rmf Redistribution Matrix File.
 * @param[in] emin Minimum observation energy.
 * @param[in] emax Maximum observation energy.
 * @return Redistribution Matrix File for stacking.
 *
 * Returns a Redistribution Matrix File for stacking that has all matrix
 * elements with true energies outside the interval [@p emin, @p emax]
 * clipped to zero. This prevents spill over from one observation into
 * another.
 *
 * To correct for the missing events at the edges, the Redistribution Matrix
 * File is renormalized so that for each reconstructed energy the sum over
 * the true energies is the same as before the clipping.
 ***************************************************************************/
GRmf GCTAOnOffObservation::rmf_stacked(const GRmf&    rmf,
                                       const GEnergy& emin,
                                       const GEnergy& emax) const
{
    // Copy RMF
    GRmf rmf_stacked = rmf;

    // Get true energy boundaries of RMF
    const GEbounds& etrue = rmf_stacked.etrue();
    int             nreco = rmf_stacked.emeasured().size();

    // Get RMF totals for all reconstructed energy bins
    std::vector<double> sums(nreco, 0.0);
    for (int ireco = 0; ireco < nreco; ++ireco) {
        for (int itrue = 0; itrue < etrue.size(); ++itrue) {
            sums[ireco] += rmf_stacked(itrue, ireco);
        }
    }

    // Set all RMF bins outside energy interval to zero (clipping)
    for (int itrue = 0; itrue < etrue.size(); ++itrue) {
        if ((etrue.emax(itrue) < emin) || (etrue.emin(itrue) > emax)) {
            for (int ireco = 0; ireco < nreco; ++ireco) {
                rmf_stacked(itrue, ireco) = 0.0;
            }
        }
    }

    // Renormalize RMF
    for (int ireco = 0; ireco < nreco; ++ireco) {
        double sum = 0.0;
        for (int itrue = 0; itrue < etrue.size(); ++itrue) {
            sum += rmf_stacked(itrue, ireco);
        }
        if (sum > 0.0) {
            double norm = sums[ireco] / sum;
            for (int itrue = 0; itrue < etrue.size(); ++itrue) {
                rmf_stacked(itrue, ireco) *= norm;
            }
        }
    }

    // Return RMF
    return rmf_stacked;
}


/***********************************************************************//**
 * @brief Set On/Off observation from a CTA observation
 *
 * @param[in] obs CTA observation.
 * @param[in] models Models including source and background model.
 * @param[in] srcname Name of soucre in models.
 * @param[in] on On regions.
 * @param[in] off Off regions.
 * @param[in] use_model_bkg Use model background.
 *
 * @exception GException::invalid_value
 *            No CTA event list found in CTA observation.
 *
 * Sets an On/Off observation from a CTA observation by filling the events
 * that fall in the On and Off regions into the PHA spectra and by computing
 * the corresponding ARF and RMF response functions.
 *
 * The @p use_model_bkg flags controls whether the background model should
 * be used for computations or not. This impacts the computations of the
 * `BACKSCAL` column in the On spectrum and the `BACKRESP` column in the Off
 * spectrum. See the compute_alpha() and compute_bgd() methods for more
 * information on the applied formulae.
 ***************************************************************************/
void GCTAOnOffObservation::set(const GCTAObservation& obs,
                               const GModels&         models,
                               const std::string&     srcname,
                               const GSkyRegions&     on,
                               const GSkyRegions&     off,
                               const bool&            use_model_bkg)
{
    // Get CTA event list pointer
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
    if (events == NULL) {
        std::string msg = "No event list found in CTA observation \""+
                          obs.name()+"\" (ID="+obs.id()+"). ON/OFF observation "
                          "can only be filled from event list.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Get spatial component of source model
    const GModelSky* model = dynamic_cast<const GModelSky*>(models[srcname]);
    if (model == NULL) {
        std::string msg = "Model \""+srcname+"\" is not a sky model. Please "
                          "specify the name of a sky model.";
        throw GException::invalid_value(G_SET, msg);
    }
    const GModelSpatial& spatial = *(model->spatial());

    // If background model is needed then extract background model components
    // from model container
    GModels bkg_models;
    if (use_model_bkg) {
        for (int i = 0; i < models.size(); ++i) {
            if ((dynamic_cast<const GModelData*>(models[i]) != NULL) &&
                (models[i]->classname().substr(0,4) == "GCTA")) {
                bkg_models.append(*models[i]);
            }
        }
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

    // Get effective area radius cut
    double rad_max = arf_rad_max(obs, on);

	// Compute response components
    if (rad_max > 0.0) {
        compute_arf_cut(obs, spatial, reg_on);
    }
    else {
        compute_arf(obs, spatial, reg_on);
    }
	compute_bgd(obs, reg_off, bkg_models, use_model_bkg);
	compute_alpha(obs, reg_on, reg_off, bkg_models, use_model_bkg);
	compute_rmf(obs, reg_on);

    // Apply reconstructed energy boundaries
    apply_ebounds(obs);

    // Set observation energy band
    m_on_spec.emin_obs(obs.ebounds().emin());
    m_on_spec.emax_obs(obs.ebounds().emax());
    m_off_spec.emin_obs(obs.ebounds().emin());
    m_off_spec.emax_obs(obs.ebounds().emax());

    // Get mission, instrument and response name
    std::string mission;
    std::string instrument;
    std::string rspname;
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(obs.response());
    if (rsp != NULL) {
        mission    = gammalib::toupper(rsp->caldb().mission());
        instrument = gammalib::toupper(rsp->caldb().instrument());
        rspname    = gammalib::toupper(rsp->rspname());
    }

    // Set header keywords
    GFitsHeader arf_header;
    GFitsHeader pha_header;
    GFitsHeader rmf_header;
    if (!mission.empty()) {
        arf_header.append(GFitsHeaderCard("TELESCOP", mission, "Telescope"));
        pha_header.append(GFitsHeaderCard("TELESCOP", mission, "Telescope"));
        rmf_header.append(GFitsHeaderCard("TELESCOP", mission, "Telescope"));
    }
    if (!instrument.empty()) {
        arf_header.append(GFitsHeaderCard("INSTRUME", instrument, "Instrument"));
        pha_header.append(GFitsHeaderCard("INSTRUME", instrument, "Instrument"));
        rmf_header.append(GFitsHeaderCard("INSTRUME", instrument, "Instrument"));
    }
    if (!rspname.empty()) {
        arf_header.append(GFitsHeaderCard("RESPNAME", rspname, "Response name"));
        pha_header.append(GFitsHeaderCard("RESPNAME", rspname, "Response name"));
        rmf_header.append(GFitsHeaderCard("RESPNAME", rspname, "Response name"));
    }
    if (rad_max > 0.0) {
        arf_header.append(GFitsHeaderCard("RAD_MAX", rad_max, "[deg] Applied radius cut"));
    }

    // Store header keywords
    m_on_spec.header(pha_header);
    m_off_spec.header(pha_header);
    m_arf.header(arf_header);
    m_rmf.header(rmf_header);

    // Set instrument name
    m_instrument = obs.instrument() + "OnOff";

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute ARF of On/Off observation
 *
 * @param[in] obs CTA observation.
 * @param[in] spatial Spatial source model.
 * @param[in] on On regions.
 *
 * @exception GException::invalid_argument
 *            No CTA response found in CTA observation.
 *
 * Computes the ARF for an On/Off observation by integration over the IRF
 * for the specified @p spatial source model over the @p on regions.
 *
 * All On regions contained in @p on are expected to be sky region maps,
 * i.e. of type GSkyRegionMap.
 ***************************************************************************/
void GCTAOnOffObservation::compute_arf(const GCTAObservation& obs,
                                       const GModelSpatial&   spatial,
                                       const GSkyRegions&     on)
{
    // Get reconstructed energy boundaries from on ARF
    const GEbounds& etrue = m_arf.ebounds();
    int             ntrue = etrue.size();

    // Continue only if there are ARF bins
    if (ntrue > 0) {

        // Get CTA IRF response
        const GCTAResponseIrf& rsp = gammalib::cta_rsp_irf(G_COMPUTE_ARF, obs);

        // Set dummy time
        const GTime time;

        // Save original energy dispersion application status
        bool save_edisp = rsp.apply_edisp();

        // Switch-off application of energy dispersion
        const_cast<GCTAResponseIrf&>(rsp).apply_edisp(false);

        // Loop over true energies
        for (int i = 0; i < ntrue; ++i) {

            // Get mean energy of bin
            GEnergy energy = etrue.elogmean(i);

            // Set source
            GSource source("", const_cast<GModelSpatial*>(&spatial), energy, time);

            // Initialize effective area for this bin
            m_arf[i] = 0.0;

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
                    GCTAInstDir pixdir = GCTAInstDir(on_map->map().inx2dir(pixidx));

                    // Get solid angle subtended by this pixel
                    double pixsolid = on_map->map().solidangle(pixidx);

                    // Set event
                    GCTAEventBin event;
                    event.dir(pixdir);
                    event.energy(energy);
                    event.time(time);

                    // Get ARF value. We need to devide by the deadtime
                    // correction since the IRF method multiplies with it.
                    double irf = rsp.irf(event, source, obs) / m_deadc;

                    // Add up effective area
                    m_arf[i] += irf * pixsolid;

                } // endfor: looped over all pixels in region map

            } // endfor: looped over all regions

        } // endfor: looped over true energies

        // Put back original energy dispersion application status
        const_cast<GCTAResponseIrf&>(rsp).apply_edisp(save_edisp);

	} // endif: there were energy bins

	// Return
	return;
}


/***********************************************************************//**
 * @brief Compute ARF of On/Off observation for a IRF with radius cut
 *
 * @param[in] obs CTA observation.
 * @param[in] spatial Spatial source model.
 * @param[in] on On regions.
 *
 * @exception GException::invalid_argument
 *            No CTA response found in CTA observation.
 *
 * Computes the ARF for an On/Off observation by computing the average
 * effective area over the @p on regions for the specified @p spatial source
 * model.
 *
 * All On regions contained in @p on are expected to be sky region maps,
 * i.e. of type GSkyRegionMap.
 ***************************************************************************/
void GCTAOnOffObservation::compute_arf_cut(const GCTAObservation& obs,
                                           const GModelSpatial&   spatial,
                                           const GSkyRegions&     on)
{
    // Get reconstructed energy boundaries from on ARF
    const GEbounds& etrue = m_arf.ebounds();
    int             ntrue = etrue.size();

    // Continue only if there are ARF bins
    if (ntrue > 0) {

        // Get CTA IRF response
        const GCTAResponseIrf& rsp = gammalib::cta_rsp_irf(G_COMPUTE_ARF, obs);

        // Get CTA observation pointing direction, zenith, and azimuth
        GCTAPointing obspnt  = obs.pointing();
        GSkyDir      obsdir  = obspnt.dir();
        double       zenith  = obspnt.zenith();
        double       azimuth = obspnt.azimuth();

        // Loop over true energies
        for (int i = 0; i < ntrue; ++i) {

            // Get mean energy of bin
            double logEtrue = etrue.elogmean(i).log10TeV();

            // Initialize effective area and solid angle sum for this bin
            m_arf[i]   = 0.0;
            double sum = 0.0;

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
                    m_arf[i] += rsp.aeff(theta, phi,
                                         zenith, azimuth,
                                         logEtrue) * pixsolid;

                    // Sum up solid angles
                    sum += pixsolid;

                } // endfor: looped over all pixels in region map

            } // endfor: looped over all regions

            // Divide effective area by solid angle so that ARF contains the
            // average effective area over the On region
            if (sum > 0.0) {
                m_arf[i] /= sum;
            }
            else {
                m_arf[i] = 0.0;
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
 * @param[in] models Models including background model.
 * @param[in] use_model_bkg Use model background.
 *
 * Computes the background rate in units of
 *           \f${\rm events} \, {\rm MeV}^{-1} \, {\rm s}^{-1}\f$
 * in the Off regions and stores the result as additional column with name
 * `BACKRESP` in the Off spectrum.
 *
 * All Off regions contained in @p off are expected to be sky region maps,
 * i.e. of type GSkyRegionMap.
 *
 * If @p use_model_bkg is @c true, the IRF background template will be used
 * for the computation, and `BACKRESP` will be computed using
 *
 * \f[
 *    {\tt BACKRESP}(E_{\rm reco}) = \sum_{\rm off} \sum_p
 *                                   {\tt BKG}_{{\rm off},p}(E_{\rm reco})
 *                                   \times \Omega_{{\rm off},p}
 * \f]
 *
 * where \f${\rm off}\f$ is the index of the Off region and \f$p\f$ is the
 * pixel in the Off region (note that each Off region is transformed into a
 * region map and \f$p\f$ indicates the pixels of this map).
 * \f${\tt BKG}_{{\rm off},p}(E_{\rm reco})\f$ is the background rate as
 * provided by the IRF background template in units of
 *  \f${\rm events} \, {\rm MeV}^{-1} \, {\rm s}^{-1} \, {\rm sr}^{-1}\f$
 * for a reconstructed energy \f$E_{\rm reco}\f$ and a pixel index \f$p\f$
 * in the Off region \f${\rm off}\f$.
 * \f$\Omega_{{\rm off},p}\f$ is the solid angle in units of \f${\rm sr}\f$
 * of the pixel index \f$p\f$ in the Off region \f${\rm off}\f$.
 *
 * If @p use_model_bkg is @c false, `BACKRESP` will be computed using
 *
 * \f[
 *    {\tt BACKRESP}(E_{\rm reco}) = \sum_{\rm off} \sum_p \Omega_{{\rm off},p}
 * \f]
 *
 * which actually assumes that the background rate is constant over the
 * field of view.
 *
 * @todo Integrate background rate over energy bin instead of computing the
 *       rate at the energy bin centre.
 ***************************************************************************/
void GCTAOnOffObservation::compute_bgd(const GCTAObservation& obs,
                                       const GSkyRegions&     off,
                                       const GModels&         models,
                                       const bool&            use_model_bkg)
{
    // Get reconstructed energies for the background rates
	const GEbounds& ereco = m_off_spec.ebounds();
    int             nreco = ereco.size();

    // Continue only if there are energy bins
    if (nreco > 0) {

		// Initialise background rates to zero
        std::vector<double> background(nreco, 0.0);

        // Get CTA observation pointing direction
        GCTAPointing obspnt = obs.pointing();

        // Continue only if the IRF background template should be used
        if (use_model_bkg) {

            // Loop over regions
            for (int k = 0; k < off.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* off_map =
                      static_cast<const GSkyRegionMap*>(off[k]);

                // Loop over pixels in Off region map and integrate
                // background rate
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

                        // Set event
                        GCTAEventAtom event(pixinstdir, ereco.elogmean(i), GTime());

                        // Get background rate in events/s/MeV
                        background[i] += models.eval(event, obs) * pixsolid;

                    } // endfor: looped over energy bins

                } // endfor: looped over all pixels in map

            } // endfor: looped over all regions

        } // endif: IRF background template was used

        // ... otherwise
        else {

            // Loop over regions
            for (int k = 0; k < off.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* off_map =
                      static_cast<const GSkyRegionMap*>(off[k]);

                // Loop over pixels in Off region map and integrate
                // background rate
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

                        // Get background rate in events/s/MeV
                        background[i] += pixsolid;

                    } // endfor: looped over energy bins

                } // endfor: looped over all pixels in map

            } // endfor: looped over all regions

        } // endelse: no IRF background template was used

        // Append background vector to Off spectrum
        m_off_spec.append("BACKRESP", background);

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
 * @param[in] models Models including background model.
 * @param[in] use_model_bkg Use model background.
 *
 * @exception GException::invalid_argument
 *            Observation does not contain relevant response or background
 *            information
 *
 * Compute the \f$\alpha\f$ parameters for all reconstructed energy bins.
 * The \f$\alpha\f$ parameter gives the ratio between the On and Off region
 * background acceptance multiplied by the On and Off region solid angles.
 *
 * If @p use_model_bkg is @c true, the IRF background template will be used
 * for the computation, and \f$\alpha(E_{\rm reco})\f$ is given by
 *
 * \f[
 *    \alpha(E_{\rm reco}) =
 *    \frac{\sum_{\rm on} \sum_p {\tt BKG}_{{\rm on},p}(E_{\rm reco})
 *          \times \Omega_{{\rm on},p}}
 *         {\sum_{\rm off} \sum_p {\tt BKG}_{{\rm off},p}(E_{\rm reco})
 *          \times \Omega_{{\rm off},p}}
 * \f]
 *
 * where the nominator sums over the On regions, indicated by the index
 * \f${\rm on}\f$, and the denominator sums over the Off regions, indicated
 * by the index \f${\rm off}\f$. Each On or Off region is defined by a
 * region sky map of type GSkyRegionMap, and the pixels of these maps
 * are index by \f$p\f$.
 * \f${\tt BKG}_{{\rm on/off},p}(E_{\rm reco})\f$ is the background rate as
 * provided by the IRF background template in units of
 *  \f${\rm events} \, {\rm MeV}^{-1} \, {\rm s}^{-1} \, {\rm sr}^{-1}\f$
 * for a reconstructed energy \f$E_{\rm reco}\f$ and a pixel index \f$p\f$
 * in the On or Off region \f${\rm on/off}\f$.
 * \f$\Omega_{{\rm on/off},p}\f$ is the solid angle in units of
 * \f${\rm sr}\f$ of the pixel index \f$p\f$ in the On or Off region
 * \f${\rm on/off}\f$.
 *
 * If @p use_model_bkg is @c false, the background acceptance is assumed
 * constant and hence cancels out, which makes the \f$\alpha\f$ parameter
 * independent of reconstructed energy \f$E_{\rm reco}\f$.
 * The \f$\alpha\f$ parameter is then given by
 *
 * \f[
 *    \alpha = \frac{\sum_{\rm on} \sum_p \Omega_{{\rm on},p}}
 *                  {\sum_{\rm off} \sum_p \Omega_{{\rm off},p}}
 * \f]
 *
 * @todo Compute alpha by integrating the background rate over the energy
 *       bins and not at the bin centre.
 ***************************************************************************/
void GCTAOnOffObservation::compute_alpha(const GCTAObservation& obs,
                                         const GSkyRegions&     on,
                                         const GSkyRegions&     off,
                                         const GModels&         models,
                                         const bool&            use_model_bkg)
{
    // Get reconstructed energy boundaries from RMF
	const GEbounds& ereco = m_rmf.emeasured();
    int             nreco = ereco.size();

    // Continue only if there are reconstructed energy bins
    if (nreco > 0) {

        // Get CTA observation pointing direction, zenith, and azimuth
        GCTAPointing obspnt = obs.pointing();

        // If IRF background templates shall be used then compute the
        // energy dependent alpha factors
        if (use_model_bkg) {

            // Loop over reconstructed energies
            for (int i = 0; i < nreco; ++i) {

                // Get mean log10 energy in TeV of bin
                //double logEreco = ereco.elogmean(i).log10TeV();

                // Initialise background rate totals
                double aon  = 0.0;
                double aoff = 0.0;

                // Loop over On regions
                for (int k = 0; k < on.size(); ++k) {

                    // Get pointer on sky region map
                    const GSkyRegionMap* on_map =
                          static_cast<const GSkyRegionMap*>(on[k]);

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

                        // Set event
                        GCTAEventAtom event(pixinstdir, ereco.elogmean(i), GTime());

                        // Add up acceptance
                        aon += models.eval(event, obs) * pixsolid;

                    } // endfor: looped over all pixels in map

                } // endfor: looped over regions

                // Loop over Off regions
                for (int k = 0; k < off.size(); ++k) {

                    // Get pointer on sky region map
                    const GSkyRegionMap* off_map =
                          static_cast<const GSkyRegionMap*>(off[k]);

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

                        // Set event
                        GCTAEventAtom event(pixinstdir, ereco.elogmean(i), GTime());

                        // Add up acceptance
                        aoff += models.eval(event, obs) * pixsolid;

                    } // endfor: looped over all pixels in map

                } // endfor: looped over all regions

                // Compute alpha for this energy bin
                double alpha = (aoff > 0.0) ? aon/aoff : 1.0;

                // Set background scaling in On spectra
                m_on_spec.backscal(i, alpha);

            } // endfor: looped over reconstructed energies

        } // endif: IRF background templates were used

        // ... otherwise compute energy independent alpha factor
        else {

            // Initialise background rate totals
            double aon  = 0.0;
            double aoff = 0.0;

            // Loop over On regions
            for (int k = 0; k < on.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* on_map =
                      static_cast<const GSkyRegionMap*>(on[k]);

                // Loop over pixels in On region map and integrate acceptance
                for (int j = 0; j < on_map->nonzero_indices().size(); ++j) {

                    // Get pixel index
                    int pixidx = on_map->nonzero_indices()[j];

                    // Get direction to pixel center
                    GSkyDir pixdir = on_map->map().inx2dir(pixidx);

                    // Translate sky direction into instrument direction
                    GCTAInstDir pixinstdir = obspnt.instdir(pixdir);

                    // Add up solid angle subtended by this pixel
                    aon += on_map->map().solidangle(pixidx);

                } // endfor: looped over all pixels in map

            } // endfor: looped over regions

            // Loop over Off regions
            for (int k = 0; k < off.size(); ++k) {

                // Get pointer on sky region map
                const GSkyRegionMap* off_map =
                      static_cast<const GSkyRegionMap*>(off[k]);

                // Loop over pixels in Off region map and integrate acceptance
                for (int j = 0; j < off_map->nonzero_indices().size(); ++j) {

                    // Get pixel index
                    int pixidx = off_map->nonzero_indices()[j];

                    // Get direction to pixel center
                    GSkyDir pixdir = off_map->map().inx2dir(pixidx);

                    // Translate sky direction into instrument direction
                    GCTAInstDir pixinstdir = obspnt.instdir(pixdir);

                    // Add up solid angle subtended by this pixel
                    aoff += off_map->map().solidangle(pixidx);

                } // endfor: looped over all pixels in map

            } // endfor: looped over all regions

            // Compute alpha for this energy bin
            double alpha = (aoff > 0.0) ? aon/aoff : 1.0;

            // Set background scaling in On spectra
            for (int i = 0; i < nreco; ++i) {
                m_on_spec.backscal(i, alpha);
            }

        } // endelse: computed energy independent alpha factor

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
 * @exception GException::invalid_argument
 *            Observation does not contain IRF response
 *
 * Compute the energy redistribution matrix for an On/Off observation. The
 * method requires that the RMF energy axes have been defined before.
 ***************************************************************************/
void GCTAOnOffObservation::compute_rmf(const GCTAObservation& obs,
                                       const GSkyRegions&     on)
{
    // Get true and reconstructed energy boundaries from RMF
    const GEbounds& etrue = m_rmf.etrue();
    const GEbounds& ereco = m_rmf.emeasured();
    int             ntrue = etrue.size();
    int             nreco = ereco.size();

    // Continue only if there are RMF bins
    if (ntrue > 0 && nreco > 0) {

        // Get CTA IRF response
        const GCTAResponseIrf& rsp = gammalib::cta_rsp_irf(G_COMPUTE_RMF, obs);

        // Get CTA observation pointing direction, zenith, and azimuth
        GCTAPointing obspnt  = obs.pointing();
        GSkyDir      obsdir  = obspnt.dir();
        double       zenith  = obspnt.zenith();
        double       azimuth = obspnt.azimuth();

        // Initialise RMF matrix
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
                    double aeff = rsp.aeff(theta, phi,
                                           zenith, azimuth,
                                           logEtrue);

                    // Loop over reconstructed energy
                    for (int ireco = 0; ireco < nreco; ++ireco) {

                        // Get RMF value
                        double value = rsp.edisp()->prob_erecobin(ereco.emin(ireco),
                                                                  ereco.emax(ireco),
                                                                  etrue.elogmean(itrue),
                                                                  theta);

                        // Update RMF value and weight
                        m_rmf(itrue, ireco)  += value * aeff;
                        weight(itrue, ireco) += aeff;

                    } // endfor: looped over reconstructed energy

                } // endfor: looped over true energy

            } // endfor: looped over all pixels in map

        } // endfor: looped over all regions

        // Normalise RMF matrix
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


/***********************************************************************//**
 * @brief Apply CTA observation energy boundaries to On/Off observation
 *
 * @param[in] obs CTA observation.
 *
 * Applies CTA energy boundaries to all histograms used for the On/Off
 * analysis in.
 *
 * For the PHA On and Off spectra, all bins are set to zero that do not fully
 * overlap with the CTA observation energy boundaries. Specifically,
 * partially overlapping bins are excluded. Some margin is applied that
 * effectively reduces the width of the PHA energy bin, which should cope
 * with any rounding errors.
 *
 * For the background response, stored in the Off spectrum, all bins are set
 * to zero that do not fully overlap with the CTA observation energy
 * boundaries.
 *
 * For the RMF, all reconstructued energy bins are set to zero that do not
 * fully overlap with the CTA observation energy boundaries. True energy
 * bins remain unchanged to properly account for energy migration.
 ***************************************************************************/
void GCTAOnOffObservation::apply_ebounds(const GCTAObservation& obs)
{
    // Set energy margin
    const GEnergy energy_margin(0.01, "GeV");

    // Get true and reconstructed energy boundaries from RMF
    const GEbounds& etrue = m_rmf.etrue();
    const GEbounds& ereco = m_rmf.emeasured();
    int             ntrue = etrue.size();
    int             nreco = ereco.size();

    // Get energy boundaries of observations
    GEbounds eobs = obs.ebounds();

    // Get reference to background response
    std::vector<double>& background = m_off_spec["BACKRESP"];

    // Apply energy boundaries in reconstructed energy
    for (int ireco = 0; ireco < nreco; ++ireco) {
        if (!eobs.contains(ereco.emin(ireco) + energy_margin,
                           ereco.emax(ireco) - energy_margin)) {
            m_on_spec[ireco]  = 0.0;
            m_off_spec[ireco] = 0.0;
            background[ireco] = 0.0;
            for (int itrue = 0; itrue < ntrue; ++itrue) {
                m_rmf(itrue, ireco) = 0.0;
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if effective area IRF has a radius cut
 *
 * @param[in] obs CTA observation.
 * @param[in] on On regions.
 * @return Radius cut value in degrees (zero for no radius cut).
 *
 * @exception GException::invalid_argument
 *            IRF has a radius cut that is incompatible with the On regions
 *
 * Checks if the effective area IRF has a radius cut. If a radius cut is
 * found the cut value is checked against the radii of the On regions. If
 * the On regions are not circular, or if they have a radius that differs
 * from the IRF cut value, an exception is thrown.
 ***************************************************************************/
double GCTAOnOffObservation::arf_rad_max(const GCTAObservation& obs,
                                         const GSkyRegions&     on) const
{
    // Initialise radius cut value
    double rad_max = 0.0;

    // Get pointer on CTA IRF response. Continue only if a valid response
    // was found
    const GCTAResponseIrf* rsp =
          dynamic_cast<const GCTAResponseIrf*>(obs.response());
    if (rsp != NULL) {

        // Get pointer to CTA 2D effective area. Continue only if a 2D
        // effective area was found
        const GCTAAeff2D* aeff = dynamic_cast<const GCTAAeff2D*>(rsp->aeff());
        if (aeff != NULL) {

            // Get cut radius. Continue only if cut radius is positive
            rad_max = aeff->rad_max();
            if (rad_max > 0.0) {

                // Verify that all On regions are circular regions and make
                // sure that they have, within some margin, the same radius
                // than the cut radius.
                for (int i = 0; i < on.size(); ++i) {

                    // Check that region is a circular region
                    const GSkyRegionCircle* region =
                          dynamic_cast<const GSkyRegionCircle*>(on[i]);
                    if (region == NULL) {
                        std::string msg = "An effective area IRF with a theta "
                                          "cut was specified, but the On region "
                                          "is not a circle, hence the "
                                          "consistency of the event selection "
                                          "could not be checked. Please specify "
                                          "a circular On region if an effective "
                                          "area with a theta cut is used.";
                        throw GException::invalid_argument(G_ARF_RADIUS_CUT, msg);
                    }

                    // Check that the cut radius is consistent
                    if (std::abs(region->radius()-rad_max) > 1.0e-3) {
                        std::string msg = "An effective area IRF with a theta "
                                          "cut of "+gammalib::str(rad_max)+ " deg "
                                          "was specified but an On region with "
                                          "a radius of "+
                                          gammalib::str(region->radius())+" deg "
                                          "was encountered. Please specify On "
                                          "regions with a radius of "+
                                          gammalib::str(rad_max)+" deg for this "
                                          "IRF.";
                        throw GException::invalid_argument(G_ARF_RADIUS_CUT, msg);
                    }

                } // endfor: looped over On regions

            } // endif: cut radius was positive

        } // endif: 2D effective area was found

    } // endif: CTA response was found

    // Return radius cut flag
    return rad_max;
}


/***********************************************************************
 * @brief Compute \f$N_{\gamma}\f$ value and model parameter gradients
 *
 * @param[in] models Model container.
 * @param[in] ibin Energy bin number.
 * @param[in,out] grad Model gradient vector.
 * @returns Predicted number of source events.
 *
 * @exception GException::invalid_value
 *            Source model is not a point source model
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

            // Spectral component (the useful one)
            GModelSpectral* spectral = sky->spectral();
            if (spectral != NULL)  {

                // Debug code
                #if defined(G_N_GAMMA_DEBUG)
                double rmf_sum = 0.0;
                #endif

                // Set instrument scale factor
                double scale = (sky->has_scales()) ? sky->scale(instrument()).value() : 1.0;

                // Loop over true energy bins
                for (int itrue = 0; itrue < m_arf.size(); ++itrue) {

                    // Get ARF value. Continue only if it is positive
                    double arf = m_arf[itrue];
                    if (arf <= 0.0) {
                        continue;
                    }

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
                    double norm_flux = arf * exposure * rmf * scale;
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

                    // Compute the parameter gradients for all model parameters
                    for (int k = 0; k < mptr->size(); ++k)  {
                        const GModelPar& par = (*mptr)[k];
                        if (par.is_free() && ipar+k < npars)  {
                            (*grad)[ipar+k] += par.factor_gradient() * norm_grad;
                        }
                    }

                } // endfor: looped over true energy bins

                // Debug code
                #if defined(G_N_GAMMA_DEBUG)
                std::cout << "sum(RMF) = " << rmf_sum << std::endl;
                #endif

            } // endif: pointer to spectral component was not NULL

            // Increment the parameter index
            ipar += mptr->size();

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

        // Get reconstructed energy bin mean and width
        GEnergy emean  = m_on_spec.ebounds().elogmean(ibin);
        double  ewidth = m_on_spec.ebounds().ewidth(ibin).MeV();

        // Perform log-log interpolation of background rate (events/MeV/s)
        // at reconstructed energy
        double background = m_off_spec["BACKRESP"][ibin];

        // Continue only if background rate is positive
        if (background > 0.0) {

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

                // Fall through if model is not a background component
                const GModelData* bgd = dynamic_cast<const GModelData*>(mptr);
                if (bgd == NULL) {
                    ipar += mptr->size();
                    continue;
                }

                // Extract model dependent spectral component
                const GModelSpectral* spectral = gammalib::cta_model_spectral(*bgd);

                // Get value from spectral component
                if (spectral != NULL)  {

                    // Determine the number of background events in model by
                    // computing the model normalization at the mean value of
                    // the energy bin and multiplying the normalisation with
                    // the number of background events. The eval() method needs
                    // a time in case that the spectral model has a time
                    // dependence. We simply use a dummy time here.
                    value += spectral->eval(emean, GTime(), true) * norm;

                    // Compute the parameter gradients for all model parameters
                    for (int k = 0; k < mptr->size(); ++k)  {
                        const GModelPar& par = (*mptr)[k];
                        if (par.is_free() && ipar+k < npars)  {
                            (*grad)[ipar+k] += par.factor_gradient() * norm;
                        }
                    }

                } // endif: pointer to spectral component was not NULL

                // Increment the parameter index
                ipar += mptr->size();

            } // endfor: looped over model components

        } // endif: background rate was positive

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

        // Compute and updated predicted number of events
        double nonpred = ngam + alpha * nbgd;
        *npred        += nonpred;

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if ((nbgd <= minmod) || (nonpred <= minmod)) {
            #if defined(G_LIKELIHOOD_DEBUG)
            n_small_model++;
            #endif
            continue;
        }

        // Now we have all predicted gamma and background events for
        // current energy bin. Update the log(likelihood)
        value += -non * std::log(nonpred) + nonpred - noff * std::log(nbgd) + nbgd;

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
 * Computes the log-likelihood value for the On/Off observation. The method
 * loops over the energy bins to update the function value, its derivatives
 * and the curvature matrix.
 *
 * The number of On counts \f$n_{\rm on}\f$ and Off counts \f$n_{\rm off}\f$
 * are taken from the On and Off spectra, the expected number of gamma-ray
 * events \f$\mu_s\f$ is computed from the spectral models of the relevant
 * components in the model container (spatial and temporal components
 * are ignored so far). See the N_gamma() method for details about the
 * model computations.
 *
 * The estimated number of background counts \f$\mu_b\f$ is derived based
 * on the measurement in the Off region by analytical minimization of
 * the Poisson likelihood, i.e., it is treated as a nuisance parameter.
 * See the wstat_value function for the derivation.
 ***********************************************************************/
double GCTAOnOffObservation::likelihood_wstat(const GModels& models,
                                              GVector*       gradient,
                                              GMatrixSparse* curvature,
                                              double*        npred) const
{
    // Debug option: dump header
    #if defined(G_LIKELIHOOD_DEBUG)
    std::cout << "GCTAOnOffObservation::likelihood_wstat: enter" << std::endl;
    #endif

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

        // Get number of gamma events (and corresponding spectral model
        // gradients)
        double ngam = N_gamma(models, i, &sky_grad);

        // Initialise variables for likelihood computation
        double nonpred     = 0.0;
        double nbgd        = 0.0;
        double dlogLdsky   = 0.0;
        double d2logLdsky2 = 0.0;

        // Perform likelihood profiling and derive likelihood value
        double logL = wstat_value(non, noff, alpha, ngam,
                                  nonpred, nbgd, dlogLdsky, d2logLdsky2);

        // Update global log-likelihood and Npred
        value  += logL;
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
                (*gradient)[j] += dlogLdsky * sky_grad[j];

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

    // Dump statistics
    #if defined(G_LIKELIHOOD_DEBUG)
    std::cout << "Number of parameters: " << npars << std::endl;
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
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
    std::cout << "GCTAOnOffObservation::optimizer::likelihood_wstat: CPU usage = "
              << t_elapse << " sec" << std::endl;
    #endif

    // Debug option: dump trailer
    #if defined(G_LIKELIHOOD_DEBUG)
    std::cout << "GCTAOnOffObservation::likelihood_wstat: exit" << std::endl;
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood value in energy bin for On/Off analysis
 * by profiling over the number of background counts
 *
 * @param[in] non number of On counts
 * @param[in] noff number of Off counts
 * @param[in] alpha background scaling rate
 * @param[in] ngam number of predicted gamma-ray counts
 * @param[in,out] nonpred number of predicted On counts
 * @param[in,out] nbgd number of predicted background counts 
 * @param[in,out] dlogLdsky first derivative of log-like w.r.t. sky pars
 * @param[in,out] d2logLdsky2 second derivative of log-like w.r.t. sky pars
 * @return Log-likelihood value.
 *
 * Computes the log-likelihood value for the On/Off observation in an energy bin
 * by treating the number of background counts as nuisance parameter. The method
 * performs an analytical minimisation of the Poisson likelihood and updates 
 * the relevant values.
 * In the general case, the log-likelihood function is computed using
 *
 * \f[
 *   L = \mu_s + (1+\alpha) \mu_b - n_{\rm on} - n_{\rm off} -
 *       n_{\rm on} \left( \ln(\mu_s + \alpha \mu_b) - \ln(n_{\rm on})
 *       \right) -
 *       n_{\rm off} \left( \ln(\mu_b) - \ln(n_{\rm off}) \right)
 * \f]
 *
 * where
 *
 * \f[
 *    \mu_b = \frac{C+D}{2\alpha(1+\alpha)}
 * \f]
 *
 * are the estimated number of background counts with
 *
 * \f[
 *    C = \alpha (n_{\rm on} + n_{\rm off}) - (1 + \alpha) \mu_s
 * \f]
 *
 * and
 *
 * \f[
 *    D = \sqrt{C^2 + 4 (1 + \alpha) \, \alpha \, n_{\rm off} \, \mu_s}
 * \f]
 *
 * \
 *
 * The first derivative of the log-likelihood function is given by
 *
 * \f[
 *    \frac{\delta L}{\delta \mu_s} =
 *    1 + (1 + \alpha) \frac{\delta \mu_b}{\delta \mu_s} -
 *    \frac{n_{\rm on}}{\mu_s + \alpha \mu_b} \left( 1 + \alpha
 *    \frac{\delta \mu_b}{\delta \mu_s} \right) -
 *    \frac{n_{\rm off}}{\mu_b} \frac{\delta \mu_b}{\delta \mu_s}
 * \f]
 *
 * with
 *
 * \f[
 *    \frac{\delta \mu_b}{\delta \mu_s} =
 *    \frac{n_{\rm off} - C}{D} - \frac{1}{2 \alpha}
 * \f]
 *
 * The second derivative of the log-likelihood function is given by
 *
 * \f[
 *    \frac{\delta^2 L}{\delta \mu_s^2} =
 *    \frac{n_{\rm on}}{(\mu_s + \alpha \mu_b)^2} \left( 1 + \alpha
 *      \frac{\delta \mu_b}{\delta \mu_s} \right) +
 *    \frac{\delta^2 \mu_b}{\delta \mu_s^2} \left( (1 + \alpha) -
 *      \frac{\alpha n_{\rm on}}{\mu_s + \alpha \mu_b} -
 *      \frac{n_{\rm off}}{\mu_b} \right) +
 *    \frac{\delta \mu_b}{\delta \mu_s} \left(
 *      \frac{\alpha n_{\rm on}}{(\mu_s + \alpha \mu_b)^2} \left(
 *        1 + \alpha \frac{\delta \mu_b}{\delta \mu_s} \right) +
 *      \frac{n_{\rm off}}{\mu_b^2} \frac{\delta \mu_b}{\delta \mu_s}
 *      \right)
 * \f]
 *
 * with
 *
 * \f[
 *    \frac{\delta^2 \mu_b}{\delta \mu_s^2} =
 *    \frac{-1}{2 \alpha} \left(
 *    \frac{1}{D} \frac{\delta C}{\delta \mu_s} +
 *    \frac{2 \alpha \, n_{\rm off} - C}{D^2} \frac{\delta D}{\delta \mu_s}
 *    \right)
 * \f]
 *
 * where
 *
 * \f[
 *    \frac{\delta C}{\delta \mu_s} = -(1 + \alpha)
 * \f]
 *
 * and
 *
 * \f[
 *    \frac{\delta D}{\delta \mu_s} =
 *    \frac{4 (1 + \alpha) \, \alpha \, n_{\rm off} - 2 \, C \, (1 + \alpha)}
 *         {2D}
 * \f]
 *
 * In the special case \f$n_{\rm on}=n_{\rm off}=0\f$ the formal
 * background estimate becomes negative, hence we set \f$\mu_b=0\f$ and
 * the log-likelihood function becomes
 *
 * \f[
 *    L = \mu_s
 * \f]
 *
 * the first derivative
 *
 * \f[
 *    \frac{\delta L}{\delta \mu_s} = 1
 * \f]
 *
 * and the second derivative
 *
 * \f[
 *    \frac{\delta^2 L}{\delta \mu_s^2} = 0
 * \f]
 *
 * In the special case \f$n_{\rm on}=0\f$ and \f$n_{\rm off}>0\f$
 * the log-likelihood function becomes
 *
 * \f[
 *    L = \mu_s + n_{\rm off} \ln(1 + \alpha)
 * \f]
 *
 * the background estimate
 *
 * \f[
 *    \mu_b = \frac{n_{\rm off}}{1+\alpha}
 * \f]
 *
 * the first derivative
 *
 * \f[
 *    \frac{\delta L}{\delta \mu_s} = 1
 * \f]
 *
 * and the second derivative
 *
 * \f[
 *    \frac{\delta^2 L}{\delta \mu_s^2} = 0
 * \f]
 *
 * In the special case \f$n_{\rm on}>0\f$ and \f$n_{\rm off}=0\f$
 * the background estimate becomes
 *
 * \f[
 *    \mu_b = \frac{n_{\rm on}}{1+\alpha} - \frac{\mu_s}{\alpha}
 * \f]
 *
 * which is positive for
 *
 * \f[
 *    \mu_s < n_{\rm on} \frac{\alpha}{1+\alpha}
 * \f]
 *
 * For positive \f$\mu_b\f$ the log-likelihood function is given by
 *
 * \f[
 *    L = -\frac{\mu_s}{\alpha}
 *        - n_{\rm on} \ln \left(\frac{\alpha}{1 + \alpha} \right)
 * \f]
 *
 * the first derivative
 *
 * \f[
 *    \frac{\delta L}{\delta \mu_s} = -\frac{1}{\alpha}
 * \f]
 *
 * and the second derivative
 *
 * \f[
 *    \frac{\delta^2 L}{\delta \mu_s^2} = 0
 * \f]
 *
 * For negative \f$\mu_b\f$ we set \f$\mu_b=0\f$ and the log-likelihood
 * function becomes
 *
 * \f[
 *    L = \mu_s - n_{\rm on} -
 *        n_{\rm on} \left( \ln(\mu_s) - \ln(n_{\rm on}) \right)
 * \f]
 *
 * the first derivative
 *
 * \f[
 *    \frac{\delta L}{\delta \mu_s} = 1 - \frac{n_{\rm on}}{\mu_s}
 * \f]
 *
 * and the second derivative
 *
 * \f[
 *    \frac{\delta^2 L}{\delta \mu_s^2} = \frac{1}{\mu_s^2}
 * \f]
 *
 * Note that the fit results may be biased and the statistical errors
 * overestimated if for some bins \f$n_{\rm on}=0\f$ and/or
 * \f$n_{\rm off}=0\f$ (i.e. if the special cases are encountered).
 ***********************************************************************/
double GCTAOnOffObservation::wstat_value(const double& non,
                                         const double& noff,
                                         const double& alpha,
                                         const double& ngam,
                                         double&       nonpred,
                                         double&       nbgd,
                                         double&       dlogLdsky,
                                         double&       d2logLdsky2) const
{
    // Initialise log-likelihood value
    double logL;

    // Precompute some values
    double alphap1  = alpha + 1.0;
    double alpharat = alpha / alphap1;

    // Calculate number of background events, profile likelihood value
    // and likelihood derivatives
    // Special case noff = 0
    if (noff == 0.0) {

        // Case A: non = 0. In this case nbgd < 0 hence we set nbgd = 0
        if (non == 0.0) {
            nbgd        = 0;
            nonpred     = ngam;
            logL        = ngam;
            dlogLdsky   = 1.0;
            d2logLdsky2 = 0.0;
        }

        // Case B: nbgd is positive
        else if (ngam < non * alpharat) {
            nbgd        = non / alphap1 - ngam / alpha;
            nonpred     = ngam + alpha * nbgd;
            logL        = -ngam / alpha - non * std::log(alpharat);
            dlogLdsky   = -1.0/alpha;
            d2logLdsky2 = 0.0;
        }

        // Case C: nbgd is zero or negative, hence set nbgd = 0
        else {
            nbgd        = 0;
            nonpred     = ngam;
            logL        = ngam + non * (std::log(non) - std::log(ngam) - 1.0);
            dlogLdsky   = 1.0 - non / ngam;
            d2logLdsky2 = non / (ngam * ngam);
        }

    } // endif: noff = 0

    // Special case non = 0
    else if (non == 0.0) {
        nbgd        = noff / alphap1;
        nonpred     = ngam + alpharat * noff;
        logL        = ngam + noff * std::log(alphap1);
        dlogLdsky   = 1.0;
        d2logLdsky2 = 0.0;
    } // endif: non = 0

    // General case
    else {

        // Compute log-likelihood value
        double alphat2 = 2.0 * alpha;
        double C       = alpha * (non + noff) - alphap1 * ngam;
        double D       = std::sqrt(C*C + 4.0 * alpha * alphap1 * noff * ngam);
        nbgd           = (C + D) / (alphat2 * alphap1);
        nonpred        = ngam + alpha * nbgd;
        logL           = ngam + alphap1 * nbgd - non - noff -
                         non * (std::log(nonpred) - std::log(non)) -
                         noff * (std::log(nbgd) - std::log(noff));

        // Compute derivatives
        double f0     = alphat2 * noff - C;
        double dCds   = -alphap1;
        double dDds   = (C * dCds + 2.0 * alphap1 * alpha * noff) / D;
        double dbds   = (f0 / D - 1.0) / alphat2;
        double d2bds2 = (-dCds / D - f0 / (D*D) * dDds) / alphat2;
        double f1     = alphap1 - alpha*non/nonpred - noff/nbgd;
        double f2     = nonpred * nonpred;
        double dpds   = 1.0 + alpha * dbds;
        double f3     = non / f2 * dpds;
        dlogLdsky     = 1.0 - non / nonpred + dbds * f1;
        d2logLdsky2   = f3 + d2bds2 * f1 +
                        dbds * (alpha * f3 + noff / (nbgd*nbgd) * dbds);

    } // endelse: general case

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(logL)    || gammalib::is_infinite(logL)    ||
        gammalib::is_notanumber(nonpred) || gammalib::is_infinite(nonpred) ||
        gammalib::is_notanumber(nbgd)    || gammalib::is_infinite(nbgd)) {
        std::cout << "*** ERROR: GCTAOnOffObservation::wstat_value";
        std::cout << "(noff=" << noff;
        std::cout << ", alpha=" << alpha;
        std::cout << ", ngam=" << ngam << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (logL=" << logL;
        std::cout << ", nonpred=" << nonpred;
        std::cout << ", nbgd=" << nbgd;
        std::cout << ", dlogLdsky=" << dlogLdsky;
        std::cout << ", d2logLdsky2=" << d2logLdsky2;
        std::cout << " )" << std::endl;
    }
    #endif

    // Return
    return logL;
}
