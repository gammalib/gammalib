/***************************************************************************
 *       GCTAModelCubeBackground.cpp - CTA cube background model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Michael Mayer                               *
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
 * @file GCTAModelCubeBackground.cpp
 * @brief CTA cube background model class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GIntegral.hpp"
#include "GCTAModelCubeBackground.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialDiffuseCube.hpp"
#include "GCTAObservation.hpp"
#include "GCTAPointing.hpp"
#include "GCTAInstDir.hpp"
#include "GCTARoi.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"
#include "GCTAResponseTable.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelCubeBackground g_cta_model_background_seed;
const GModelRegistry      g_cta_model_background_registry(&g_cta_model_background_seed);

/* __ Method name definitions ____________________________________________ */
#define G_SET_SPATIAL "GCTAModelCubeBackground::set_spatial(GCTAObservation&, "\
                                            "std::string&, int&, int&, int&)"
#define G_EVAL        "GCTAModelCubeBackground::eval(GEvent&, GObservation&)"
#define G_EVAL_GRADIENTS   "GCTAModelCubeBackground::eval_gradients(GEvent&,"\
                                                            " GObservation&)"
#define G_NPRED            "GCTAModelCubeBackground::npred(GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC              "GCTAModelCubeBackground::mc(GObservation&, GRan&)"
#define G_XML_SPATIAL    "GCTAModelCubeBackground::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL  "GCTAModelCubeBackground::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL  "GCTAModelCubeBackground::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_NPRED_CACHE      //!< Use Npred cache in npred_diffuse method
#define G_NPRED_AROUND_ROI        //!< Perform Npred integration around ROI

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_NPRED                       //!< Debug npred_diffuse method
//#define G_DUMP_MC                                  //!< Dump MC information


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty CTA cube background model.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(void) : GModelData()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA background model from the information that is found in a
 * XML element. Please refer to the method GCTAModelCubeBackground::read to
 * learn more about the information that is expected in the XML element.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GXmlElement& xml) :
                     GModelData(xml)
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from spatial and spectral components
 *
 * @param[in] spatial Spatial background model component.
 * @param[in] spectral Spectral model component.
 *
 * Constructs a CTA background model from a spatial and a spectral
 * model component. The temporal component is assumed to be constant.
 * Please refer to the classes GModelSpatial and GModelSpectral to learn
 * more about the definition of the spatial and spectral components.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GModelSpatial&  spatial,
                                                 const GModelSpectral& spectral) :
                         GModelData()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_spatial  = spatial.clone();
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model CTA background model.
 *
 * Constructs a CTA background model by copying information from an
 * existing model. Note that the copy is a deep copy, so the original object
 * can be destroyed after the copy without any loss of information.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GCTAModelCubeBackground& model) :
                         GModelData(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Construct from filename and observation
 *
 * @param[in] obs CTA observation.
 * @param[in] filename Background table filename.
 * @param[in] spectral Spectral model component.
 * @param[in] nx_sky Number of bins in x direction (optional).
 * @param[in] ny_sky Number of bins in y direction (optional).
 * @param[in] n_energy Number of energy bins (optional).
 *
 * Constructs a CTA background model from a spatial and a spectral
 * model component. The temporal component is assumed to be constant.
 * Please refer to the classes GModelSpatial and GModelSpectral to learn
 * more about the definition of the spatial and spectral components.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GCTAObservation& obs,
                                                 const std::string&     filename,
                                                 const GModelSpectral&  spectral,
                                                 const int&             nx_sky,
                                                 const int&             ny_sky,
                                                 const int&             n_energy) :
                         GModelData()
{
    // Initialise private members for clean destruction
	init_members();

	// Clone spectral model
	m_spectral = spectral.clone();

	// Create spatial cube from background file
	set_spatial(obs, filename, nx_sky, ny_sky, n_energy);

	// Set the temporal model
	GModelTemporalConst temporal;
	m_temporal = temporal.clone();

	// Set parameter pointers
	set_pointers();

	// Return
	return;
}



/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys a CTA background model.
 ***************************************************************************/
GCTAModelCubeBackground::~GCTAModelCubeBackground(void)
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
 * @param[in] model CTA background model.
 *
 * Assigns the information from a CTA background model to the actual
 * object. Note that a deep copy of the information is performed, so the
 * original object can be destroyed after the assignment without any loss of
 * information.
 ***************************************************************************/
GCTAModelCubeBackground& GCTAModelCubeBackground::operator=(const GCTAModelCubeBackground& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelData::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members (this method also sets the parameter pointers)
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
 * @brief Clear instance
 *
 * Resets the object to a clean initial state. All information that resided
 * in the object will be lost.
 ***************************************************************************/
void GCTAModelCubeBackground::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelData::free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    this->GModelData::init_members();
    init_members();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of CTA background model.
 *
 * Clone a CTA background model. Cloning performs a deep copy of the
 * information, so the original object can be destroyed after cloning without
 * any loss of information.
 ***************************************************************************/
GCTAModelCubeBackground* GCTAModelCubeBackground::clone(void) const
{
    return new GCTAModelCubeBackground(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            No CTA instrument direction found in event.
 *
 * Evaluates the CTA background model which is a factorization of a
 * spatial, spectral and temporal model component. This method also applies
 * a deadtime correction factor, so that the normalization of the model is
 * a real rate (counts/exposure time).
 *
 * @todo Add bookkeeping of last value and evaluate only if argument 
 *       changed
 ***************************************************************************/
double GCTAModelCubeBackground::eval(const GEvent& event,
                                     const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Extract CTA instrument direction
    const GCTAInstDir* dir  = dynamic_cast<const GCTAInstDir*>(&(event.dir()));
    if (dir == NULL) {
        std::string msg = "No CTA instrument direction found in event.";
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Create a Photon from the event.
    // We need the GPhoton to evaluate the spatial model.
    // For the background, GEvent and GPhoton are identical
    // since the IRFs are not folded in
    GPhoton ev_photon(dir->dir(), event.energy(), event.time());

    // Evaluate function and gradients
    double spat = (spatial() != NULL)
                  ? spatial()->eval(ev_photon) : 1.0;
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time()) : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time()) : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // Apply deadtime correction
    value *= obs.deadc(event.time());

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            No CTA instrument direction found in event.
 *
 * Evaluates tha CTA background model and parameter gradients. The CTA
 * background model is a factorization of a spatial, spectral and
 * temporal model component. This method also applies a deadtime correction
 * factor, so that the normalization of the model is a real rate
 * (counts/exposure time).
 *
 * @todo Add bookkeeping of last value and evaluate only if argument 
 *       changed
 ***************************************************************************/
double GCTAModelCubeBackground::eval_gradients(const GEvent& event,
                                               const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Extract CTA instrument direction
    const GCTAInstDir* dir  = dynamic_cast<const GCTAInstDir*>(&(event.dir()));
    if (dir == NULL) {
        std::string msg = "No CTA instrument direction found in event.";
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Create a Photon from the event
    // We need the photon to evaluate the spatial model
    // For the background, GEvent and GPhoton are identical
    // since the IRFs are not folded in
    GPhoton photon = GPhoton(dir->dir(), event.energy(),event.time());

    // Evaluate function and gradients
    double spat = (spatial() != NULL)
                  ? spatial()->eval_gradients(photon) : 1.0;
    double spec = (spectral() != NULL)
                  ? spectral()->eval_gradients(event.energy(), event.time()) : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval_gradients(event.time()) : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // Apply deadtime correction
    double deadc = obs.deadc(event.time());
    value       *= deadc;

    // Multiply factors to spatial gradients
    if (spatial() != NULL) {
        double fact = spec * temp * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < spatial()->size(); ++i)
                (*spatial())[i].factor_gradient( (*spatial())[i].factor_gradient() * fact );
        }
    }

    // Multiply factors to spectral gradients
    if (spectral() != NULL) {
        double fact = spat * temp * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < spectral()->size(); ++i)
                (*spectral())[i].factor_gradient( (*spectral())[i].factor_gradient() * fact );
        }
    }

    // Multiply factors to temporal gradients
    if (temporal() != NULL) {
        double fact = spat * spec * deadc;
        if (fact != 1.0) {
            for (int i = 0; i < temporal()->size(); ++i)
                (*temporal())[i].factor_gradient( (*temporal())[i].factor_gradient() * fact );
        }
    }

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated data model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated model.
 *
 * @exception GException::invalid_argument
 *            No CTA event list found in observation.
 *            No CTA pointing found in observation.
 *
 * Spatially integrates the data model for a given measured event energy and
 * event time. This method also applies a deadtime correction factor, so that
 * the normalization of the model is a real rate (counts/exposure time).
 ***************************************************************************/
double GCTAModelCubeBackground::npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const
{
    // Initialise result
    double npred     = 0.0;
    bool   has_npred = false;

    // Build unique identifier
    std::string id = obs.instrument() + "::" + obs.id();

    // Check if Npred value is already in cache
    #if defined(G_USE_NPRED_CACHE)
    if (!m_npred_names.empty()) {

        // Search for unique identifier, and if found, recover Npred value
		// and break
		for (int i = 0; i < m_npred_names.size(); ++i) {
			if (m_npred_names[i] == id && m_npred_energies[i] == obsEng) {
				npred     = m_npred_values[i];
				has_npred = true;
				#if defined(G_DEBUG_NPRED)
				std::cout << "GCTAModelCubeBackground::npred:";
				std::cout << " cache=" << i;
				std::cout << " npred=" << npred << std::endl;
				#endif
				break;
			}
		}

    } // endif: there were values in the Npred cache
    #endif

    // Continue only if no Npred cache value was found
    if (!has_npred) {

        // Evaluate only if model is valid
        if (valid_model()) {

            // Get CTA event list
			const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
            if (events == NULL) {
                std::string msg = "No CTA event list found in observation.\n" +
                                  obs.print();
                throw GException::invalid_argument(G_NPRED, msg);
            }

            #if !defined(G_NPRED_AROUND_ROI)
			// Get CTA pointing direction
			GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing());
            if (pnt == NULL) {
                std::string msg = "No CTA pointing found in observation.\n" +
                                  obs.print();
                throw GException::invalid_argument(G_NPRED, msg);
            }
            #endif

            // Get reference to ROI centre
            const GSkyDir& roi_centre = events->roi().centre().dir();

			// Get ROI radius in radians
			double roi_radius = events->roi().radius() * gammalib::deg2rad;

			// Get distance from ROI centre in radians
            #if defined(G_NPRED_AROUND_ROI)
			double roi_distance = 0.0;
            #else
			double roi_distance = roi_centre.dist(pnt->dir());
            #endif

			// Initialise rotation matrix to transform from ROI system to
            // celestial coordinate system
			GMatrix ry;
			GMatrix rz;
			ry.eulery(roi_centre.dec_deg() - 90.0);
			rz.eulerz(-roi_centre.ra_deg());
			GMatrix rot = (ry * rz).transpose();

			// Compute position angle of ROI centre with respect to model
			// centre (radians)
            #if defined(G_NPRED_AROUND_ROI)
            double omega0 = 0.0;
            #else
			double omega0 = pnt->dir().posang(events->roi().centre().dir());
            #endif

			// Setup integration function
			GCTAModelCubeBackground::npred_roi_kern_theta integrand(spatial(),
                                                                    obsEng,
                                                                    obsTime,
                                                                    rot,
                                                                    roi_radius,
                                                                    roi_distance,
                                                                    omega0);

			// Setup integrator
			GIntegral integral(&integrand);
			integral.eps(1.0e-3);

			// Setup integration boundaries
            #if defined(G_NPRED_AROUND_ROI)
			double rmin = 0.0;
			double rmax = roi_radius;
            #else
			double rmin = (roi_distance > roi_radius) ? roi_distance-roi_radius : 0.0;
			double rmax = roi_radius + roi_distance;
            #endif

			// Spatially integrate spatial component
			npred = integral.romberg(rmin, rmax);

	        // Store result in Npred cache
	        #if defined(G_USE_NPRED_CACHE)
	        m_npred_names.push_back(id);
	        m_npred_energies.push_back(obsEng);
	        m_npred_times.push_back(obsTime);
	        m_npred_values.push_back(npred);
	        #endif

	        // Debug: Check for NaN
	        #if defined(G_NAN_CHECK)
	        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
	            std::cout << "*** ERROR: GCTAModelCubeBackground::npred:";
	            std::cout << " NaN/Inf encountered";
	            std::cout << " (npred=" << npred;
	            std::cout << ", roi_radius=" << roi_radius;
	            std::cout << ")" << std::endl;
	        }
	        #endif

        } // endif: model was valid

    } // endif: Npred computation required

	// Multiply in spectral and temporal components
	npred *= spectral()->eval(obsEng, obsTime);
	npred *= temporal()->eval(obsTime);

	// Apply deadtime correction
	npred *= obs.deadc(obsTime);

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return Pointer to list of simulated events (needs to be de-allocated by
 *         client)
 *
 * @exception GException::invalid_argument
 *            No CTA event list found in observation.
 *
 * Draws a sample of events from the background model using a Monte
 * Carlo simulation. The pointing information, the energy boundaries and the
 * good time interval for the sampling will be extracted from the observation
 * argument that is passed to the method. The method also requires a random
 * number generator of type GRan which is passed by reference, hence the
 * state of the random number generator will be changed by the method.
 *
 * The method also applies a deadtime correction using a Monte Carlo process,
 * taking into account temporal deadtime variations. For this purpose, the
 * method makes use of the time dependent GObservation::deadc method.
 ***************************************************************************/
GCTAEventList* GCTAModelCubeBackground::mc(const GObservation& obs, GRan& ran) const
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Extract event list to access the ROI, energy boundaries and GTIs
        const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
        if (events == NULL) {
            std::string msg = "No CTA event list found in observation.\n" +
                              obs.print();
            throw GException::invalid_argument(G_MC, msg);
        }

        // Get simulation region
        const GCTARoi&  roi     = events->roi();
        const GEbounds& ebounds = events->ebounds();
        const GGti&     gti     = events->gti();

        // Set simulation region for result event list
        list->roi(roi);
        list->ebounds(ebounds);
        list->gti(gti);

        // Loop over all energy boundaries
        for (int ieng = 0; ieng < ebounds.size(); ++ieng) {

            // Initialise de-allocation flag
            bool free_spectral = false;

            // Set pointer to spectral model
            GModelSpectral* spectral = m_spectral;

            // If the spectral model is a diffuse cube then create a node
            // function spectral model that is the product of the diffuse
            // cube node function and the spectral model evaluated at the
            // energies of the node function
            GModelSpatialDiffuseCube* cube =
                dynamic_cast<GModelSpatialDiffuseCube*>(m_spatial);
            if (cube != NULL) {

			   // Set MC simulation cone based on ROI
			   cube->set_mc_cone(roi.centre().dir(), roi.radius());

			   // Allocate node function to replace the spectral component
			   GModelSpectralNodes* nodes = new GModelSpectralNodes(cube->spectrum());
			   for (int i = 0; i < nodes->nodes(); ++i) {
				   GEnergy energy    = nodes->energy(i);
				   double  intensity = nodes->intensity(i);
				   double  norm      = m_spectral->eval(energy, events->tstart());
				   nodes->intensity(i, norm*intensity);
			   }

			   // Signal that node function needs to be de-allocated later
			   free_spectral = true;

			   // Set the spectral model pointer to the node function
			   spectral = nodes;

            } // endif: spatial model was a diffuse cube

            // Compute the background rate in model within the energy boundaries
            // from spectral component (units: cts/s).
            // Note that the time here is ontime. Deadtime correction will be done
            // later.
            double rate = spectral->flux(ebounds.emin(ieng), ebounds.emax(ieng));

            // Debug option: dump rate
            #if defined(G_DUMP_MC)
            std::cout << "GCTAModelCubeBackground::mc(\"" << name() << "\": ";
            std::cout << "rate=" << rate << " cts/s)" << std::endl;
            #endif

            // Loop over all good time intervals
            for (int itime = 0; itime < gti.size(); ++itime) {

                // Get Monte Carlo event arrival times from temporal model
                GTimes times = m_temporal->mc(rate,
                                              gti.tstart(itime),
                                              gti.tstop(itime),
                                              ran);

                // Get number of events
                int n_events = times.size();

                // Reserve space for events
                if (n_events > 0) {
                    list->reserve(n_events);
                }

                // Loop over events
                for (int i = 0; i < n_events; ++i) {

                    // Apply deadtime correction
                    double deadc = obs.deadc(times[i]);
                    if (deadc < 1.0) {
                        if (ran.uniform() > deadc) {
                            continue;
                        }
                    }

                    // Get Monte Carlo event energy from spectral model
                    GEnergy energy = spectral->mc(ebounds.emin(ieng),
                                                  ebounds.emax(ieng),
                                                  times[i],
                                                  ran);

                    // Get Monte Carlo event direction from spatial model
                    GSkyDir dir = spatial()->mc(energy, times[i], ran);

                    // Allocate event
                    GCTAEventAtom event;

                    // Set event attributes
                    event.dir(GCTAInstDir(dir));
                    event.energy(energy);
                    event.time(times[i]);

                    // Append event to list if it falls in ROI
                    if (events->roi().contains(event)) {
                        list->append(event);
                    }

                } // endfor: looped over all events

            } // endfor: looped over all GTIs

            // Free spectral model if required
            if (free_spectral) delete spectral;

        } // endfor: looped over all energy boundaries

    } // endif: model was valid

    // Return
    return list;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * The model is composed of a spectrum component ('spectrum'), a spatial
 * component ('spatialModel'), and, optionally, of a temporal component
 * ('temporalModel'). If no temporal component is found a constant model is
 * assumed.
 *
 * The method also handles a instrumental background models that are defined
 * using information found in the "BACKGROUND" extension of the instrument
 * response functions. For this purpose a special spatialModel of type
 * "Instrumental" has been implemented, following the format
 *
 *     <spatialModel type="Instrumental" instrument="CTA">
 *       <parameter name="EventList"           file="events.fits"/>
 *       <parameter name="EffectiveArea"       file=""/>
 *       <parameter name="PointSpreadFunction" file=""/>
 *       <parameter name="EnergyDispersion"    file=""/>
 *       <parameter name="Background"          file="background.fits"/>
 *     </spatialModel>
 *
 * The format is identical to the format used for an observation definition,
 * and the XML element is in fact parsed using the GCTAObservation::read
 * method. This is for simplicity for the moment, yet the instrument response
 * related parameters are in fact not used. It should also be noted that the
 * actual code requires that the model identifier is set to exactly the same
 * identifier that is used for the corresponding CTA observation. This is
 * needed because the source models and the observation container are not
 * necessarily linked, hence there is no general way to pass this information
 * directly.
 ***************************************************************************/
void GCTAModelCubeBackground::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Initialise XML elements
    const GXmlElement* spatial  = NULL;
    const GXmlElement* spectral = NULL;
    const GXmlElement* temporal = NULL;

    // Get pointers on spectrum and spatial model
    spatial  = xml.element("spatialModel", 0);
    spectral = xml.element("spectrum", 0);

    // Handle special case of instrumental background model
    std::string type = spatial->attribute("type");
    if (type == "Instrumental") {

        // Create CTA observation from information in XML file
        GCTAObservation obs;
        obs.read(*spatial);

        // Set observation ID from id attribute
        obs.id(xml.attribute("id"));

        // Set spatial model
        set_spatial(obs, obs.m_bgdfile);

        // Set spectral model
        m_spectral = xml_spectral(*spectral);

    } // endif: handled special case of instrumental background

    // ... otherwise clone spatial and spectral models
    else {
        m_spatial  = xml_spatial(*spatial);
        m_spectral = xml_spectral(*spectral);
    }

    // Optionally get temporal model
    if (xml.elements("temporalModel")) {
        temporal   = xml.element("temporalModel", 0);
        m_temporal = xml_temporal(*temporal);
    }
    else {
        GModelTemporalConst constant;
        m_temporal = constant.clone();
    }

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set observation identifiers
    ids(xml.attribute("id"));

    // Check flag if TS value should be computed
    bool tscalc = (xml.attribute("tscalc") == "1") ? true : false;

    // Set flag if TS value should be computed
    this->tscalc(tscalc);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelCubeBackground::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = xml.element("source", k);
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (spatial()  != NULL) src->append(GXmlElement("spatialModel"));
        //if (temporal() != NULL) src->append(GXmlElement("temporalModel"));
    }

    // Set model type, name and optionally instruments
    src->attribute("name", name());
    src->attribute("type", type());
    if (instruments().length() > 0) {
        src->attribute("instrument", instruments());
    }
    std::string identifiers = ids();
    if (identifiers.length() > 0) {
        src->attribute("id", identifiers);
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write spatial model
    if (spatial() != NULL) {
        GXmlElement* spat = src->element("spatialModel", 0);
        spatial()->write(*spat);
    }

    // Write temporal model
    /*
    if (temporal() != NULL) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("temporalModel", 0);
            temporal()->write(*temp);
        }
    }
    */

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GCTAModelCubeBackground::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelCubeBackground ===");

        // Determine number of parameters per type
        int n_spatial  = (spatial()  != NULL) ? spatial()->size()  : 0;
        int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
        int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type"));
        if (n_spatial > 0) {
            result.append("\""+spatial()->type()+"\"");
            if (n_spectral > 0 || n_temporal > 0) {
                result.append(" * ");
            }
        }
        if (n_spectral > 0) {
            result.append("\""+spectral()->type()+"\"");
            if (n_temporal > 0) {
                result.append(" * ");
            }
        }
        if (n_temporal > 0) {
            result.append("\""+temporal()->type()+"\"");
        }

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of spatial par's") +
                      gammalib::str(n_spatial));
        for (int i = 0; i < n_spatial; ++i) {
            result.append("\n"+(*spatial())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of spectral par's") +
                      gammalib::str(n_spectral));
        for (int i = 0; i < n_spectral; ++i) {
            result.append("\n"+(*spectral())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of temporal par's") +
                      gammalib::str(n_temporal));
        for (int i = 0; i < n_temporal; ++i) {
            result.append("\n"+(*temporal())[i].print());
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Create spatial cube model from background file
 *
 * @param[in] obs CTA observation.
 * @param[in] filename Background table filename.
 * @param[in] spectral Spectral model component.
 * @param[in] nx_sky Number of bins in x direction (optional).
 * @param[in] ny_sky Number of bins in y direction (optional).
 * @param[in] n_energy Number of energy bins (optional).
 *
 * @exception GException::invalid_argument
 *            Invalid argument encountered.
 *
 * The method assumes that energies are stored in units of TeV in the
 * response table.
 *
 * @todo Document method.
 * @todo We do not really need a CTA observation, we only need the pointing
 *       direction.
 ***************************************************************************/
void GCTAModelCubeBackground::set_spatial(const GCTAObservation& obs,
                                          const std::string&     filename,
                                          const int&             nx_sky,
                                          const int&             ny_sky,
                                          const int&             n_energy)
{
    // Free any existing spatial model
    if (m_spatial  != NULL) delete m_spatial;
    m_spatial = NULL;

    // Tie model to observation by assigning the same id
    ids(obs.id());

    // Extract pointing information from CTA observation
    const GCTAPointing& pointing = obs.pointing();

    // Read the fits file with the background information
    GFits fits(filename);

    // TODO: Could have different extensions
    const GFitsTable& table = *fits.table("BACKGROUND");

    // Get the content as GCTAResponseTable
    GCTAResponseTable background(table);

    // Retrieve spatial bounds from response table
    double xlow  = background.axis_lo(0,0);
    double xhigh = background.axis_hi(0,background.axis(0)-1);
    double ylow  = background.axis_lo(1,0);
    double yhigh = background.axis_hi(1,background.axis(1)-1);

    // Read the length of the axes
    int nx;
    int ny;
    if (nx_sky > 0 && ny_sky > 0) {
        nx = nx_sky;
        ny = ny_sky;
    }
    else if (nx_sky == 0 && ny_sky == 0) {
        nx = background.axis(0);
        ny = background.axis(1);
    }
    else {
        std::string msg = "Specified number of ";
        if (nx_sky < 0) {
            msg += "x bins "+gammalib::str(nx_sky)+" ";
        }
        if (ny_sky < 0) {
            if (nx_sky < 0) {
                msg += "and ";
            }
            msg += "y bins "+gammalib::str(ny_sky)+" ";
        }
        msg += "is smaller than zero.";
        throw GException::invalid_argument(G_SET_SPATIAL, msg);
    }

    // Retrieve DETX and DETY units ("radians" if undefined)
    std::string unit_detx("rad");
    std::string unit_dety("rad");
    if (!background.axis_lo_unit(0).empty()) {
        unit_detx = background.axis_lo_unit(0);
    }
    if (!background.axis_lo_unit(1).empty()) {
        unit_dety = background.axis_lo_unit(1);
    }

    // Retrieve energy units ("TeV" if undefined)
    std::string unit_lo("TeV");
    std::string unit_hi("TeV");
    if (!background.axis_lo_unit(2).empty()) {
        unit_lo = background.axis_lo_unit(2);
    }
    if (!background.axis_hi_unit(2).empty()) {
        unit_hi = background.axis_hi_unit(2);
    }

    // Set number of energy bins
    GEnergies energies; 
    int       n_energies;
    if (n_energy > 0) {
        n_energies = n_energy;
    }
    else if (n_energy == 0) {
        n_energies = background.axis(2);
    }
    else {
        std::string msg = "Specified number of energy bins "+
                          gammalib::str(n_energy)+" is smaller than zero.";
        throw GException::invalid_argument(G_SET_SPATIAL, msg);
    }

    // Create energies as the log means of the energy bins
    GEnergy  emin(background.axis_lo(2,0), unit_lo);
    GEnergy  emax(background.axis_hi(2,background.axis(2)-1), unit_hi);
    GEbounds ebounds(n_energies, emin, emax);
    for (int i = 0; i < n_energies; ++i) {
        energies.append(ebounds.elogmean(i));
    }

    // Set interpolation unit for DETX and DETY to radians. This is only needed
    // if DETX and/or DETY are given in degrees. We also set conversion factors
    // for the x and y binsize which is needed in degrees.
    double x_convert = gammalib::rad2deg;
    double y_convert = gammalib::rad2deg;
    if (gammalib::tolower(unit_detx) == "deg") {
        background.axis_radians(0); // Switch degrees to radians
        x_convert = 1.0;
    }
    if (gammalib::tolower(unit_dety) == "deg") {
        background.axis_radians(1); // Switch degrees to radians
        y_convert = 1.0;
    }

    // Set interpolation to logscale for energies
    background.axis_log10(2);

    // Creating the sky map
    double  bin_size_x = (xhigh - xlow) / nx * x_convert;
    double  bin_size_y = (yhigh - ylow) / ny * y_convert;
    GSkymap cube       = GSkymap("TAN", "CEL",
                                 pointing.dir().ra_deg(),
                                 pointing.dir().dec_deg(),
                                 -bin_size_x,
                                 bin_size_y,
                                 nx,
                                 ny,
                                 n_energies);

    // Loop over skymap pixel
    for (int i = 0; i < cube.npix(); ++i) {

        // Get sky direction from pixel number
        GSkyDir pix_dir = cube.inx2dir(i);

        // Get instrument coordinates for this pixel
        const GCTAInstDir instdir = pointing.instdir(pix_dir);

        // loop on the energy map
        for (int j = 0; j < energies.size(); ++j) {

            // Determine background value in instrument system
            double value = background(0, instdir.detx(),
                                         instdir.dety(),
                                         energies[j].log10TeV());

            // Make sure that value is positive
            if (value < 0.0) {
                value = 0.0;
            }

            // Set skymap value
            cube(i,j) = value;

        } // endfor: loop over energies

    } // endfor: loop over sky bins

    // Create the GModelSpatialDiffuseCube
    m_spatial = new GModelSpatialDiffuseCube(cube, energies);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise class members
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelCubeBackground::init_members(void)
{
    // Initialise members
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Initialise Npred cache
    m_npred_names.clear();
    m_npred_energies.clear();
    m_npred_times.clear();
    m_npred_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelCubeBackground::copy_members(const GCTAModelCubeBackground& model)
{
    // Clone spatial, spectral and temporal model components
    m_spatial  = (model.m_spatial  != NULL) ? model.m_spatial->clone()  : NULL;
    m_spectral = (model.m_spectral != NULL) ? model.m_spectral->clone() : NULL;
    m_temporal = (model.m_temporal != NULL) ? model.m_temporal->clone() : NULL;

    // Copy cache
    m_npred_names    = model.m_npred_names;
    m_npred_energies = model.m_npred_energies;
    m_npred_times    = model.m_npred_times;
    m_npred_values   = model.m_npred_values;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelCubeBackground::free_members(void)
{
    // Free memory
    if (m_spatial  != NULL) delete m_spatial;
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModelData base class.
 ***************************************************************************/
void GCTAModelCubeBackground::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Determine the number of parameters
    int n_spatial  = (spatial()  != NULL) ? spatial()->size()  : 0;
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spatial + n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather spatial parameter pointers
        for (int i = 0; i < n_spatial; ++i) {
            m_pars.push_back(&((*spatial())[i]));
        }

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i) {
            m_pars.push_back(&((*spectral())[i]));
        }

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i) {
            m_pars.push_back(&((*temporal())[i]));
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 *
 * Returns 'true' if models has a spatial, a spectral and a temporal
 * component. Otherwise returns 'false'.
 ***************************************************************************/
bool GCTAModelCubeBackground::valid_model(void) const
{
    // Set result
    bool result = ((spatial()  != NULL) &&
                   (spectral() != NULL) &&
                   (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Construct spatial model from XML element
 *
 * @param[in] spatial XML element containing spatial model information.
 *
 * @exception GException::model_invalid_spatial
 *            Invalid spatial model type encountered.
 *
 * Returns pointer to a spatial model that is defined in an XML element.
 ***************************************************************************/
GModelSpatial* GCTAModelCubeBackground::xml_spatial(const GXmlElement& spatial) const
{
    // Get spatial model type
    std::string type = spatial.attribute("type");

    // Get spatial model
    GModelSpatialRegistry registry;
    GModelSpatial*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(spatial);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct spectral model from XML element
 *
 * @param[in] spectral XML element containing spectral model information.
 *
 * @exception GException::model_invalid_spectral
 *            Invalid spectral model type encountered.
 *
 * Returns pointer to a spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelCubeBackground::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model type
    std::string type = spectral.attribute("type");

    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(spectral);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct temporal model from XML element
 *
 * @param[in] temporal XML element containing temporal model information.
 *
 * @exception GException::model_invalid_temporal
 *            Invalid temporal model type encountered.
 *
 * Returns pointer to a temporal model that is defined in an XML element.
 ***************************************************************************/
GModelTemporal* GCTAModelCubeBackground::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model type
    std::string type = temporal.attribute("type");

    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(temporal);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_temporal(G_XML_TEMPORAL, type);
    }

    // Return pointer
    return ptr;
}



/***********************************************************************//**
 * @brief Kernel for zenith angle Npred integration of background model
 *
 * @param[in] theta Offset angle from ROI centre (radians).
 *
 * Computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \theta \times
 *                     \int_{\phi_{\rm min}}^{\phi_{\rm max}}
 *                     S_{\rm p}(\theta,\phi | E, t) \,
 *                     N_{\rm pred}(\theta,\phi) d\phi
 * \f]
 *
 * The azimuth angle integration range
 * \f$[\phi_{\rm min}, \phi_{\rm max}\f$
 * is limited to an arc around the vector connecting the model centre to
 * the ROI centre. This limitation assures proper converges properly.
 ***************************************************************************/
double GCTAModelCubeBackground::npred_roi_kern_theta::eval(const double& theta)
{
    // Initialise value
    double value = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Compute sine of offset angle
        double sin_theta = std::sin(theta);
        #if defined(G_NPRED_AROUND_ROI)
        double dphi = gammalib::pi;
        #else
        double dphi = 0.5 * gammalib::cta_roi_arclength(theta,
                                                        m_dist,
                                                        m_cosdist,
                                                        m_sindist,
                                                        m_roi,
                                                        m_cosroi);
        #endif

        // Continue only if arc length is positive
        if (dphi > 0.0) {

		   // Compute phi integration range
		   double phi_min = m_omega0 - dphi;
		   double phi_max = m_omega0 + dphi;

			// Setup phi integration kernel
            GCTAModelCubeBackground::npred_roi_kern_phi integrand(m_model,
												                  m_obsEng,
												                  m_obsTime,
												                  m_rot,
												                  theta,
												                  sin_theta);

			// Integrate over phi
			GIntegral integral(&integrand);
	        integral.eps(1e-3);
			value = integral.romberg(phi_min, phi_max) * sin_theta;

			// Debug: Check for NaN
			#if defined(G_NAN_CHECK)
			if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
				std::cout << "*** ERROR: GCTAModelCubeBackground::npred_roi_kern_theta::eval";
				std::cout << "(theta=" << theta << "):";
				std::cout << " NaN/Inf encountered";
				std::cout << " (value=" << value;
				std::cout << ", sin_theta=" << sin_theta;
				std::cout << ")" << std::endl;
			}
			#endif

        } // endif: arc length was positive

    } // endif: offset angle was positive

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for Npred azimuth angle integration of background model
 *
 * @param[in] phi Azimuth angle around ROI centre (radians).
 *
 * Computes
 *
 * \f[
 *    S_{\rm p}(\theta, \phi | E, t) \, N_{\rm pred}(\theta, \phi)
 * \f]
 *
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double GCTAModelCubeBackground::npred_roi_kern_phi::eval(const double& phi)
{
    // Initialise value
    double value = 0.0;

    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir obsDir;
    obsDir.celvector(cel);

    // Set Photon
    GPhoton photon(obsDir, m_obsEng, m_obsTime);

    // Get sky intensity for this photon
    value = m_model->eval(photon);

	// Debug: Check for NaN
	#if defined(G_NAN_CHECK)
	if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
		std::cout << "*** ERROR: GCTAModelCubeBackground::npred_roi_kern_phi::eval";
		std::cout << "(phi=" << phi << "):";
		std::cout << " NaN/Inf encountered";
		std::cout << " (value=" << value;
		std::cout << ", cos_phi=" << cos_phi;
		std::cout << ", sin_phi=" << sin_phi;
		std::cout << ")" << std::endl;
	}
	#endif

    // Return Npred
    return value;
}
