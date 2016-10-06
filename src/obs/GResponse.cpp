/***************************************************************************
 *                GResponse.cpp - Abstract response base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2016 by Juergen Knoedlseder                         *
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
 * @file GResponse.hpp
 * @brief Abstract response base class interface definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <unistd.h>           // access() function
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GIntegral.hpp"
#include "GResponse.hpp"
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GSource.hpp"        // will become obsolete
#include "GEbounds.hpp"       // will become obsolete
#include "GObservation.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialComposite.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF_RADIAL               "GResponse::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_ELLIPTICAL       "GResponse::irf_elliptical(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE             "GResponse::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"

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
GResponse::GResponse(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp Response.
 ***************************************************************************/
GResponse::GResponse(const GResponse& rsp)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GResponse::~GResponse(void)
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
 * @param[in] rsp Response.
 * @return Response.
 ***************************************************************************/
GResponse& GResponse::operator=(const GResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(rsp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Convolve sky model with the instrument response
 *
 * @param[in] model Sky model.
 * @param[in] event Event.
 * @param[in] obs Observation.
 * @param[in] grad Should model gradients be computed? (default: true)
 * @return Event probability.
 *
 * Computes the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * without taking into account any time dispersion. Energy dispersion is
 * correctly handled by this method. If time dispersion is indeed needed,
 * an instrument specific method needs to be provided.
 ***************************************************************************/
double GResponse::convolve(const GModelSky&    model,
                           const GEvent&       event,
                           const GObservation& obs,
                           const bool&         grad) const
{
    // Set number of iterations for Romberg integration.
    static const int iter = 6;

    // Initialise result
    double prob = 0.0;

    // Continue only if the model has a spatial component
    if (model.spatial() != NULL) {

        // Get source time (no dispersion)
        GTime srcTime = event.time();

        // Case A: Integration
        if (use_edisp()) {
    
            // Retrieve true energy boundaries
            GEbounds ebounds = this->ebounds(event.energy());
    
            // Loop over all boundaries
            for (int i = 0; i < ebounds.size(); ++i) {

                // Get boundaries in MeV
                double emin = ebounds.emin(i).MeV();
                double emax = ebounds.emax(i).MeV();

                // Continue only if valid
                if (emax > emin) {

                    // Setup integration function
                    edisp_kern integrand(this, &obs, &model, &event, srcTime, grad);
                    GIntegral  integral(&integrand);

                    // Set number of iterations
                    integral.fixed_iter(iter);

                    // Do Romberg integration
                    emin  = std::log(emin);
                    emax  = std::log(emax);
                    prob += integral.romberg(emin, emax);
    
                } // endif: interval was valid

            } // endfor: looped over intervals

        }

        // Case B: No integration (assume no energy dispersion)
        else {
    
            // Get source energy (no dispersion)
            GEnergy srcEng  = event.energy();

            // Evaluate probability
            prob = eval_prob(model, event, srcEng, srcTime, obs, grad);

        }

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(prob) || gammalib::is_infinite(prob)) {
            std::cout << "*** ERROR: GResponse::convolve:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (prob=" << prob;
            std::cout << ", event=" << event;
            std::cout << ", srcTime=" << srcTime;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: spatial component valid

    // Return probability
    return prob;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GResponse::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response.
 ***************************************************************************/
void GResponse::copy_members(const GResponse& rsp)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GResponse::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Convolve sky model with the instrument response
 *
 * @param[in] model Sky model.
 * @param[in] event Event.
 * @param[in] srcEng Source energy.
 * @param[in] srcTime Source time.
 * @param[in] obs Observation.
 * @param[in] grad Should model gradients be computed? (default: true)
 * @return Event probability.
 *
 * Computes the event probability
 *
 * \f[
 *    P(p',E',t'|E,t) = \int S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * for a given true energy \f$E\f$ and time \f$t\f$.
 ***************************************************************************/
double GResponse::eval_prob(const GModelSky&    model,
                            const GEvent&       event,
                            const GEnergy&      srcEng,
                            const GTime&        srcTime,
                            const GObservation& obs,
                            const bool&         grad) const
{
    // Initialise result
    double prob = 0.0;

    // Continue only if the model has a spatial component
    if (model.spatial() != NULL) {

        // Initialise IRF value
        double irf = 0.0;

        // Try to retrieve a composite spectral model
        const GModelSpatialComposite* composite = dynamic_cast<const
                GModelSpatialComposite*>(model.spatial());

        // Compute IRF using GSource if not a composite model
        if (composite == NULL) {

            // Set source
            GSource source(model.name(), model.spatial(), srcEng, srcTime);

            // Get IRF value. This method returns the spatial component of the
            // source model.
            irf = this->irf(event, source, obs);
        }
        else {

            // Loop over composite model components
            for(int i = 0; i < composite->components(); ++i) {

                // Circumvent function returning const component
                GModelSpatial* spat = const_cast<GModelSpatial*>(composite->component(i));

                // Set source
                GSource source(model.name(), spat, srcEng, srcTime);

                // Get IRF value. This method returns the spatial component of the
                // source model.
                irf += this->irf(event, source, obs);
            }
        }

        // If required, apply instrument specific model scaling
        if (model.has_scales()) {
            irf *= model.scale(obs.instrument()).value();
        }

        // Evaluate spectral and temporal components
        double spec = (model.spectral() != NULL)
                      ? model.spectral()->eval(srcEng, srcTime, grad) : 1.0;
        double temp = (model.temporal() != NULL)
                      ? model.temporal()->eval(srcTime, grad) : 1.0;

        // Compute probability
        prob = spec * temp * irf;

        // Optionally compute partial derivatives
        if (grad) {

            // Multiply factors to spectral gradients
            if (model.spectral() != NULL) {
                double fact = temp * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < model.spectral()->size(); ++i) {
                        (*model.spectral())[i].factor_gradient((*model.spectral())[i].factor_gradient() * fact);
                    }
                }
            }

            // Multiply factors to temporal gradients
            if (model.temporal() != NULL) {
                double fact = spec * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < model.temporal()->size(); ++i) {
                        (*model.temporal())[i].factor_gradient((*model.temporal())[i].factor_gradient() * fact);
                    }
                }
            }

        } // endif: computed partial derivatives

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(prob) || gammalib::is_infinite(prob)) {
            std::cout << "*** ERROR: GResponse::eval_prob:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (prob=" << prob;
            std::cout << ", spec=" << spec;
            std::cout << ", temp=" << temp;
            std::cout << ", irf=" << irf;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: Gamma-ray source model had a spatial component

    // Return event probability
    return prob;
}


/***********************************************************************//**
 * @brief Integration kernel for edisp_kern() method
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel needed for the edisp_kern()
 * method.
 ***************************************************************************/
double GResponse::edisp_kern::eval(const double& x)
{
    // Set energy
    GEnergy eng;
    double expx = std::exp(x);
    eng.MeV(expx);

    // Get function value
    double value = m_parent->eval_prob(*m_model, *m_event, eng, m_srcTime, *m_obs, m_grad);

    // Save value if needed
    #if defined(G_NAN_CHECK)
    double value_out = value;
    #endif

    // Correct for variable substitution
    value *= expx;

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GResponse::edisp_kern::eval";
        std::cout << "(x=" << x << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value_out;
        std::cout << " exp(x)=" << expx;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}
