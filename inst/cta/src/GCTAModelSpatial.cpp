/***************************************************************************
 *         GCTAModelSpatial.cpp - Spatial model abstract base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Jurgen Knodlseder                                *
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
 * @file GCTAModelSpatial.cpp
 * @brief Abstract spatial model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GObservation.hpp"
#include "GCTASupport.hpp"
#include "GCTAEventList.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAModelSpatial.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                        "GCTAModelSpatial::operator[](int&)"
#define G_ACCESS2                "GCTAModelSpatial::operator[](std::string&)"
#define G_NPRED    "GCTAModelSpatial::npred(GEnergy&, GTime&, GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_NPRED_CACHE

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_NPRED                           //!< Dump Npred information


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelSpatial::GCTAModelSpatial(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spatial model.
 ***************************************************************************/
GCTAModelSpatial::GCTAModelSpatial(const GCTAModelSpatial& model)
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
GCTAModelSpatial::~GCTAModelSpatial(void)
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
 * @param[in] model Spatial model.
 ***************************************************************************/
GCTAModelSpatial& GCTAModelSpatial::operator=(const GCTAModelSpatial& model)
{
    // Execute only if object is not identical
    if (this != &model) {

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


/***********************************************************************//**
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar& GCTAModelSpatial::operator[](const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ACCESS1, "Spatial parameter index",
                                       index, size());
    }
    #endif

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
const GModelPar& GCTAModelSpatial::operator[](const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ACCESS1, "Spatial parameter index",
                                       index, size());
    }
    #endif

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found in container.
 ***************************************************************************/
GModelPar& GCTAModelSpatial::operator[](const std::string& name)
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name)
            break;
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Parameter \""+name+"\" not found in spatial "
                          "component of background model.";
        throw GException::invalid_argument(G_ACCESS2, msg);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter (const version)
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found in container.
 ***************************************************************************/
const GModelPar& GCTAModelSpatial::operator[](const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name)
            break;
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Parameter \""+name+"\" not found in spatial "
                          "component of background model.";
        throw GException::invalid_argument(G_ACCESS2, msg);
    }

    // Return reference
    return *(m_pars[index]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return integral of spatial model component
 *
 * @param[in] energy Measured event energy.
 * @param[in] time Measured event time.
 * @param[in] obs Observation.
 * @return Integral of spatial model component.
 *
 * Spatially integrates the spatial background model component for a given
 * measured event energy and event time over the region of interest (RoI)
 * of a given observation.
 *
 * The method uses a 2D Romberg integration to numerically integrate the
 * spatial background model component.
 *
 * The method implements an integration cache so that for a given
 * observation, energy and time the values are cached and not recomputed.
 *
 * @todo The method currently assumes that the RoI is centred on DETX=DETY=0
 ***************************************************************************/
double GCTAModelSpatial::npred(const GEnergy&      energy,
                               const GTime&        time,
                               const GObservation& obs) const
{
    // Set number of iterations for Romberg integration.
    static const int iter_theta = 6;
    static const int iter_phi   = 6;

    // Initialise result
    double npred = 0.0;

    // Initialise flag
    bool has_npred = false;

    // Build unique identifier
    std::string id = obs.instrument() + "::" + obs.id();

    // Check if Npred value is already in cache
    #if defined(G_USE_NPRED_CACHE)
    if (!m_npred_names.empty()) {

        // Search for unique identifier, and if found, recover Npred
        // value and break
        for (int i = 0; i < m_npred_names.size(); ++i) {
            if ((m_npred_names[i]    == id) &&
                (m_npred_energies[i] == energy) &&
                (m_npred_times[i]    == time)) {
                npred     = m_npred_values[i];
                has_npred = true;
                #if defined(G_DEBUG_NPRED)
                std::cout << "GCTAModelSpatial::npred:";
                std::cout << " cache=" << i;
                std::cout << " npred=" << npred << std::endl;
                #endif
                break;
            }
        }

    } // endif: there were values in the Npred cache
    #endif

    // Continue only if no Npred cache value has been found
    if (!has_npred) {

        // Retrieve CTA event list
        const GCTAEventList& events = gammalib::cta_event_list(G_NPRED, obs);

        // Get reference to RoI centre
        //const GSkyDir& roi_centre = events.roi().centre().dir();

        // Get RoI radius in radians
        double roi_radius = events.roi().radius() * gammalib::deg2rad;

        // Setup integration function
        GCTAModelSpatial::npred_roi_kern_theta integrand(this,
                                                         energy,
                                                         time,
                                                         iter_phi);

        // Setup integration
        GIntegral integral(&integrand);

        // Set fixed number of iterations
        integral.fixed_iter(iter_theta);

        // Spatially integrate radial component (assumes that RoI centre is
        // at DETX=DETY=0)
        npred = integral.romberg(0.0, roi_radius);

        // Store result in Npred cache
        #if defined(G_USE_NPRED_CACHE)
        m_npred_names.push_back(id);
        m_npred_energies.push_back(energy);
        m_npred_times.push_back(time);
        m_npred_values.push_back(npred);
        #endif

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
            std::string origin  = "GCTAModelSpatial::npred";
            std::string message = " NaN/Inf encountered (npred=" +
                                  gammalib::str(npred) + ", roi_radius=" +
                                  gammalib::str(roi_radius) + ")";
            gammalib::warning(origin, message);
        }
        #endif

    } // endif: Npred computation required

    // Return Npred
    return npred;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelSpatial::init_members(void)
{
    // Initialise members
    m_pars.clear();

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
 * @param[in] model Spatial background model component.
 ***************************************************************************/
void GCTAModelSpatial::copy_members(const GCTAModelSpatial& model)
{
    // Copy members
    m_pars = model.m_pars;

    // Copy cache
    m_npred_names    = model.m_npred_names;
    m_npred_energies = model.m_npred_energies;
    m_npred_times    = model.m_npred_times;
    m_npred_values   = model.m_npred_values;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelSpatial::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Kernel for offset angle integration of spatial component
 *
 * @param[in] theta Offset angle from RoI centre (radians).
 *
 * Computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \theta \times
 *                     \int_{0}^{2\pi}
 *                     B(\theta,\phi | E, t) d\phi
 * \f]
 *
 * where \f$B(\theta,\phi | E, t)\f$ is the spatial component of the
 * background model for a specific observed energy \f$E\f$ and time \f$t\f$.
 ***************************************************************************/
double GCTAModelSpatial::npred_roi_kern_theta::eval(const double& theta)
{
    // Initialise value
    double value = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Setup phi integration kernel
        GCTAModelSpatial::npred_roi_kern_phi integrand(m_spatial,
                                                       m_energy,
                                                       m_time,
                                                       theta);

        // Setup integration
        GIntegral integral(&integrand);

        // Set fixed number of iterations
        integral.fixed_iter(m_iter);

        // Integrate over phi
        value = integral.romberg(0.0, gammalib::twopi) * std::sin(theta);

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::string origin  = "GCTAModelSpatial::npred_roi_kern_theta::eval"
                                  "(" + gammalib::str(theta) + ")";
            std::string message = " NaN/Inf encountered (value=" +
                                  gammalib::str(value) + ")";
            gammalib::warning(origin, message);
        }
        #endif

    } // endif: offset angle was positive

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle integration of spatial component
 *
 * @param[in] phi Azimuth angle around RoI centre (radians).
 *
 * Computes
 *
 * \f[
 *    B(\theta, \phi | E, t)
 * \f]
 *
 * using
 *
 * \f[ {\rm detx} = \theta \cos \phi \f]
 * \f[ {\rm dety} = \theta \sin \phi \f]
 *
 * @todo Verify correct orientation of detx and dety with respect to phi
 ***************************************************************************/
double GCTAModelSpatial::npred_roi_kern_phi::eval(const double& phi)
{
    // Compute detx and dety in radians
    double detx(0.0);
    double dety(0.0);
    if (m_theta > 0.0 ) {
        detx = m_theta * std::cos(phi);
        dety = m_theta * std::sin(phi);
    }

    // Setup CTA instrument direction
    GCTAInstDir dir(detx, dety);

    // Get background value
    double value = m_spatial->eval(dir, m_energy, m_time);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::string origin  = "GCTAModelSpatial::npred_roi_kern_phi::eval"
                              "(" + gammalib::str(phi) + ")";
        std::string message = " NaN/Inf encountered (value=" +
                              gammalib::str(value) + ", detx=" +
                              gammalib::str(detx) + ", dety=" +
                              gammalib::str(dety) + ")";
        gammalib::warning(origin, message);
    }
    #endif

    // Return Npred
    return value;
}
