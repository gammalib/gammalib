/***************************************************************************
 *          GModelSpatial.cpp - Abstract spatial model base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2023 by Juergen Knoedlseder                         *
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
 * @file GModelSpatial.cpp
 * @brief Abstract spatial model base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialPointSource.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                    "GModelSpatial::operator[](std::string&)"
#define G_AT                             "GModelPar& GModelSpatial::at(int&)"
#define G_FLUX     "GModelSpatial::flux(GSkyRegionCircle&, GEnergy&, GTime&)"

/* __ Constants __________________________________________________________ */
const double g_kludge_radius = 1.0e-12;            //!< Tiny angle (radians)

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
GModelSpatial::GModelSpatial(void)
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
GModelSpatial::GModelSpatial(const GModelSpatial& model)
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
GModelSpatial::~GModelSpatial(void)
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
 * @return Spatial model.
 ***************************************************************************/
GModelSpatial& GModelSpatial::operator=(const GModelSpatial& model)
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
 * @param[in] name Parameter name.
 * @return Model parameter reference.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found.
 *
 * Returns reference to the model parameter of specified @p name.
 ***************************************************************************/
GModelPar& GModelSpatial::operator[](const std::string& name)
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Model parameter \""+name+"\" not found in model. "
                          "Please specify a valid model parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] name Parameter name.
 * @return Model parameter reference.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found.
 *
 * Returns reference to the model parameter of specified @p name.
 ***************************************************************************/
const GModelPar& GModelSpatial::operator[](const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Model parameter \""+name+"\" not found in model. "
                          "Please specify a valid model parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
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
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()[.
 * @return Model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Returns model parameter with @p index range checking.
 ***************************************************************************/
GModelPar& GModelSpatial::at(const int& index)
{
    // Compile option: raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()[.
 * @return Model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Returns model parameter with @p index range checking.
 ***************************************************************************/
const GModelPar& GModelSpatial::at(const int& index) const
{
    // Compile option: raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Checks if parameter name exists
 *
 * @param[in] name Parameter name.
 * @return True if parameter with specified @p name exists.
 *
 * Searches all parameter names for a match with the specified @p name. If
 * the specified name has been found, true is returned.
 ***************************************************************************/
bool GModelSpatial::has_par(const std::string& name) const
{
    // Default found flag to false
    bool found = false;

    // Search for parameter name
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i]->name() == name) {
            found = true;
            break;
        }
    }

    // Return
    return found;
}


/***********************************************************************//**
 * @brief Checks if the spatial model has free parameters
 *
 * @return True if spatial model has free parameters.
 ***************************************************************************/
bool GModelSpatial::has_free_pars(void) const
{
    // Initialise result flag
    bool has_free_pars = false;

    // Search for free parameters
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i]->is_free()) {
            has_free_pars = true;
            break;
        }
    }

    // Return result
    return has_free_pars;
}


/***********************************************************************//**
 * @brief Autoscale parameters
 *
 * Sets the scale factors for all parameters so that the values are unity.
 ***************************************************************************/
void GModelSpatial::autoscale(void)
{
    // Loop over all parameters
    for (int i = 0; i < m_pars.size(); ++i) {
        if (m_pars[i] != NULL) {
            m_pars[i]->autoscale();
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns model flux within sky region
 *
 * @param[in] region Sky region.
 * @param[in] srcEng Energy.
 * @param[in] srcTime Time.
 * @return Flux (adimensional or ph/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} \sin \rho \times
 *    \int_{\omega} M(\rho, \omega | E, t) d\omega d\rho
 * \f]
 *
 * where
 * \f$M(\rho, \omega | E, t)\f$ is the spatial model,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the azimuth angle around the model centre.
 *
 * The integration is performed over the circular region that bounds the
 * spatial model.
 ***************************************************************************/
double GModelSpatial::flux(const GSkyRegion& region,
                           const GEnergy&    srcEng,
                           const GTime&      srcTime) const
{
    // Initialise flux
    double flux = 0.0;

    // Make sure that m_region is set
    set_region();

    // Continue only if region overlaps with model
    if (m_region.overlaps(region)) {

        // Retrieve model centre and radial integration boundaries in radians
        GSkyDir centre  = m_region.centre();
        double  rho_min = 0.0;
        double  rho_max = m_region.radius() * gammalib::deg2rad;

        // Setup integration kernel
        GModelSpatial::circle_int_kern_rho integrand(this,
                                                     region,
                                                     centre,
                                                     srcEng,
                                                     srcTime);
        GIntegral integral(&integrand);

        // Suppress integration warnings
        integral.silent(true);

        // Perform integration
        flux = integral.romberg(rho_min, rho_max);

    } // endif: model overlaps with region

    // Return
    return flux;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatial::init_members(void)
{
    // Initialise members
    m_type.clear();
    m_region.clear();
    m_pars.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spatial model.
 ***************************************************************************/
void GModelSpatial::copy_members(const GModelSpatial& model)
{
    // Copy members
    m_type   = model.m_type;
    m_region = model.m_region;
    m_pars   = model.m_pars;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatial::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Kernel for spatial model radial integration
 *
 * @param[in] rho Radial distance from model centre (radians).
 * @return Integration kernel defined as
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho | E, t) d\rho
 * \f]
 *
 * of a spatial model over a circular region. The eval() method computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     \int_{\omega} M(\rho, \omega | E, t) d\omega
 * \f]
 *
 * where
 * \f$M(\rho, \omega | E, t)\f$ is the spatial model,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the azimuth angle around the model centre.
 ***************************************************************************/
double GModelSpatial::circle_int_kern_rho::eval(const double& rho)
{
    // Initialise model value
    double flux = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Set omega integration range
        const double omega_min = 0.0;
        const double omega_max = gammalib::twopi;

        // Setup integration kernel for azimuth integration
        GModelSpatial::circle_int_kern_omega integrand(m_model,
                                                       m_region,
                                                       m_centre,
                                                       rho,
                                                       m_srcEng,
                                                       m_srcTime);

        // Setup integrator
        GIntegral integral(&integrand);

        // Compute sine term for radial integration
        double sin_rho = std::sin(rho);

        // Integrate over omega
        flux = integral.romberg(omega_min, omega_max) * sin_rho;

    } // endif: rho was positive

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Kernel for spatial model azimuth angle integration
 *
 * @param[in] omega Azimuth angle (radians).
 * @return Integration kernel defined as
 *
 * \f[
 *    K(\omega | \rho, E, t) = M(\omega | \rho)
 * \f]
 *
 * where
 * \f$M(\omega | \rho)\f$ is the spatial model,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the azimuth angle around the model centre.
 ***************************************************************************/
double GModelSpatial::circle_int_kern_omega::eval(const double& omega)
{
    // Initialise flux
    double flux = 0.0;

    // Compute sky direction corresponding to given rho, omega by rotating
    // model centre
    GSkyDir dir = m_centre;
    dir.rotate(omega, m_rho);

    // Evaluate model if sky direction is contained in region
    if (m_region.contains(dir)) {

        // Set photon for this sky direction
        GPhoton photon = GPhoton(dir, m_srcEng, m_srcTime);

        // Evaluate model for this sky direction
        flux = m_model->eval(photon);

    } // endif: sky direction was contained in region

    // Return
    return flux;
}
