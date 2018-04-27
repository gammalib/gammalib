/***************************************************************************
 *   GModelSpatialRadialProfile.cpp - Radial profile source model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialRadialProfile.cpp
 * @brief Radial profile model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"
#include "GModelSpatialRadialProfile.hpp"
#include <iomanip>

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_PRECOMPUTATION


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty radial profile
 ***************************************************************************/
GModelSpatialRadialProfile::GModelSpatialRadialProfile(void) :
                            GModelSpatialRadial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs radial profile model by extracting information from an XML
 * elements.
 ***************************************************************************/
GModelSpatialRadialProfile::GModelSpatialRadialProfile(const GXmlElement& xml) :
                            GModelSpatialRadial()
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
 * @param[in] model Radial profile model.
 *
 * Copies radial profile model from another radial profile model.
 ***************************************************************************/
GModelSpatialRadialProfile::GModelSpatialRadialProfile(const GModelSpatialRadialProfile& model) :
                            GModelSpatialRadial(model)
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
 *
 * Destructs radial profile model.
 ***************************************************************************/
GModelSpatialRadialProfile::~GModelSpatialRadialProfile(void)
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
 * @param[in] model Radial profile model.
 * @return Radial profile model.
 *
 * Assigns radial profile model.
 ***************************************************************************/
GModelSpatialRadialProfile& GModelSpatialRadialProfile::operator=(const GModelSpatialRadialProfile& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialRadial::operator=(model);

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
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] theta Angular distance from model centre (radians).
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * Evaluate the radial profile model for a given angular distance @p theta
 * the model centre, a given photon @p energy and a given @p time. The
 * method evaluates the model by linear interpolation.
 ***************************************************************************/
double GModelSpatialRadialProfile::eval(const double&  theta,
                                        const GEnergy& energy,
                                        const GTime&   time,
                                        const bool&    gradients) const
{
    // Get pre-computation cache index
    int icache = cache_index();
    
    // Get interpolation value
    double value = m_profile[icache].nodes.interpolate(theta, m_profile[icache].values);

    // Make sure that value is not-negative
    if (value < 0.0) {
        value = 0.0;
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpatialRadialProfile::eval";
        std::cout << "(theta=" << theta << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction from the radial profile using a rejection
 * method.
 ***************************************************************************/
GSkyDir GModelSpatialRadialProfile::mc(const GEnergy& energy,
                                       const GTime&   time,
                                       GRan&          ran) const
{
    // Get pre-computation cache index
    int icache = cache_index();

    // Simulate theta angle (degrees) using a rejection method
    double theta = 0.0;
    double theta_width = theta_max() - theta_min() ;
    while (true) {
        theta     = theta_min() + ( ran.uniform() * theta_width ) ;
        double mc = m_profile[icache].nodes.interpolate(theta, m_profile[icache].mc);
        if ((ran.uniform() * m_profile[icache].mc_max) <= mc) {
            break;
        }
    }
    theta *= gammalib::rad2deg;

    // Simulate uniform azmiuth angle (degrees)
    double phi = 360.0 * ran.uniform();

    // Rotate sky direction by offset
    GSkyDir sky_dir = dir();
    sky_dir.rotate_deg(phi, theta);

    // Return sky direction
    return sky_dir;
}


/***********************************************************************//**
 * @brief Checks where model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the radial disk model.
 ***************************************************************************/
bool GModelSpatialRadialProfile::contains(const GSkyDir& dir,
                                          const double&  margin) const
{
    // Compute distance to centre (radians)
    double distance = dir.dist(this->dir());

    // Return flag
    return ( ( distance <= theta_max() + margin*gammalib::deg2rad ) && 
             ( distance >= theta_min() + margin*gammalib::deg2rad ) ) ;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialRadialProfile::init_members(void)
{
    // Initialise members
    m_coord_indep = false;
    m_num_nodes   = 100;
    m_region.clear();

    // Initialise pre-computation cache
    m_profile.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial disk model.
 ***************************************************************************/
void GModelSpatialRadialProfile::copy_members(const GModelSpatialRadialProfile& model)
{
    // Copy members
    m_coord_indep = model.m_coord_indep;
    m_num_nodes   = model.m_num_nodes;
    m_profile     = model.m_profile;
    m_region      = model.m_region;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialProfile::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return index to pre-computation cache
 *
 * @return Index to pre-computation cache
 *
 * Returns the index to the pre-computation cache. If no pre-computation
 * cache was found the method will create one and return the index to that
 * cache.
 ***************************************************************************/
int GModelSpatialRadialProfile::cache_index(void) const
{
    // Initialise index
    int index = -1;

    // Search for index
    for (int i = 0; i < m_profile.size(); ++i) {

        // Initialise found flag
        bool found = true;

        // Loop over all spatial parameters and check if they are equal
        // to the parameters in the pre-computation cache. Break and set
        // the found flag to false on non-equality.
        for (int k = 0; k < m_pars.size(); ++k) {
            
            // Skip if model is coordinate independent and par is RA or DEC
            if (m_coord_indep && (m_pars[k]->name()=="RA" || m_pars[k]->name()=="DEC")) {
                continue;
            } 
            // Otherwise ...
            else if (m_pars[k]->value() != m_profile[i].pars[k]) {
                found = false;
                break;
            }
        }

        // If profile was found then set index
        if (found) {
            index = i;
            break;
        }

    } // endfor: looped over index

    // If index was not found then allocate a new pre-computation cache
    if (index == -1) {

        // Initialise profile
        profile prf;
        prf.nodes.clear();
        prf.values.clear();
        prf.mc.clear();
        prf.mc_max = 0.0;

        // Set profile parameters
        std::cout << "CACHE:" << std::endl;
        for (int k = 0; k < m_pars.size(); ++k) {
            prf.pars.push_back(m_pars[k]->value());
            std::cout << "  pars[" << m_pars[k]->name() << "] = " << m_pars[k]->value() << std::endl;
        }

        // Pre-compute radial values, MC values, and mc_max. Compute
        // also the normalization.
        double rmax = theta_max();
        double dr   = rmax / m_num_nodes;
        std::cout << "  rmax=" << rmax << "  dr=" << dr << std::endl;
        double r    = theta_min();
        double norm = 0.0;
        for (int j = 0; j < m_num_nodes; ++j) {
            double value = profile_value(r);
            double mc    = value * std::sin(r) * dr;
            std::cout << "  j=" << j << "  r=" << r << "  (deg=" ;
            std::cout.width(4) ;
            std::cout << r * gammalib::rad2deg << ")  value=" ;
            std::cout.width(7);
            std::cout << value << "  mc=" ;
            std::cout.width(7) ;
            std::cout << mc << std::endl;
            norm        += mc;
            if (mc > prf.mc_max) {
                prf.mc_max = mc;
            }
            prf.nodes.append(r);
            prf.values.push_back(value);
            prf.mc.push_back(mc);
            r += dr;
        }
        norm *= gammalib::twopi;
        std::cout << "  prf.mc_max=" << prf.mc_max << "  norm=" << norm << std::endl;

        // Normalize radial profile
        if (norm > 0.0) {
            double inv_norm = 1.0 / norm;
            for (int j = 0; j < m_num_nodes; ++j) {
                prf.values[j] *= inv_norm;
            }
        }

        // Push profile in pre-computation cache
        m_profile.push_back(prf);

        // Get index to last profile
        index = m_profile.size() - 1;

        // Log pre-computation
        #if defined(G_DEBUG_PRECOMPUTATION)
        std::cout << "GModelSpatialRadialProfile::cache_index";
        std::cout << std::endl;
        std::cout << "  dr=" << dr << std::endl;
        for (int k = 0; k < m_pars.size(); ++k) {
            std::cout << "  par[" << k << "]=" << m_pars[k]->value();
            std::cout << " (" << m_pars[k]->name();
            std::cout << ")" << std::endl;
        }
        std::cout << "  norm=" << norm << std::endl;
        #endif

    } // endif: allocated new pre-computation cache

    // Return index
    return (index);
}


/***********************************************************************//**
 * @brief Set boundary sky region
 ***************************************************************************/
void GModelSpatialRadialProfile::set_region(void) const
{
    // Set sky region centre to disk centre
    m_region.centre(m_ra.value(), m_dec.value());

    // Set sky region radius to maximum theta angle
    m_region.radius(theta_max()*gammalib::rad2deg);

    // Return
    return;
}
