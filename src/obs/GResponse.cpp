/***************************************************************************
 *                GResponse.cpp - Abstract response base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
#include "GResponse.hpp"
#include "GObservation.hpp"
#include "GIntegral.hpp"
#include "GVector.hpp"
#include "GSkyDir.hpp"
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialDiffuse.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF_RADIAL               "GResponse::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_ELLIPTICAL       "GResponse::irf_elliptical(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE             "GResponse::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_NPRED_RADIAL     "GResponse::npred_radial(GSource&, GObservation&)"
#define G_NPRED_ELLIPTICAL            "GResponse::npred_elliptical(GSource&,"\
                                                            " GObservation&)"
#define G_NPRED_DIFFUSE   "GResponse::npred_diffuse(GSource&, GObservation&)"
#define G_EBOUNDS_SRC                      "GResponse::ebounds_src(GEnergy&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_NPRED_RADIAL                        //!< Debug npred_radial


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
GResponse& GResponse::operator= (const GResponse& rsp)
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
 * @brief Return value of instrument response function
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * Returns the instrument response function for a given event, source and
 * observation.
 ***************************************************************************/
double GResponse::irf(const GEvent&       event,
                      const GSource&      source,
                      const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Select IRF depending on the spatial model type
    switch (source.model()->code()) {
        case GMODEL_SPATIAL_POINT_SOURCE:
            irf = irf_ptsrc(event, source, obs);
            break;
        case GMODEL_SPATIAL_RADIAL:
            irf = irf_radial(event, source, obs);
            break;
        case GMODEL_SPATIAL_ELLIPTICAL:
            irf = irf_elliptical(event, source, obs);
            break;
        case GMODEL_SPATIAL_DIFFUSE:
            irf = irf_diffuse(event, source, obs);
            break;
        default:
            break;
    }

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of point source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Value of instrument response function for a point source.
 *
 * This method returns the value of the instrument response function for a
 * point source. The method assumes that source.model() is of type
 * GModelSpatialPointSource.
 ***************************************************************************/
double GResponse::irf_ptsrc(const GEvent&       event,
                            const GSource&      source,
                            const GObservation& obs) const
{
    // Get point source spatial model
    const GModelSpatialPointSource* src =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Set Photon
    GPhoton photon(src->dir(), source.energy(), source.time());
    
    // Compute IRF
    double irf = this->irf(event, photon, obs);

    // Return IRF
    return irf;
}


/***********************************************************************//**
 * @brief Return IRF value for radial source model
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_radial(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_RADIAL,
          "IRF computation not implemented for radial models.");

    // Return IRF
    return 0.0;
}


/***********************************************************************//**
 * @brief Return IRF value for elliptical source model
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_elliptical(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_ELLIPTICAL,
          "IRF computation not implemented for elliptical models.");

    // Return IRF
    return 0.0;
}


/***********************************************************************//**
 * @brief Return value of diffuse source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_diffuse(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_DIFFUSE,
          "IRF computation not implemented for diffuse models.");

    // Return IRF
    return 0.0;
}


/***********************************************************************//**
 * @brief Return data space integral of instrument response function
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * Returns the data space integral of the instrument response function for
 * a given source and a particular observation. This method is needed for
 * an unbinned maximum likelihood analysis.
 *
 * The method applies the deadtime correction, so that the result can be
 * directly multiplied by the exposure time (also known as ontime).
 ***************************************************************************/
double GResponse::npred(const GSource& source, const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Select NPRED depending on the spatial model type
    switch (source.model()->code()) {
        case GMODEL_SPATIAL_POINT_SOURCE:
            npred = npred_ptsrc(source, obs);
            break;
        case GMODEL_SPATIAL_RADIAL:
            npred = npred_radial(source, obs);
            break;
        case GMODEL_SPATIAL_ELLIPTICAL:
            npred = npred_elliptical(source, obs);
            break;
        case GMODEL_SPATIAL_DIFFUSE:
            npred = npred_diffuse(source, obs);
            break;
        default:
            break;
    }

    // Return response value
    return npred;
}


/***********************************************************************//**
 * @brief Return ROI integral of point source model
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Integral of point source model over ROI.
 *
 * This method returns the spatial integral of a point source model over the
 * region of interest. The method assumes that source.model() is of type
 * GModelSpatialPointSource.
 ***************************************************************************/
double GResponse::npred_ptsrc(const GSource& source,
                              const GObservation& obs) const
{
    // Get point source spatial model
    const GModelSpatialPointSource* src =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Set Photon
    GPhoton photon(src->dir(), source.energy(), source.time());

    // Compute Npred
    double npred = this->npred(photon, obs);

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return spatial integral of radial source model
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Integral of radial source model over ROI.
 *
 * This method returns the spatial integral of a radial source model over
 * the region of interest. The method assumes that source.model() is of type
 * GModelSpatialRadial.
 ***************************************************************************/
double GResponse::npred_radial(const GSource& source,
                               const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Get radial spatial model
    const GModelSpatialRadial* spatial =
          static_cast<const GModelSpatialRadial*>(source.model());

    // Compute rotation matrix to convert from native coordinates given
    // by (theta,phi) into celestial coordinates.
    GMatrix ry;
    GMatrix rz;
    ry.eulery(spatial->dec() - 90.0);
    rz.eulerz(-spatial->ra());
    GMatrix rot = (ry * rz).transpose();

    // Set offset angle integration range
    double theta_min = 0.0;
    double theta_max = spatial->theta_max();

    // Perform offset angle integration if interval is valid
    if (theta_max > theta_min) {

        // Setup integration kernel
        GResponse::npred_radial_kern_theta integrand(*this,
                                                     *spatial,
                                                     source.energy(),
                                                     source.time(),
                                                     obs,
                                                     rot);

        // Integrate over theta
        GIntegral integral(&integrand);
        npred = integral.romberg(theta_min, theta_max);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NPRED_RADIAL)
        std::cout << "GResponse::npred_radial:";
        std::cout << " theta_min=" << theta_min;
        std::cout << " theta_max=" << theta_max;
        std::cout << " npred=" << npred << std::endl;
        #endif

    } // endif: offset angle range was valid

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_radial:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", theta_min=" << theta_min;
        std::cout << ", theta_max=" << theta_max;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return spatial integral of elliptical source model
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * This method returns the spatial integral of a radial source model over
 * the region of interest. The method assumes that source.model() is of type
 * GModelSpatialElliptical.
 ***************************************************************************/
double GResponse::npred_elliptical(const GSource& source,
                                   const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Get elliptical spatial model
    const GModelSpatialElliptical* spatial =
          static_cast<const GModelSpatialElliptical*>(source.model());

    // Compute rotation matrix to convert from native coordinates given
    // by (theta,phi) into celestial coordinates.
    GMatrix ry;
    GMatrix rz;
    ry.eulery(spatial->dec() - 90.0);
    rz.eulerz(-spatial->ra());
    GMatrix rot = (ry * rz).transpose();

    // Set offset angle integration range
    double theta_min = 0.0;
    double theta_max = spatial->theta_max();

    // Perform offset angle integration if interval is valid
    if (theta_max > theta_min) {

        // Setup integration kernel
        GResponse::npred_elliptical_kern_theta integrand(*this,
                                                         *spatial,
                                                         source.energy(),
                                                         source.time(),
                                                         obs,
                                                         rot);

        // Integrate over theta
        GIntegral integral(&integrand);
        npred = integral.romberg(theta_min, theta_max);

        // Compile option: Show integration results
        #if defined(G_DEBUG_NPRED_ELLIPTICAL)
        std::cout << "GResponse::npred_elliptical:";
        std::cout << " theta_min=" << theta_min;
        std::cout << " theta_max=" << theta_max;
        std::cout << " npred=" << npred << std::endl;
        #endif

    } // endif: performed offset angle integration

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_elliptical:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", theta_min=" << theta_min;
        std::cout << ", theta_max=" << theta_max;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return spatial integral of diffuse source model
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented.
 ***************************************************************************/
double GResponse::npred_diffuse(const GSource& source,
                                const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_NPRED_DIFFUSE,
          "Npred computation not implemented for diffuse models.");

    // Return Npred
    return 0.0;
}


/***********************************************************************//**
 * @brief Return true energy boundaries for a specific observed energy
 *
 * @param[in] obsEnergy Observed Energy.
 * @return True energy boundaries for given observed energy.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented.
 ***************************************************************************/
GEbounds GResponse::ebounds_src(const GEnergy& obsEnergy) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_EBOUNDS_SRC,
          "Npred computation not implemented for diffuse models.");

    // Allocate dummy energy boundaries
    GEbounds ebounds;

    // Return energy boundaries
    return (ebounds);
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
 * @brief Kernel for offset angle Npred integration of radial model
 *
 * @param[in] theta Radial model offset angle (radians).
 ***************************************************************************/
double GResponse::npred_radial_kern_theta::eval(const double& theta)
{
    // Initialise Npred value
    double npred = 0.0;

    // Get radial model value
    double model = m_spatial.eval(theta, m_srcEng, m_srcTime);

    // Compute sine of offset angle
    double sin_theta = std::sin(theta);

    // Setup phi integration kernel
    GResponse::npred_radial_kern_phi integrand(m_rsp,
                                               m_srcEng,
                                               m_srcTime,
                                               m_obs,
                                               m_rot,
                                               theta,
                                               sin_theta);

    // Integrate over phi
    GIntegral integral(&integrand);
    npred = integral.romberg(0.0, gammalib::twopi) * sin_theta * model;

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_radial_kern_theta::eval";
        std::cout << "(theta=" << theta << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", model=" << model;
        std::cout << ", sin_theta=" << sin_theta;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Npred integration of radial model
 *
 * @param[in] phi Azimuth angle (radians).
 ***************************************************************************/
double GResponse::npred_radial_kern_phi::eval(const double& phi)
{
    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Compute Npred for this sky direction
    double npred = m_rsp.npred(photon, m_obs);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_radial_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", cos_phi=" << cos_phi;
        std::cout << ", sin_phi=" << sin_phi;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for offset angle Npred integration of elliptical model
 *
 * @param[in] theta Radial model offset angle (radians).
 ***************************************************************************/
double GResponse::npred_elliptical_kern_theta::eval(const double& theta)
{
    // Initialise Npred value
    double npred = 0.0;

    // Compute sine of offset angle
    double sin_theta = std::sin(theta);

    // Setup phi integration kernel
    GResponse::npred_elliptical_kern_phi integrand(m_rsp,
                                                   m_spatial,
                                                   m_srcEng,
                                                   m_srcTime,
                                                   m_obs,
                                                   m_rot,
                                                   theta,
                                                   sin_theta);

    // Integrate over phi
    GIntegral integral(&integrand);
    npred = integral.romberg(0.0, gammalib::twopi) * sin_theta;

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_radial_kern_theta::eval";
        std::cout << "(theta=" << theta << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", sin_theta=" << sin_theta;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Npred integration of elliptical model
 *
 * @param[in] phi Azimuth angle (radians).
 ***************************************************************************/
double GResponse::npred_elliptical_kern_phi::eval(const double& phi)
{
    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Get elliptical model value
    double model = m_spatial.eval(m_theta, phi, m_srcEng, m_srcTime);

    // Compute Npred for this sky direction
    double npred = m_rsp.npred(photon, m_obs) * model;

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_radial_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", model=" << model;
        std::cout << ", cos_phi=" << cos_phi;
        std::cout << ", sin_phi=" << sin_phi;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}
