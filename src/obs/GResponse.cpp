/***************************************************************************
 *                GResponse.cpp - Response abstract base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
#include "GModelSpatialPointSource.hpp"
#include "GModelRadial.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF_EXTENDED         "GResponse::irf_extended(GInstDir&, GEnergy&,"\
           " GTime&, GModelExtendedSource&, GEnergy&, GTime&, GObservation&)"
#define G_IRF_DIFFUSE           "GResponse::irf_diffuse(GInstDir&, GEnergy&,"\
            " GTime&, GModelDiffuseSource&, GEnergy&, GTime&, GObservation&)"
#define G_NPRED_EXTENDED   "GResponse::npred_extended(GModelExtendedSource&,"\
                                          " GEnergy&, GTime&, GObservation&)"
#define G_NPRED_DIFFUSE      "GResponse::npred_diffuse(GModelDiffuseSource&,"\
                                          " GEnergy&, GTime&, GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_NPRED_EXTENDED                    //!< Debug npred_extended


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
 *
 * Note that this method applies the deadtime correction, so that the
 * response function can be directly multiplied by the exposure time (also
 * known as ontime).
 ***************************************************************************/
double GResponse::irf(const GEvent&       event,
                      const GSource&      source,
                      const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Is model a point source?
    if (dynamic_cast<const GModelSpatialPointSource*>(source.model()) != NULL) {
        irf = irf_ptsrc(event, source, obs);
    }

    // Is model an extended source?
    else if (dynamic_cast<const GModelRadial*>(source.model()) != NULL) {
        irf = irf_extended(event, source, obs);
    }

    // ... otherwise we have a diffuse model
    else {
        irf = irf_diffuse(event, source, obs);
    }

    // Apply deadtime correction
    irf *= obs.deadc(source.time());

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
 * point source. If the source is not a point source, the method returns 0.
 ***************************************************************************/
double GResponse::irf_ptsrc(const GEvent&       event,
                            const GSource&      source,
                            const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get point source spatial model
    const GModelSpatialPointSource* src =
          dynamic_cast<const GModelSpatialPointSource*>(source.model());

    // Continue only if model is valid
    if (src != NULL) {
    
        // Set Photon
        GPhoton photon(src->dir(), source.energy(), source.time());
    
        // Compute IRF
        irf = this->irf(event, photon, obs);

    }

    // Return IRF
    return irf;
}


/***********************************************************************//**
 * @brief Return value of extended source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_extended(const GEvent&       event,
                               const GSource&      source,
                               const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_EXTENDED,
          "Extended IRF not implemented for this class.");

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
          "Diffuse IRF not implemented for this class.");

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
 * a given source and a particular observation.
 *
 * The method applies the deadtime correction, so that the result can be
 * directly multiplied by the exposure time (also known as ontime).
 ***************************************************************************/
double GResponse::npred(const GSource& source, const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Is model a point source?
    if (dynamic_cast<const GModelSpatialPointSource*>(source.model()) != NULL) {
        npred = npred_ptsrc(source, obs);
    }

    // Is model an extended source?
    else if (dynamic_cast<const GModelRadial*>(source.model()) != NULL) {
        npred = npred_extended(source, obs);
    }

    // ... otherwise we have a diffuse model
    else {
        npred = npred_diffuse(source, obs);
    }

    // Apply deadtime correction
    npred *= obs.deadc(source.time());

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
 * region of interest. If the source is not a point source, the method
 * returns 0.
 ***************************************************************************/
double GResponse::npred_ptsrc(const GSource& source,
                              const GObservation& obs) const
{
    // Initialise Npred
    double npred = 0.0;

    // Get point source spatial model
    const GModelSpatialPointSource* src =
          dynamic_cast<const GModelSpatialPointSource*>(source.model());

    // Continue only if model is valid
    if (src != NULL) {
    
        // Set Photon
        GPhoton photon(src->dir(), source.energy(), source.time());

        // Compute Npred
        npred = this->npred(photon, obs);

    }

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return ROI integral of extended source model
 *
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Integral of extended source model over ROI.
 *
 * This method returns the spatial integral of an extended source model over
 * the region of interest. If the model is not a radial model, the method
 * returns 0.
 ***************************************************************************/
double GResponse::npred_extended(const GSource& source,
                                 const GObservation& obs) const
{
    // Initialise Npred value
    double npred = 0.0;

    // Get radial model
    const GModelRadial* src = dynamic_cast<const GModelRadial*>(source.model());

    // Continue only if model is valid
    if (src != NULL) {
    
        // Compute rotation matrix to convert from native coordinates given
        // by (theta,phi) into celestial coordinates.
        GMatrix ry;
        GMatrix rz;
        GMatrix rot;
        ry.eulery(src->dec() - 90.0);
        rz.eulerz(-src->ra());
        rot = transpose(ry * rz);

        // Set offset angle integration range
        double theta_min = 0.0;
        double theta_max = src->theta_max();

        // Perform offset angle integration if interval is valid
        if (theta_max > theta_min) {

            // Setup integration kernel
            GResponse::npred_kern_theta integrand(this,
                                                  src,
                                                  &(source.energy()),
                                                  &(source.time()),
                                                  &obs,
                                                  &rot);

            // Integrate over theta
            GIntegral integral(&integrand);
            npred = integral.romb(theta_min, theta_max);

            // Compile option: Show integration results
            #if defined(G_DEBUG_NPRED_EXTENDED)
            std::cout << "GResponse::npred_extended:";
            std::cout << " theta_min=" << theta_min;
            std::cout << " theta_max=" << theta_max;
            std::cout << " npred=" << npred << std::endl;
            #endif

        } // endif: offset angle range was valid

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (isnotanumber(npred) || isinfinite(npred)) {
            std::cout << "*** ERROR: GResponse::npred_extended:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", theta_min=" << theta_min;
            std::cout << ", theta_max=" << theta_max;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: radial model was valid

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
          "Method for extended source not implemented for this class.");

    // Return Npred
    return 0.0;
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
 * @brief Kernel for offset angle Npred integration
 *
 * @param[in] theta Radial model offset angle (radians).
 ***************************************************************************/
double GResponse::npred_kern_theta::eval(double theta)
{
    // Initialise Npred value
    double npred = 0.0;

    // Get radial model value
    double model = m_radial->eval(theta);

    // Compute sine of offset angle
    double sin_theta = std::sin(theta);

    // Setup phi integration kernel
    GResponse::npred_kern_phi integrand(m_rsp,
                                        m_srcEng,
                                        m_srcTime,
                                        m_obs,
                                        m_rot,
                                        theta,
                                        sin_theta);

    // Integrate over phi
    GIntegral integral(&integrand);
    npred = integral.romb(0.0, twopi) * sin_theta * model;

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_kern_theta::eval";
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
 * @brief Kernel for azimuth angle Npred integration
 *
 * @param[in] phi Azimuth angle (radians).
 ***************************************************************************/
double GResponse::npred_kern_phi::eval(double phi)
{
    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    GPhoton photon(srcDir, *m_srcEng, *m_srcTime);

    // Compute Npred for this sky direction
    double npred = m_rsp->npred(photon, *m_obs);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GResponse::npred_kern_phi::eval";
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
