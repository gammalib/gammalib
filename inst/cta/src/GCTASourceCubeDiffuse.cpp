/***************************************************************************
 *        GCTASourceCubeDiffuse.cpp - CTA diffuse source cube class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTASourceCubeDiffuse.cpp
 * @brief CTA diffuse source cube class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GCTASourceCubeDiffuse.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GObservation.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GIntegral.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAResponseCube.hpp"
#include "GCTAResponse_helpers.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Method name definitions ____________________________________________ */
#define G_SET     "GCTASourceCubeDiffuse::set(GModelSpatial&, GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_PSF_INTEGRATE               //!< Integrate map over PSF (default)

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_SET

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTASourceCubeDiffuse::GCTASourceCubeDiffuse(void) : GCTASourceCube()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Diffuse source cube.
 ***************************************************************************/
GCTASourceCubeDiffuse::GCTASourceCubeDiffuse(const GCTASourceCubeDiffuse& cube) :
                       GCTASourceCube(cube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTASourceCubeDiffuse::~GCTASourceCubeDiffuse(void)
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
 * @param[in] cube Diffuse source cube.
 * @return Diffuse source cube.
 ***************************************************************************/
GCTASourceCubeDiffuse& GCTASourceCubeDiffuse::operator=(const GCTASourceCubeDiffuse& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Copy base class members
        this->GCTASourceCube::operator=(cube);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(cube);

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
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCTASourceCubeDiffuse::clear(void)
{
    // Free class members
    free_members();
    this->GCTASourceCube::free_members();

    // Initialise members
    this->GCTASourceCube::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of diffuse source cube.
 ***************************************************************************/
GCTASourceCubeDiffuse* GCTASourceCubeDiffuse::clone(void) const
{
    // Return deep copy
    return new GCTASourceCubeDiffuse(*this);
}


/***********************************************************************//**
 * @brief Set diffuse source cube for a given observation
 *
 * @param[in] model Spatial model.
 * @param[in] obs Observation.
 *
 * Sets the diffuse source cube assuming no energy dispersion and assuming
 * a negligible variation of the effective area over the size of the
 * point spread function.
 ***************************************************************************/
void GCTASourceCubeDiffuse::set(const std::string&   name,
                                const GModelSpatial& model,
                                const GObservation&  obs)
{
    // Debug option: initialise statistics
    #if defined(G_DEBUG_SET)
    int n_pixels_computed = 0;
    std::cout << "GCTASourceCubeDiffuse::set entred." << std::endl;
    #ifdef _OPENMP
    double t_start = omp_get_wtime();
    #else
    clock_t t_start = clock();
    #endif
    #endif

    // Get pointer on CTA event cube
    const GCTAEventCube* cube = dynamic_cast<const GCTAEventCube*>(obs.events());
    if (cube == NULL) {
        std::string msg = "Observation does not contain a CTA event cube.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Get pointer on CTA response cube
    const GCTAResponseCube* rsp = dynamic_cast<const GCTAResponseCube*>(obs.response());
    if (rsp == NULL) {
        std::string msg = "Observation does not contain a CTA response cube.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Get diffuse source attributes
    m_name          = name;
    GTime   obsTime = cube->time();

    // Setup empty skymap
    m_cube = cube->map();
    m_cube = 0.0;

    // Loop over all spatial bins
    for (int pixel = 0; pixel < cube->npix(); ++pixel) {

        // Get cube pixel sky direction
        GSkyDir obsDir = cube->map().inx2dir(pixel);

        // Loop over all energy layers
        for (int iebin = 0; iebin < cube->ebins(); ++iebin) {

            // Get cube layer energy
            const GEnergy& obsEng = cube->energy(iebin);

            // Determine exposure. We assume here that the exposure does
            // not vary significantly over the PSF and just compute it at
            // the pixel centre. We furthermore assume no energy dispersion,
            // and thus compute exposure using the observed energy.
            double aeff = rsp->exposure()(obsDir, obsEng);

            // Continue only if effective area is positive
            if (aeff > 0.0) {

                // Revover effective area from exposure
                aeff /= obs.livetime();

                // Compute product of PSF and diffuse map, integrated
                // over the relevant PSF area. We assume no energy
                // dispersion and thus compute the product using the
                // observed energy.
                #if defined(G_PSF_INTEGRATE)
                double psf = this->psf(rsp, &model, obsDir, obsEng, obsTime);
                #else
                double psf = model.eval(GPhoton(obsDir, obsEng, obsTime));
                #endif

                // Set cube value
                m_cube(pixel, iebin) = aeff * psf * obs.deadc(obsTime);

            } // endif: effective area was positive
        
        } // endfor: looped over all energy layers

        // Debug option: update statistics
        #if defined(G_DEBUG_SET)
        n_pixels_computed++;
        #endif
        
    } // endfor: looped over all spatial pixels

    // Debug option: show statistics
    #if defined(G_DEBUG_SET)
    #ifdef _OPENMP
    double t_elapse = omp_get_wtime()-t_start;
    #else
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    #endif
    std::cout << "  Maximum theta ................: " << theta_max << " radians" << std::endl;
    std::cout << "  Number of spatial pixels used : " << n_pixels_computed << std::endl;
    std::cout << "  CPU usage ....................: " << t_elapse << " sec" << std::endl;
    std::cout << "GCTASourceCubeDiffuse::set exit." << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Integrate PSF over diffuse model
 *
 * @param[in] rsp Response function.
 * @param[in] model Diffuse spatial model.
 * @param[in] obsDir Observed event direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * Computes the integral
 * 
 * \f[
 *    \int_0^{\delta_{\rm max}}
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm Map}(\delta, \phi) \sin \delta
 *    {\rm d}\phi {\rm d}\delta
 * \f]
 *
 * where \f${\rm Map}(\delta, \phi)\f$ is the diffuse map in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 * \f${\rm PSF}(\delta)\f$ is the azimuthally symmetric point spread
 * function.
 ***************************************************************************/
double GCTASourceCubeDiffuse::psf(const GCTAResponseCube* rsp,
                                  const GModelSpatial*    model,
                                  const GSkyDir&          obsDir,
                                  const GEnergy&          srcEng,
                                  const GTime&            srcTime) const
{
    // Set integration precision and Romberg order. These values have been
    // determined after careful testing, see
    // https://cta-redmine.irap.omp.eu/issues/1291
    const double eps_delta   = 1.0e-2;
    const int    order_delta =      5;
    const double eps_phi     = 1.0e-2;
    const int    order_phi   =      5;

    // Initialise PSF
    double psf = 0.0;

    // Compute rotation matrix to convert from PSF centred coordinate system
    // spanned by delta and phi into the reference frame of the observed
    // arrival direction given in Right Ascension and Declination.
    GMatrix ry;
    GMatrix rz;
    ry.eulery(obsDir.dec_deg() - 90.0);
    rz.eulerz(-obsDir.ra_deg());
    GMatrix rot = (ry * rz).transpose();

    // Get offset angle integration interval in radians
    double delta_min = 0.0;
    double delta_max = 1.1 * rsp->psf().delta_max();

    // Setup integration kernel. We take here the observed photon arrival
    // direction as the true photon arrival direction because the PSF does
    // not vary significantly over a small region.
    cta_psf_diffuse_kern_delta integrand(rsp, model, obsDir, srcEng, srcTime,
                                         rot, eps_phi, order_phi);

    // Integrate over PSF delta angle
    GIntegral integral(&integrand);
    integral.eps(eps_delta);
    psf = integral.romb(delta_min, delta_max, order_delta);

    // Return PSF
    return psf;
}


/***********************************************************************//**
 * @brief Print diffuse source cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing diffuse source cube information.
 ***************************************************************************/
std::string GCTASourceCubeDiffuse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTASourceCubeDiffuse ===");
        result.append("\n"+gammalib::parformat("Source name") + name());

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
 * @brief Initialise class members
 ***************************************************************************/
void GCTASourceCubeDiffuse::init_members(void)
{
    // Initialise members
    m_cube.clear();
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Point diffuse cube.
 ***************************************************************************/
void GCTASourceCubeDiffuse::copy_members(const GCTASourceCubeDiffuse& cube)
{
    // Copy members
    m_cube = cube.m_cube;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTASourceCubeDiffuse::free_members(void)
{
    // Return
    return;
}
