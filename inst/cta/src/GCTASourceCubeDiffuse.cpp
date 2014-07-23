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

/* __ Method name definitions ____________________________________________ */
#define G_SET     "GCTASourceCubeDiffuse::set(GModelSpatial&, GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_PSF_INTEGRATE               //!< Integrate map over PSF (default)

/* __ Debug definitions __________________________________________________ */

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
    GSkyDir centre;
    double  theta_max = 99.0;

    // Store spatial model parameter values
    m_pars.clear();
    for (int i = 0; i < model.size(); ++i) {
        m_pars.push_back(model[i].value());
    }

    // Get attributes for radial or elliptical model
    const GModelSpatialRadial* radial = dynamic_cast<const GModelSpatialRadial*>(&model);
    if (radial != NULL) {
        centre    = radial->dir();
        theta_max = radial->theta_max();
    }
    else {
        const GModelSpatialElliptical* elliptical = dynamic_cast<const GModelSpatialElliptical*>(&model);
        if (radial != NULL) {
            centre    = elliptical->dir();
            theta_max = elliptical->theta_max();
        }
    }

    // Setup empty skymap
    m_cube = cube->map();
    m_cube = 0.0;

    // Loop over all spatial bins
    for (int pixel = 0; pixel < cube->npix(); ++pixel) {

        // Get cube pixel sky direction
        GSkyDir obsDir = cube->map().inx2dir(pixel);

        // Compute distance to model centre in radians
        double theta = centre.dist(obsDir);
        if (theta <= theta_max) {

            // Loop over all energy bins
            for (int iebin = 0; iebin < cube->ebins(); ++iebin) {

                // Get cube layer energy
                const GEnergy& obsEng = cube->energy(iebin);

                // Determine deadtime corrected effective area. We assume
                // here that the effective area does not vary significantly
                // over the PSF and just compute it at the pixel centre. We
                // furthermore assume no energy dispersion, and thus compute
                // effective area using the observed energy.
                double aeff = rsp->exposure()(obsDir, obsEng) *
                              obs.deadc(obsTime) /
                              obs.ontime();

                // Continue only if effective area is positive
                if (aeff > 0.0) {

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
                    m_cube(pixel, iebin) = aeff * psf;

                } // endif: effective area was positive
            
            } // endfor: looped over all energy layers
        
        } // endif: pixel was relevant

    } // endfor: looped over all spatial pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return parameter value
 *
 * @param[in] index Parameter index.
 * @return Parameter value.
 *
 * Returns the parameter value for the specified parameter @p index..
 ***************************************************************************/
double GCTASourceCubeDiffuse::par(const int& index) const
{
    //TODO: Range checking
    
    // Return
    return (m_pars[index]);
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
    psf_kern_delta integrand(rsp, model, obsDir, srcEng, srcTime, rot);

    // Integrate over PSF delta angle
    GIntegral integral(&integrand);
    integral.eps(1.0e-2);
    psf = integral.romb(delta_min, delta_max);

    // Return PSF
    return psf;
}


/***********************************************************************//**
 * @brief Kernel for PSF integration of diffuse model
 *
 * @param[in] delta PSF offset angle (radians).
 * @return Azimuthally integrated product between PSF and model.
 *
 * Computes the azimuthally integrated product of point spread function and
 * diffuse model intensity. As the PSF is azimuthally symmetric, it is
 * not included in the azimuthally integration, but just multiplied on the
 * azimuthally integrated map. The method returns thus
 *
 * \f[
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm Map}(\delta, \phi) \sin \delta {\rm d}\phi
 * \f]
 *
 * where \f${\rm Map}(\delta, \phi)\f$ is the diffuse map in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 *
 * The method adjust the azimuth integration precision as function of the PSF
 * offset angle \f$\delta\f$, with a higher precision required for regions
 * where the PSF value is larger, and lower precision required in the PSF
 * tails.
 *
 * @todo Perform a detailed study to adjust the integration precision
 *       evolution as function of @p delta argument.
 ***************************************************************************/
double GCTASourceCubeDiffuse::psf_kern_delta::eval(const double& delta)
{
    // Get PSF for this delta
    double value = m_rsp->psf()(m_srcDir, delta, m_srcEng);

    // Initialize spatially integrated map value
    double map = 0.0;

    // If we're at the PSF peak the model is zero (due to the sin(delta)
    // term. We thus only integrate for positive deltas, and of course only
    // for positive PSF values.
    if (value > 0.0 && delta > 0.0) {

        // Compute sine and cosine of delta
        double sin_delta = std::sin(delta);
        double cos_delta = std::cos(delta);

        // Setup kernel for azimuthal integration of the spatial map
        psf_kern_phi integrand(m_model, m_srcEng, m_srcTime, m_rot,
                               sin_delta, cos_delta);

        // Set the requested integration precision. This goes from 1.0e-2
        // for regions where the PSF is small to 1.0e-4 for regions where
        // the PSF is high
        //double eps = 1.0e-1 / (990.0 * value / m_psf_max + 10.0);
        double eps = 1.0e-2;

        // Azimuthally integrate spatial map
        GIntegral integral(&integrand);
        integral.eps(eps);
        map = integral.romb(0.0, gammalib::twopi) * sin_delta;

    }

    // Multiply map with PSF value
    value *= map;

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTASourceCubeDiffuse::psf_kern_delta::eval";
        std::cout << "(delta=" << delta << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return kernel value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for map integration of diffuse model
 *
 * @param[in] phi Azimuth angle (radians).
 * @return Diffuse model value.
 *
 * Computes the value of the diffuse model at the position (delta,phi) given
 * in point spread function coordinates. The transformation from point
 * spread function coordinates into sky coordinates is done using a rotation
 * matrix that is pre-computed on entry.
 ***************************************************************************/
double GCTASourceCubeDiffuse::psf_kern_phi::eval(const double& phi)
{
    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_delta, sin_phi*m_sin_delta, m_cos_delta);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Compute map value this sky direction
    double value = m_model->eval(photon); 

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTASourceCubeDiffuse::psf_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return kernel value
    return value;
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
    m_pars.clear();
   
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
    m_pars = cube.m_pars;

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
