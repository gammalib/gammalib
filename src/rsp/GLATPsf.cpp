/***************************************************************************
 *           GLATPsf.cpp  -  GLAST LAT Response class Psf methods          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GTools.hpp"
#include "GException.hpp"
#include "GVector.hpp"
#include "GLATResponse.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsDblImage.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_INIT_PSF "GLATResponse::init_psf()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const double scale_factor_log = 39.810717;
const double scale_base_log   = 0.15848932;
const double scale_front_c0   = 3.77e-4;
const double scale_front_c1   = 5.8e-2;
const double scale_back_c0    = 1.3e-3;
const double scale_back_c1    = 9.6e-2;
const int    angle_num        = 5000;    // Number of angles for PSF
const double angle_min        = 0.0001;  // Minimum angular separation (rad)
const double angle_bin        = 0.01;    // Angular separation binning (rad)
const double max_sep          = pihalf;  // Maximum angular separation (rad)


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of point spread function
 *
 * @param[in] obsDir Observed photon direction
 * @param[in] obsEng Observed energy of photon
 * @param[in] srcDir True photon direction
 * @param[in] srcEng True energy of photon
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
 * @param[in] instPosAng Instrument position angle
 * @param[in] time Photon arrival time
 ***************************************************************************/
double GLATResponse::psf(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time)
{
    // Return PSF value
    return 0.0;
}


/***********************************************************************//**
 * @brief Return value of point spread function
 *
 * @param[in] delta Angular separation between true and observed photon 
 *            direction (radians).
 * @param[in] logE Logarithm of the true photon energy (MeV).
 * @param[in] ctheta Cosine of the zenith angle with respect to the pointing 
 *            axis.
 *
 * Returns the GLAST LAT PSF value for a given angular separation between
 * true and observed photon direction, true photon energy and cosine of
 * the zenith angle.
 ***************************************************************************/
double GLATResponse::psf(const double& delta, const double& logE, 
                         const double& ctheta)
{
    // Get sigma value
    double sigma = m_psf_bins.interpolate(logE, ctheta, m_sigma);

    // Get PSF scale
    double scale = psf_scale_log(logE);

    // Compute stretch factor
    double stretch = sigma * scale;
    
    // Compute maximum allowed angular separation
    double delta_max = m_psf_angle_max * stretch;
    if (delta_max >= max_sep)
        delta_max = max_sep;
    
    // If maximum angular separation is exceeded then return 0.0
    double psf;
    if (delta > delta_max)
        psf = 0.0;
        
    // ... otherwise interpolate PSF value
    else {
    
        // Get PSF value from table
        double r    = delta / stretch;  // real angular distance in radians
        m_psf_angle.set_value(r);       // interpolate array
        int    inx1 = m_psf_angle.inx_left();
        int    inx2 = m_psf_angle.inx_right();
        double wgt1 = m_psf_angle.wgt_left();
        double wgt2 = m_psf_angle.wgt_right();
        psf = wgt1 * m_psf_bins.interpolate(logE, ctheta, m_psf, inx1, angle_num) +
              wgt2 * m_psf_bins.interpolate(logE, ctheta, m_psf, inx2, angle_num);

        // Get normalization value
        double norm = m_psf_bins.interpolate(logE, ctheta, m_norm) /
                      (stretch * stretch);

        // Normalize PSF
        psf *= norm;
        
    } // endelse: PSF normalization was required

    // Return PSF value
    return psf;
}


/***********************************************************************//**
 * @brief Return vector of point spread function values
 *
 * @param[in] delta Vector of angular separation between true and observed 
 *            photon direction (radians).
 * @param[in] logE Logarithm of the true photon energy (MeV).
 * @param[in] ctheta Cosine of the zenith angle with respect to the pointing 
 *            axis.
 *
 * Returns the GLAST LAT PSF value for a given angular separation between
 * true and observed photon direction, true photon energy and cosine of
 * the zenith angle.
 ***************************************************************************/
GVector GLATResponse::psf(const GVector& delta, const double& logE, 
                          const double& ctheta)
{
    // Get sigma value
    double sigma = m_psf_bins.interpolate(logE, ctheta, m_sigma);

    // Get PSF scale
    double scale = psf_scale_log(logE);

    // Compute stretch factor
    double stretch = sigma * scale;

    // Get normalization value
    double norm = m_psf_bins.interpolate(logE, ctheta, m_norm) /
                  (stretch * stretch);
    
    // Compute maximum allowed angular separation
    double delta_max = m_psf_angle_max * stretch;
    if (delta_max >= max_sep)
        delta_max = max_sep;

    // Allocate result vector
    GVector psf(delta.size());
    
    // Loop over all angular separations
    for (int i = 0; i < delta.size(); ++i) {

        // If maximum angular separation is exceeded then stop
        if (delta(i) > delta_max)
            break;
        
        // Get PSF value from table
        double r    = delta(i) / stretch;      // real angular distance in radians
        m_psf_angle.set_value(r);              // interpolate array
        int    inx1 = m_psf_angle.inx_left();
        int    inx2 = m_psf_angle.inx_right();
        double wgt1 = m_psf_angle.wgt_left();
        double wgt2 = m_psf_angle.wgt_right();
        psf(i) = wgt1 * m_psf_bins.interpolate(logE, ctheta, m_psf, inx1, angle_num) +
                 wgt2 * m_psf_bins.interpolate(logE, ctheta, m_psf, inx2, angle_num);

        // Normalize PSF
        psf(i) *= norm;
        
    } // endfor: looped over angular separations

    // Return PSF vector
    return psf;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise PSF
 *
 * Initialise PSF for response calculation. Initialisation consists of 
 * building the PSF response vectors for a set of energy and zenith angle
 * values.
 *
 * The PSF is built using information found in the 'RPSF' table of the
 * calibration database.
 ***************************************************************************/
void GLATResponse::psf_init(void)
{
    // Build filename
    std::string filename = "psf_"  + m_rspname + "_" + m_rsptype + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer to PSF HDU
    GFitsHDU* hdu = file.hdu("RPSF");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_PSF, "RPSF");

    // Get the energy and cos(theta) bins
    m_psf_bins.load(hdu);
    
    // Get the PSF characterisation data
    GVector v_ncore = get_fits_vector(hdu, "NCORE");
    GVector v_sigma = get_fits_vector(hdu, "SIGMA");
    GVector v_gcore = get_fits_vector(hdu, "GCORE");
    GVector v_gtail = get_fits_vector(hdu, "GTAIL");

    // Free old PSF memory
    if (m_psf   != NULL) delete [] m_psf;
    if (m_norm  != NULL) delete [] m_norm;
    if (m_sigma != NULL) delete [] m_sigma;

    // Allocate memory for PSF
    int bin_size = m_psf_bins.num_energy() * m_psf_bins.num_ctheta();
    int psf_size = bin_size * angle_num;
    m_psf        = new double[psf_size];
    m_norm       = new double[bin_size];
    m_sigma      = new double[bin_size];

    // Setup PSF vector axis
    GVector xaxis(angle_num);
    GVector u(angle_num);
    double x = angle_min;
    for (int i = 0; i < angle_num; ++i, x += angle_bin) {
        xaxis(i) = x;
        u(i)     = 0.5*x*x;
    }
    
    // Setup angle nodes
    m_psf_angle.nodes(xaxis);

    // Get maximum angle
    m_psf_angle_max = xaxis(angle_num-1);

    // Loop over all energy bins
    for (int ie = 0; ie < m_psf_bins.num_energy(); ++ie) {

        // Get mean energy of bin
        double energy = m_psf_bins.energy(ie);

        // Get scale at specified energy
        double scale = psf_scale(energy);

        // Loop over all cos theta bins
        for (int ic = 0; ic < m_psf_bins.num_ctheta(); ++ic) {

            // Extract PSF parameters
            int    inx   = m_psf_bins.index(ie,ic);
            double ncore = v_ncore(inx);
            double sigma = v_sigma(inx);
            double gcore = v_gcore(inx);
            double gtail = v_gtail(inx);

            // Evaluate ntail
            double denom = psf_base_value(10.0, gtail);
            double ntail = (denom != 0.0) 
                           ? ncore * psf_base_value(10.0, gcore) / denom : 0.0;

            // Calculate PSF dependence
            GVector psf = ncore * psf_base_vector(u, gcore) +
                          ntail * psf_base_vector(u, gtail);

            // Normalize PSF vector to identical amplitudes at smallest angular
            // distance. This allows for accurate interpolation afterwards
            double norm = 1.0 / psf(0);
            psf *= norm;

            // Evaluate angular separation
            double  stretch = sigma * scale;
            GVector delta   = xaxis * stretch;

            // Calculate PSF normalization using trapezoid rule. Note that 
            // normalization to unity would required to multiply the binsize
            // 'angle_bin' with the 'stretch' factor. We would then have to 
            // write:
            //  sum *= 0.5 * twopi * stretch * angle_bin;
            //  norm = 1.0 / sum
            // However, to keep the normalization values rather independent
            // of energy, we multiply the normalization factors by
            // 'stretch^2'. In the PSF routines we then have to devide
            // 'norm' by 'stretch^2'.
            double sum = 0.0;
            for (int i = 1; i < angle_num; ++i) {
                if (delta(i) >= max_sep)
                    psf(i) = 0.0;
                else
                    sum += (psf(i)*sin(delta(i)) + psf(i-1)*sin(delta(i-1)));
            }
            sum *= 0.5 * twopi * angle_bin;
            norm = stretch / sum;

            // Store PSF information
            m_norm[inx]  = norm;
            m_sigma[inx] = sigma;

            // Store PSF vector
            int ipsf = inx * angle_num;
            for (int i = 0; i < angle_num; ++i, ++ipsf)
                m_psf[ipsf] = psf(i);

        } // endfor: looped over cos theta
    } // endfor: looped over energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append PSF HDUs tu FOTS file
 *
 * @param[in] file FITS file into which the PSF HDUs will be appended
 *
 * Append 4 HDUs to the FITS file:
 * PBOUNDS (PSF energy and zenith angle boundaries)
 * PNORM (PSF normalization array)
 * PSIGMA (PSF sigma array)
 * PSF (PSF vectors)
 ***************************************************************************/
void GLATResponse::psf_append(GFits& file) const
{
    // Get PSF boundary table
    GFitsHDU hdu_psf_bounds;
    m_psf_bins.save(&hdu_psf_bounds);
    hdu_psf_bounds.extname("PBOUNDS");
    
    // Build PSF tables and images
    int naxes_psf2[] = {m_psf_bins.num_energy(), m_psf_bins.num_ctheta()};
    int naxes_psf3[] = {angle_num, m_psf_bins.num_energy(), 
                        m_psf_bins.num_ctheta()};
    GFitsDblImage image_psf_norm(2, naxes_psf2, m_norm);
    GFitsDblImage image_psf_sigma(2, naxes_psf2, m_sigma);
    GFitsDblImage image_psf_psf(3, naxes_psf3, m_psf);

    // Construct PSF HDUs
    GFitsHDU hdu_psf_norm(image_psf_norm);
    GFitsHDU hdu_psf_sigma(image_psf_sigma);
    GFitsHDU hdu_psf_psf(image_psf_psf);
    hdu_psf_norm.extname("PNORM");
    hdu_psf_sigma.extname("PSIGMA");
    hdu_psf_psf.extname("PSF");

    // Append HDUs to FITS file
    file.append_hdu(hdu_psf_bounds);
    file.append_hdu(hdu_psf_norm);
    file.append_hdu(hdu_psf_sigma);
    file.append_hdu(hdu_psf_psf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns PSF scale
 *
 * @param[in] energy Energy at which the PSF scale should be computed (MeV).
 *
 * The PSF size scales as
 *   \f$scale=\sqrt{c_0^2 + (c_1 \times t)^2}\f$
 * with
 *   \f$t = (0.01*E)^{-0.8}\f$
 * where E is the energy in MeV.
 * The coefficients \f$c_0\f$ and  \f$c_1\f$ dependend on the section:
 * Front: \f$c_0 = 3.77 \times 10^{-4}\f$, \f$c_1= 5.8 \times 10^{-2}\f$
 * Back:  \f$c_0 = 1.3  \times 10^{-3}\f$, \f$c_1= 9.6 \times 10^{-2}\f$
 ***************************************************************************/
double GLATResponse::psf_scale(const double& energy) const
{
    // Get section dependent coefficients
    const double *c0;
    const double *c1;
    if (m_section == 0) {   // front
        c0 = &scale_front_c0;
        c1 = &scale_front_c1;
    }
    else {                  // back
        c0 = &scale_back_c0;
        c1 = &scale_back_c1;
    }

    // Angular resolution scales as a powerlaw with index -0.80
    double t = pow(0.01*energy,-0.80);

    // Get scaling
    double c1t   = *c1 * t;
    double scale = sqrt(*c0 * *c0 + c1t * c1t);

    // Return scaling
    return scale;
}


/***********************************************************************//**
 * @brief Returns PSF scale for logarithmic energy
 *
 * @param[in] logE Base 10 logarithmic energy at which the PSF scale should 
 *            be computed (MeV).
 *
 * The PSF size scales as
 *   \f$scale=\sqrt{c_0^2 + (c_1 \times t)^2}\f$
 * with
 *   \f$t = (0.01*(10^E))^{-0.8}\f$
 * where E is the energy in MeV.
 * The coefficients \f$c_0\f$ and  \f$c_1\f$ dependend on the section:
 * Front: \f$c_0 = 3.77 \times 10^{-4}\f$, \f$c_1= 5.8 \times 10^{-2}\f$
 * Back:  \f$c_0 = 1.3  \times 10^{-3}\f$, \f$c_1= 9.6 \times 10^{-2}\f$
 ***************************************************************************/
double GLATResponse::psf_scale_log(const double& logE) const
{
    // Get section dependent coefficients
    const double *c0;
    const double *c1;
    if (m_section == 0) {   // front
        c0 = &scale_front_c0;
        c1 = &scale_front_c1;
    }
    else {                  // back
        c0 = &scale_back_c0;
        c1 = &scale_back_c1;
    }

    // Angular resolution scales as a powerlaw with index -0.80
    double t = scale_factor_log * pow(scale_base_log, logE);

    // Get scaling
    double c1t   = *c1 * t;
    double scale = sqrt(*c0 * *c0 + c1t * c1t);

    // Return scaling
    return scale;
}


/***********************************************************************//**
 * @brief Returns PSF basis function
 *
 * @param[in] u angular distance squared
 * @param[in] gamma Powerlaw index
 *
 * This method returns the PSF basis function for a single value.
 * The PSF basis functions is defined by
 *   \f$psf(u,\gamma) = (1.0 - 1.0/\gamma) (1.0 + u/\gamma)^{-\gamma}\f$
 * where
 *   \f$u\f$ is an angle^2 (typically the scaled angular distance squared)
 *   \f$\gamma\f$ is a powerlaw index
 ***************************************************************************/
double GLATResponse::psf_base_value(const double& u, const double& gamma) const
{
    // Evaluate PSF basis function
    double invgamma = 1.0/gamma;
    double c0       = 1.0 - invgamma;
    double c1       = 1.0 + u * invgamma;
    double result   = c0 * pow(c1,-gamma);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns PSF basis function
 *
 * @param[in] u Vector of angular distance squared
 * @param[in] gamma Powerlaw index
 *
 * See GLATResponse::psf_base_value for more information.
 ***************************************************************************/
GVector GLATResponse::psf_base_vector(const GVector& u, const double& gamma) const
{
    // Evaluate PSF basis function
    double  invgamma = 1.0/gamma;
    double  c0       = 1.0 - invgamma;
    GVector c1       = 1.0 + u * invgamma;
    GVector result   = c0 * pow(c1,-gamma);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Initialise PSF class members
 ***************************************************************************/
void GLATResponse::psf_init_members(void)
{
    // Initialise PSF members
    m_psf_angle_max  = 0.0;
    m_psf            = NULL;
    m_norm           = NULL;
    m_sigma          = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy PSF class members
 ***************************************************************************/
void GLATResponse::psf_copy_members(const GLATResponse& rsp)
{
    // Copy attributes
    m_psf_bins      = rsp.m_psf_bins;
    m_psf_angle     = rsp.m_psf_angle;
    m_psf_angle_max = rsp.m_psf_angle_max;
    
    // Copy arrays
    if (rsp.m_psf != NULL) {
        int size = m_psf_bins.num_energy() * m_psf_bins.num_ctheta() * angle_num;
        m_psf    = new double[size];
        memcpy(m_psf, rsp.m_psf, size*sizeof(double));
    }
    if (rsp.m_norm != NULL) {
        int size = m_psf_bins.num_energy() * m_psf_bins.num_ctheta();
        m_norm   = new double[size];
        memcpy(m_norm, rsp.m_norm, size*sizeof(double));
    }
    if (rsp.m_sigma != NULL) {
        int size = m_psf_bins.num_energy() * m_psf_bins.num_ctheta();
        m_sigma  = new double[size];
        memcpy(m_sigma, rsp.m_sigma, size*sizeof(double));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
***************************************************************************/
void GLATResponse::psf_free_members(void)
{
    // Free PSF memory
    if (m_psf       != NULL) delete [] m_psf;
    if (m_norm      != NULL) delete [] m_norm;
    if (m_sigma     != NULL) delete [] m_sigma;

    // Signal that PSF memory is free
    m_psf   = NULL;
    m_norm  = NULL;
    m_sigma = NULL;

    // Return
    return;
}


