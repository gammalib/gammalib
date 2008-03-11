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
const double scale_front_c0 = 3.77e-4;
const double scale_front_c1 = 5.8e-2;
const double scale_back_c0  = 1.3e-3;
const double scale_back_c1  = 9.6e-2;

const double pi          = 3.1415926535897931159979635;
const double pihalf      = pi/2.0;
const double twopi       = 2.0*pi;
const double deg2rad     = 0.0174532925199432954743717;
const double rad2deg     = 57.295779513082322864647722;
const double ln10        = 2.3025850929940459010936138;


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of point spread function
 *
 * @param obsDir Observed photon direction
 * @param obsEng Observed energy of photon
 * @param srcDir True photon direction
 * @param srcEng True energy of photon
 * @param instPntDir Instrument pointing direction (e.g. z-axis)
 * @param instPosAng Instrument position angle
 * @param time Photon arrival time
 ***************************************************************************/
double GLATResponse::psf(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time)
{
    // Return PSF value
    return 0.0;
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

    // Set binning in angle, energy and cos theta
    m_psf_angle_num  = 2000;
    m_psf_angle_min  = 0.01;
    m_psf_angle_bin  = 0.01;

    // Free old PSF memory
    if (m_psf   != NULL) delete [] m_psf;
    if (m_norm  != NULL) delete [] m_norm;
    if (m_sigma != NULL) delete [] m_sigma;

    // Allocate memory for PSF
    int bin_size = m_psf_bins.num_energy() * m_psf_bins.num_ctheta();
    int psf_size = bin_size * m_psf_angle_num;
    m_psf        = new double[psf_size];
    m_norm       = new double[bin_size];
    m_sigma      = new double[bin_size];

    // Setup PSF vector axis
    GVector xaxis(m_psf_angle_num);
    GVector u(m_psf_angle_num);
    double x = m_psf_angle_min;
    for (int i = 0; i < m_psf_angle_num; ++i, x += m_psf_angle_bin) {
        xaxis(i) = x;
        u(i)     = 0.5*x*x;
    }

    // Get maximum angle
    m_psf_angle_max = xaxis(m_psf_angle_num-1);

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

            // Calculate PSF normalization factor
            double sum = 0.0;
            for (int i = 0; i < m_psf_angle_num; ++i) {
                if (delta(i) >= pihalf)
                    break;
                sum += psf(i) * sin(delta(i));
            }
            sum *= twopi;
            norm = scale / sum;

            // Store PSF information
            m_norm[inx]  = norm;
            m_sigma[inx] = sigma;

            // Store PSF vector
            int ipsf = inx * m_psf_angle_num;
            for (int i = 0; i < m_psf_angle_num; ++i, ++ipsf)
                m_psf[ipsf] = psf(i);

        } // endfor: looped over cos theta
    } // endfor: looped over energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save PSF
 *
 * @param file FITS file into which the PSF will be saved
 *
 * Writes 4 HDUs to the FITS file:
 * PBOUNDS (PSF energy and zenith angle boundaries)
 * PSF
 * NORM
 * SIGMA
 ***************************************************************************/
void GLATResponse::psf_save(GFits& file) const
{
    // Get PSF boundary table
    GFitsHDU hdu_psf_bounds;
    m_psf_bins.save(&hdu_psf_bounds);
    hdu_psf_bounds.extname("PBOUNDS");
    
    // Build PSF tables and images
    int           naxes_psf2[] = {m_psf_bins.num_energy(),  m_psf_bins.num_ctheta()};
    int           naxes_psf3[] = {m_psf_angle_num,  m_psf_bins.num_energy(),
                                   m_psf_bins.num_ctheta()};
    GFitsDblImage image_psf_psf(3, naxes_psf3, m_psf);
    GFitsDblImage image_psf_norm(2, naxes_psf2, m_norm);
    GFitsDblImage image_psf_sigma(2, naxes_psf2, m_sigma);

    // Construct PSF HDUs
    GFitsHDU hdu_psf_psf(image_psf_psf);
    GFitsHDU hdu_psf_norm(image_psf_norm);
    GFitsHDU hdu_psf_sigma(image_psf_sigma);
    hdu_psf_psf.extname("PSF");
    hdu_psf_norm.extname("NORM");
    hdu_psf_sigma.extname("SIGMA");

    // Append HDUs to FITS file
    file.append_hdu(hdu_psf_bounds);
    file.append_hdu(hdu_psf_psf);
    file.append_hdu(hdu_psf_norm);
    file.append_hdu(hdu_psf_sigma);

    // Save FITS file
    //file.save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get PSF value
 *
 * @param delta angular separation between true and observed photon direction
 * @param logE logarithm of the true photon energy
 * @param ctheta cosine of the zenith angle with respect to the pointing axis
 *
 * Returns the GLAST LAT PSF value
 ***************************************************************************/
double GLATResponse::psf_get(const double& delta, const double& logE, 
                             const double& ctheta)
{
    // Interpolate sigma and norm tables
    double sigma = 1.0;  // Perform bi-linear interpolation in logE and ctheta
    double norm  = 1.0;  // Perform bi-linear interpolation in logE and ctheta

    // Evaluate scale
    double energy = pow(10.0,logE);    // DUMMY: t=(0.01*(10.0^logE))^-0.8
    double scale  = psf_scale(energy);

    // Get stretch factor
    double stretch = sigma * scale;

    // Get PSF value from table
    double r   = delta * deg2rad / stretch;  // in radians
    int    inx = (int)((r - m_psf_angle_min) / m_psf_angle_bin);  // Perform linear interpolation
    double psf = 1.0;

    // Normalize PSF
    psf *= norm;

    // Return PSF value
    return psf;
}


/***********************************************************************//**
 * @brief Returns PSF scale
 *
 * @param energy
 *
 * The PSF size scales as
 *   scale = sqrt(c0^2 + (c1*t)^2)
 * with
 *   t = (0.01*E)^-0.8
 * where E is the energy in MeV.
 * The coefficients c0 and c1 dependend on the section:
 * Front: c0 = 3.77e-4, c1 = 5.8e-2
 * Back:  c0 = 1.3e-3,  c1 = 9.6e-2
 ***************************************************************************/
double GLATResponse::psf_scale(const double& energy)
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
 * @brief Returns PSF basis function
 *
 * @param u angular distance squared
 * @param gamma Powerlaw index
 *
 * This method returns the PSF basis function for a single value.
 * The PSF basis functions is defined by
 *   psf(u,gamma) = (1.0 - 1.0/gamma) * (1.0 + u/gamma)^(-gamma)
 * where
 *   u     is an angle^2 (typically the scaled angular distance squared)
 *   gamma is a powerlaw index
 ***************************************************************************/
double GLATResponse::psf_base_value(const double& u, const double& gamma)
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
 * @param u Vector of angular distance squared
 * @param gamma Powerlaw index
 *
 * This method returns the PSF basis function for a single value.
 * The PSF basis functions is defined by
 *   psf(u,gamma) = (1.0 - 1.0/gamma) * (1.0 + u/gamma)^(-gamma)
 * where
 *   u     is an angle^2 (typically the scaled angular distance squared)
 *   gamma is a powerlaw index
 ***************************************************************************/
GVector GLATResponse::psf_base_vector(const GVector& u, const double& gamma)
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
    m_psf_angle_num  = 0;
    m_psf_angle_min  = 0.0;
    m_psf_angle_max  = 0.0;
    m_psf_angle_bin  = 0.0;
    m_psf            = NULL;
    m_norm           = NULL;
    m_sigma          = NULL;

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


