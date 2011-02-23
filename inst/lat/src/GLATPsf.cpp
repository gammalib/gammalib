/***************************************************************************
 *             GLATPsf.cpp  -  Fermi LAT point spread function             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATPsf.cpp
 * @brief Fermi LAT point spread function class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATPsf.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_PSF1           "GLATPsf::operator() (double&, double&, double&)"
#define G_READ                             "GLATPsf::read(const GFits* file)"
#define G_READ_PSF                           "GLATPsf::read_psf(GFitsTable*)"
#define G_READ_SCALE                       "GLATPsf::read_scale(GFitsTable*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_CHECK_PSF_NORM 0                       //!< Check PSF normalization

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATPsf::GLATPsf(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename FITS file name.
 *
 * Construct instance by loading the point spread function from FITS file.
 ***************************************************************************/
GLATPsf::GLATPsf(const std::string filename)
{
    // Initialise class members
    init_members();

    // Load energy dispersion from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsf::GLATPsf(const GLATPsf& psf)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(psf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATPsf::~GLATPsf(void)
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
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsf& GLATPsf::operator= (const GLATPsf& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(psf);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return point spread function value
 *
 * @param[in] offset Offset angle (deg).
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * @exception GLATException::invalid_response
 *            Invalid response function version.
 *
 * This operator returns the version dependent value of the point spread
 * function.
 *
 * @todo Implement also PSF versions 2 and 3
 ***************************************************************************/
double GLATPsf::operator() (const double& offset, const double& logE,
                            const double& ctheta)
{
    // Initialise response
    double psf = 0.0;

    // Get point spread function value
    switch (m_version) {
    case 1:
        psf = psf_v1(offset, logE, ctheta);
        break;
    default:
        throw GLATException::invalid_response(G_OP_PSF1, 
              "Unsupported response function version "+str(m_version)+".");
        break;
    }

    // Return point spread function
    return psf;
}


/***********************************************************************//**
 * @brief Return point spread function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt LAT pointing.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GLATPsf::operator() (const GLATInstDir& obsDir, const GSkyDir& srcDir,
                            const GEnergy& srcEng, const GTime& srcTime,
                            const GLATPointing& pnt)
{
    // Initialise response
    double psf = 0.0;

    // Return point spread function
    return psf;
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
void GLATPsf::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATPsf* GLATPsf::clone(void) const
{
    return new GLATPsf(*this);
}


/***********************************************************************//**
 * @brief Load point spread function from FITS file
 *
 * @param[in] filename FITS file.
 ***************************************************************************/
void GLATPsf::load(const std::string filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read point spread function from file
    read(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save point spread function into FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] clobber Overwrite existing file?
 ***************************************************************************/
void GLATPsf::save(const std::string filename, bool clobber)
{
    // Open FITS file
    GFits fits(filename);

    // Write point spread function into file
    write(fits);

    // Close FITS file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read point spread function from FITS file
 *
 * @param[in] fits FITS file.
 *
 * @exception GException::fits_hdu_not_found
 *            Response tables not found in FITS file
 * @exception GException::invalid_response
 *            Invalid response type or unsupported response version found.
 *
 * Reads the version dependent point spread function from the FITS file.
 *
 * @todo Implement also PSF versions 2 and 3
 ***************************************************************************/
void GLATPsf::read(const GFits& fits)
{
    // Clear instance
    clear();

    // Get pointer to PSF parameters table
    GFitsTable* hdu_rpsf  = fits.table("RPSF");
    if (hdu_rpsf == NULL)
        throw GException::fits_hdu_not_found(G_READ, "RPSF");

    // Get pointer to PSF scaling parameters table
    GFitsTable* hdu_scale = fits.table("PSF_SCALING_PARAMS");
    if (hdu_scale == NULL)
        throw GException::fits_hdu_not_found(G_READ, "PSF_SCALING_PARAMS");

    // Determine PSF version
    try {
        m_version = hdu_rpsf->integer("PSFVER");
    }
    catch (GException::fits_key_not_found) {
        m_version = 1;
    }

    // Determine PSF type
    std::string detnam = strip_whitespace(toupper(hdu_rpsf->string("DETNAM")));
    if (detnam == "FRONT")
        m_front = true;
    else if (detnam == "BACK")
        m_front = false;
    else
        throw GLATException::invalid_response(G_READ, 
              "Unknown response type "+detnam+".");

    // Read point spread function
    switch (m_version) {
    case 1:
        read_psf_v1(hdu_rpsf);
        break;
    default:
        throw GLATException::invalid_response(G_READ, 
              "Unsupported response function version "+str(m_version)+".");
        break;
    }

    // Read scaling parameters
    read_scale(hdu_scale);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file
 *
 * @param[in] fits FITS file.
 *
 * @todo Implement also PSF versions 2 and 3
 ***************************************************************************/
void GLATPsf::write(GFits& fits) const
{
    // Write point spread function
    switch (m_version) {
    case 1:
        write_psf_v1(fits);
        break;
    default:
        throw GLATException::invalid_response(G_READ, 
              "Unsupported response function version "+str(m_version)+".");
        break;
    }

    // Write scaling parameters
    write_scale(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum cos(theta) angle for point spread function access
 *
 * @param[in] ctheta Cosine of maximum zenith angle.
 ***************************************************************************/
void GLATPsf::costhetamin(const double& ctheta)
{
    // Set minimum cos(theta) value
    m_min_ctheta = ctheta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point spread function information
 ***************************************************************************/
std::string GLATPsf::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATPsf ===");
    result.append("\n"+parformat("Version")+str(version()));
    result.append("\n"+parformat("Detector section"));
    if (isfront())
        result.append("front");
    else
        result.append("back");
    result.append("\n"+parformat("Energy scaling"));
    result.append("sqrt(");
    result.append("("+str(m_scale_par1)+"*(E/100)^"+str(m_scale_index)+")^2");
    result.append(" + ("+str(m_scale_par2)+")^2)");

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
void GLATPsf::init_members(void)
{
    // Initialise members
    m_version = 1;
    m_front   = true;
    m_rpsf_bins.clear();
    m_scale_par1  = 0.0;
    m_scale_par2  = 0.0;
    m_scale_index = 0.0;
    m_min_ctheta  = 0.0;

    // Initialise PSF version 1 members
    m_ncore.clear();
    m_sigma.clear();
    m_gcore.clear();
    m_gtail.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GLATPsf::copy_members(const GLATPsf& psf)
{
    // Copy members
    m_version     = psf.m_version;
    m_front       = psf.m_front;
    m_rpsf_bins   = psf.m_rpsf_bins;
    m_scale_par1  = psf.m_scale_par1;
    m_scale_par2  = psf.m_scale_par2;
    m_scale_index = psf.m_scale_index;
    m_min_ctheta  = psf.m_min_ctheta;

    // Copy PSF version 1 members
    m_ncore       = psf.m_ncore;
    m_sigma       = psf.m_sigma;
    m_gcore       = psf.m_gcore;
    m_gtail       = psf.m_gtail;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read PSF scale factors from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GException::fits_column_not_found
 *            Table column not found
 ***************************************************************************/
void GLATPsf::read_scale(const GFitsTable* hdu)
{
    // Get pointer to column
    GFitsTableCol* scale = ((GFitsTable*)hdu)->column("PSFSCALE");
    if (scale == NULL)
        throw GException::fits_column_not_found(G_READ_SCALE, "PSFSCALE");

    // Get scaling factors
    if (isfront()) {
        m_scale_par1 = scale->real(0,0);
        m_scale_par2 = scale->real(0,1);
    }
    else {
        m_scale_par1 = scale->real(0,2);
        m_scale_par2 = scale->real(0,3);
    }
    m_scale_index = scale->real(0,4);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write PSF scale factors
 *
 * @param[in] file FITS file.
 ***************************************************************************/
void GLATPsf::write_scale(GFits& file) const
{
    // Create new binary table
    GFitsBinTable* hdu_scale = new GFitsBinTable;

    // Set table attributes
    hdu_scale->extname("PSF_SCALING_PARAMS");

    // Allocate floating point vector column
    GFitsTableFloatCol col_scale = GFitsTableFloatCol("PSFSCALE",  1, 5);

    // Fill columns
    if (isfront()) {
        col_scale(0,0) = m_scale_par1;
        col_scale(0,1) = m_scale_par2;
        col_scale(0,2) = 0.0;
        col_scale(0,3) = 0.0;
    }
    else {
        col_scale(0,0) = 0.0;
        col_scale(0,1) = 0.0;
        col_scale(0,2) = m_scale_par1;
        col_scale(0,3) = m_scale_par2;
    }
    col_scale(0,4) = m_scale_index;

    // Append column to table
    hdu_scale->append_column(col_scale);

    // Append HDU to FITS file
    file.append(*hdu_scale);

    // Free binary table
    delete hdu_scale;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return scale factor for energy (in MeV)
 *
 * @param[in] energy Photon energy (in MeV).
 ***************************************************************************/
double GLATPsf::scale_factor(const double& energy) const
{
    // Compute scale factor
    double f1       = m_scale_par1 * pow(0.01*energy, m_scale_index);
    double scale    = sqrt(f1*f1 + m_scale_par2*m_scale_par2);

    // Return scale factor
    return scale;
}


/*==========================================================================
 =                                                                         =
 =                             PSF Version 1                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Read point spread function from FITS table (version 1)
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GException::fits_column_not_found
 *            Table column not found
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 ***************************************************************************/
void GLATPsf::read_psf_v1(const GFitsTable* hdu)
{
    // Clear arrays
    m_ncore.clear();
    m_sigma.clear();
    m_gcore.clear();
    m_gtail.clear();

    // Get energy and cos theta binning
    m_rpsf_bins.read(hdu);

    // Set minimum cos(theta)
    m_min_ctheta = m_rpsf_bins.costheta_lo(0);

    // Continue only if there are bins
    int size = m_rpsf_bins.size();
    if (size > 0) {

        // Allocate arrays
        m_ncore.reserve(size);
        m_sigma.reserve(size);
        m_gcore.reserve(size);
        m_gtail.reserve(size);

        // Get pointer to columns
        GFitsTableCol* ncore = ((GFitsTable*)hdu)->column("NCORE");
        GFitsTableCol* sigma = ((GFitsTable*)hdu)->column("SIGMA");
        GFitsTableCol* gcore = ((GFitsTable*)hdu)->column("GCORE");
        GFitsTableCol* gtail = ((GFitsTable*)hdu)->column("GTAIL");
        if (ncore == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "NCORE");
        if (sigma == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "SIGMA");
        if (gcore == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "GCORE");
        if (gtail == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "GTAIL");

        // Check consistency of columns
        if (ncore->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       ncore->number(), size);
        if (sigma->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       sigma->number(), size);
        if (gcore->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       gcore->number(), size);
        if (gtail->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       gtail->number(), size);

        // Copy data
        for (int i = 0; i < size; ++i) {
            m_ncore.push_back(ncore->real(0,i));
            m_sigma.push_back(sigma->real(0,i));
            m_gcore.push_back(gcore->real(0,i));
            m_gtail.push_back(gtail->real(0,i));
        }

    } // endif: there were bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file (version 1)
 *
 * @param[in] file FITS file.
 *
 * This method does not write anything if the instance is empty.
 ***************************************************************************/
void GLATPsf::write_psf_v1(GFits& file) const
{
    // Continue only if there are bins
    int size = m_rpsf_bins.size();
    if (size > 0) {

        // Create new binary table
        GFitsBinTable* hdu_rpsf = new GFitsBinTable;

        // Set table attributes
        hdu_rpsf->extname("RPSF");

        // Write boundaries into table
        m_rpsf_bins.write(hdu_rpsf);

        // Allocate floating point vector columns
        GFitsTableFloatCol col_ncore = GFitsTableFloatCol("NCORE",  1, size);
        GFitsTableFloatCol col_sigma = GFitsTableFloatCol("SIGMA",  1, size);
        GFitsTableFloatCol col_gcore = GFitsTableFloatCol("GCORE",  1, size);
        GFitsTableFloatCol col_gtail = GFitsTableFloatCol("GTAIL",  1, size);

        // Fill columns
        for (int i = 0; i < size; ++i) {
            col_ncore(0,i) = m_ncore[i];
            col_sigma(0,i) = m_sigma[i];
            col_gcore(0,i) = m_gcore[i];
            col_gtail(0,i) = m_gtail[i];
        }

        // Append columns to table
        hdu_rpsf->append_column(col_ncore);
        hdu_rpsf->append_column(col_sigma);
        hdu_rpsf->append_column(col_gcore);
        hdu_rpsf->append_column(col_gtail);

        // Append HDU to FITS file
        file.append(*hdu_rpsf);

        // Free binary table
        delete hdu_rpsf;

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return point spread function value (version 1)
 *
 * @param[in] offset Offset angle (deg).
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * @todo Some optimisation could be done as in many cases gcore==gtail,
 *       and for this special case ntail=ncore, hence things become a little
 *       simpler.
 ***************************************************************************/
double GLATPsf::psf_v1(const double& offset, const double& logE,
                       const double& ctheta)
{
    // Set constants
    const double ub = 10.0;

    // Initialise response
    double psf = 0.0;

    // Compute point spread function
    if (ctheta >= m_min_ctheta) {

        // Get response parameters
        double ncore = m_rpsf_bins.interpolate(logE, ctheta, m_ncore);
        double sigma = m_rpsf_bins.interpolate(logE, ctheta, m_sigma);
        double gcore = m_rpsf_bins.interpolate(logE, ctheta, m_gcore);
        double gtail = m_rpsf_bins.interpolate(logE, ctheta, m_gtail);

        // Compute energy in MeV
        double energy = pow(10.0, logE);

        // Rescale the sigma value after interpolation
        sigma *= scale_factor(energy);

        // Compute base function argument
        double r = deg2rad * offset / sigma;
        double u = 0.5 * r * r;

        // Compute normalization of tail
        double ntail = ncore * (base_fct_v1(ub, gcore) / base_fct_v1(ub, gtail));

        // Ensure that PSF integrates to unity. For small energies perform
        // a numerical integration over the solid angle, while for larger
        // energies use a small angle approximation.
        if (energy < 120.0) {
            GLATPsf::base_integrand_v1 integrand(ncore, ntail, sigma, gcore, gtail);
            GIntegral integral(&integrand);
            ncore /= integral.romb(0.0, pihalf) * twopi;
        }
        else {
            double rmax = pihalf / sigma;
            double umax = 0.5 * rmax * rmax;
            double norm = ncore*base_int_v1(umax, gcore) + ntail*base_int_v1(umax, gtail);
            ncore /= norm * twopi * sigma * sigma;
        }

        // Re-compute normalization of tail
        ntail = ncore * (base_fct_v1(ub, gcore) / base_fct_v1(ub, gtail));

        // Compute PSF value
        psf = ncore * base_fct_v1(u, gcore) + ntail * base_fct_v1(u, gtail);

        // Compile option: check PSF normalization
        #if G_CHECK_PSF_NORM
        GLATPsf::base_integrand_v1 integrand(ncore, ntail, sigma, gcore, gtail);
        GIntegral integral(&integrand);
        double sum = integral.romb(0.0, pihalf) * twopi;
        std::cout << "Energy=" << energy;
        std::cout << " Offset=" << offset;
        std::cout << " cos(theta)=" << ctheta;
        std::cout << " error=" << sum-1.0 << std::endl;
        #endif
    }

    // Return point spread function
    return psf;
}


/***********************************************************************//**
 * @brief Return point spread base function value (version 1)
 *
 * @param[in] u Function argument.
 * @param[in] gamma Index.
 *
 * The version 1 PSF base function is given by
 * \f[\left(1 - \frac{1}{\Gamma} \right)
 *    \left(1 + \frac{u}{\Gamma} \right)^{-\Gamma}\f]
 ***************************************************************************/
double GLATPsf::base_fct_v1(const double& u, const double& gamma)
{
    // Get base function value. The special case of gamma==1 is a ugly
    // kluge because of sloppy programming in handoff response when
    // setting boundaries of fit parameters for the PSF.
    double base = (gamma == 1)
                  ? (1.0 - 1.0/1.001) * pow(1.0 + u/1.001, -1.001)
                  : (1.0 - 1.0/gamma) * pow(1.0 + u/gamma, -gamma);

    // Return base function
    return base;
}


/***********************************************************************//**
 * @brief Return approximation of point spread base function integral
 *        (version 1)
 *
 * @param[in] u Function argument.
 * @param[in] gamma Index.
 *
 * The version 1 PSF base function integral is approximated by
 * \f[1 - \left(1 + \frac{u}{\Gamma} \right)^{1-\Gamma}\f]
 * which is valid for small angles \f$u\f$. For larger angles a numerical
 * integration of the base function has to be performed.
 ***************************************************************************/
double GLATPsf::base_int_v1(const double& u, const double& gamma)
{
    // Compute integral of base function
    double integral = 1.0 - pow(1.0 + u/gamma, 1.0 - gamma);

    // Return integral
    return integral;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] psf Point spread function.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATPsf& psf)
{
     // Write point spread function in output stream
    os << psf.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GLATPsf& psf)
{
    // Write point spread function into logger
    log << psf.print();

    // Return logger
    return log;
}
