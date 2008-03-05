/***************************************************************************
 *               GLATResponse.cpp  -  GLAST LAT Response class             *
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
#define G_LOAD            "load(const std::string&, const std::string&)"
#define G_INIT_AEFF       "GLATResponse::init_aeff()"
#define G_INIT_PSF        "GLATResponse::init_psf()"
#define G_INIT_EDISP      "GLATResponse::init_edisp()"
#define G_GET_FITS_VECTOR "GLATResponse::get_fits_vector(const GFitsHDU*, const std::string&, int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const double scale_front_c0 = 3.77e-4;
const double scale_front_c1 = 5.8e-2;
const double scale_back_c0  = 1.3e-3;
const double scale_back_c1  = 9.6e-2;

const double pi          = 3.1415926535897931159979635;
const double pihalf      = 1.5707963267948966;
const double twopi       = 6.2831853071795862319959269;
const double sqrt2pi     = 2.5066282746310002416123552;
const double twosqrt2ln2 = 2.3548200450309493270140138;
const double deg2rad     = 0.0174532925199432954743717;
const double rad2deg     = 57.295779513082322864647722;
const double ln10        = 2.3025850929940459010936138;


/*==========================================================================
 =                                                                         =
 =                   GLATResponse constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GLATResponse::GLATResponse() : GResponse()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param rsp Response to be copied
 ***************************************************************************/
GLATResponse::GLATResponse(const GLATResponse& rsp) : GResponse(rsp)
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
GLATResponse::~GLATResponse()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GLATResponse operators                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param rsp Response to be assigned
 ***************************************************************************/
GLATResponse& GLATResponse::operator= (const GLATResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

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
 =                        GLATResponse public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of instrument response function.
 ***************************************************************************/
double GLATResponse::irf(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time)
{
    // Return IRF value
    return 0.0;
}


/***********************************************************************//**
 * @brief Return effective area
 ***************************************************************************/
double GLATResponse::aeff(const GSkyDir& obsDir, const double& obsEng,
                          const GSkyDir& srcDir, const double& srcEng,
                          const GSkyDir& instPntDir, const double& instPosAng,
                          const double& time)
{
    // Return Aeff value
    return 0.0;
}


/***********************************************************************//**
 * @brief Return value of point spread function
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
 * @brief Return value of energy dispersion
 ***************************************************************************/
double GLATResponse::edisp(const GSkyDir& obsDir, const double& obsEng,
                           const GSkyDir& srcDir, const double& srcEng,
                           const GSkyDir& instPntDir, const double& instPosAng,
                           const double& time)
{
    // Return Edisp value
    return 0.0;
}


/***********************************************************************//**
 * @brief Set the path to the calibration database.
 *
 * @param caldb Absolute path to calibration database
 *
 * NOTE: So far no check is done on whether the path exists!
 ***************************************************************************/
void GLATResponse::set_caldb(const std::string& caldb)
{
    // Simply copy path
    // NOTE: Some check should be done on whether the path exists !!!
    m_caldb = caldb;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load a specified LAT response function.
 *
 * @param rspname Name of response
 * @param rsptype Type of response ('front' or 'back')
 *
 * Loads the specified GLAST LAT response from the calibration database and
 * performs some pre-calculations for faster response determination.
 ***************************************************************************/
void GLATResponse::load(const std::string& rspname, const std::string& rsptype)
{
    // Store response name
    m_rspname = rspname;

    // Convert response type string to section
    if (rsptype == "front")
        m_section = 0;
    else if (rsptype == "back")
        m_section = 1;
    else
        throw GException::rsp_invalid_type(G_LOAD, rsptype);

    // Store response type
    m_rsptype = rsptype;

    // Initialise effective area
    aeff_init();

    // Initialise PSF
    psf_init();

    // Initialise energy dispersion
    edisp_init();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save response
 *
 * @param rspname Filename into which the response will be saved
 ***************************************************************************/
void GLATResponse::save(const std::string& rspname) const
{
    // Open FITS file
    GFits fits;
    fits.open(rspname);

    // Build PSF tables and images
    GFitsBinTable table_psf_bounds(1);
    int           naxes_psf2[] = {m_psf_energy_num, m_psf_ctheta_num};
    int           naxes_psf3[] = {m_psf_angle_num, m_psf_energy_num,
                                  m_psf_ctheta_num};
    GFitsDblImage image_psf_psf(3, naxes_psf3, m_psf);
    GFitsDblImage image_psf_norm(2, naxes_psf2, m_norm);
    GFitsDblImage image_psf_sigma(2, naxes_psf2, m_sigma);

    // Construct PSF HDUs
    GFitsHDU hdu_psf_bounds(table_psf_bounds);
    GFitsHDU hdu_psf_psf(image_psf_psf);
    GFitsHDU hdu_psf_norm(image_psf_norm);
    GFitsHDU hdu_psf_sigma(image_psf_sigma);
    hdu_psf_bounds.extname("PBOUNDS");
    hdu_psf_psf.extname("PSF");
    hdu_psf_norm.extname("NORM");
    hdu_psf_sigma.extname("SIGMA");

    // Append HDUs to FITS file
    fits.append_hdu(hdu_psf_bounds);
    fits.append_hdu(hdu_psf_psf);
    fits.append_hdu(hdu_psf_norm);
    fits.append_hdu(hdu_psf_sigma);

    // Save FITS file
    fits.save();

    // Return
    return;
}
/*
    GVector v_energy_lo = get_fits_vector(hdu, "ENERG_LO");
    GVector v_energy_hi = get_fits_vector(hdu, "ENERG_HI");
    GVector v_ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector v_ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
*/

/*==========================================================================
 =                                                                         =
 =               GLATResponse private effective area methods               =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        Initialise effective area                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::aeff_init(void)
{
    // Build filename
    std::string filename = "aeff_"  + m_rspname + "_" + m_rsptype + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer to effective area HDU
    GFitsHDU* hdu = file.hdu("EFFECTIVE AREA");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_AEFF, "EFFECTIVE AREA");

    // Get the data
    GVector energ_lo  = get_fits_vector(hdu, "ENERG_LO");
    GVector energ_hi  = get_fits_vector(hdu, "ENERG_HI");
    GVector ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
    GVector effarea   = get_fits_vector(hdu, "EFFAREA");

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                      GLATResponse private PSF methods                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Initialise PSF                              *
 * ----------------------------------------------------------------------- *
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

    // Get the data
    GVector v_energy_lo = get_fits_vector(hdu, "ENERG_LO");
    GVector v_energy_hi = get_fits_vector(hdu, "ENERG_HI");
    GVector v_ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector v_ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
    GVector v_ncore     = get_fits_vector(hdu, "NCORE");
    GVector v_sigma     = get_fits_vector(hdu, "SIGMA");
    GVector v_gcore     = get_fits_vector(hdu, "GCORE");
    GVector v_gtail     = get_fits_vector(hdu, "GTAIL");

    // Set binning in angle, energy and cos theta
    m_psf_angle_num  = 2000;
    m_psf_angle_min  = 0.01;
    m_psf_angle_bin  = 0.01;
    m_psf_energy_num = v_energy_lo.size();
    m_psf_ctheta_num = v_ctheta_lo.size();

    // Free old PSF memory
    if (m_psf   != NULL) delete [] m_psf;
    if (m_norm  != NULL) delete [] m_norm;
    if (m_sigma != NULL) delete [] m_sigma;

    // Allocate memory for PSF
    int psf_size = m_psf_angle_num * m_psf_energy_num * m_psf_ctheta_num;
    int bin_size = m_psf_energy_num * m_psf_ctheta_num;
    m_psf   = new double[psf_size];
    m_norm  = new double[bin_size];
    m_sigma = new double[bin_size];

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
    for (int ie = 0; ie < m_psf_energy_num; ++ie) {

        // Get mean energy of bin
        // NOTE: This should be replaced by a logarithmic mean
        double energy = 0.5 * (v_energy_lo(ie) + v_energy_hi(ie));

        // Get scale at specified energy
        double scale = psf_scale(energy);

        // Loop over all cos theta bins
        for (int ic = 0; ic < m_psf_ctheta_num; ++ic) {

            // Extract PSF parameters
            int    inx   = ic * m_psf_energy_num + ie;
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


/***************************************************************************
 *                              Get PSF value                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/// Returns the GLAST LAT PSF value
/// @param delta angular separation between true and observed photon direction
/// @param logE logarithm of the true photon energy
/// @param ctheta cosine of the zenith angle with respect to the pointing axis
double GLATResponse::psf_get(const double& delta, const double& logE, const double& ctheta)
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


/***************************************************************************
 *                            Returns PSF scale                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/// PSF size scale as
///   scale = sqrt(c0^2 + (c1*t)^2)
/// with
///   t = (0.01*E)^-0.8
/// where E is the energy in MeV
/// The coefficients c0 and c1 dependend on the section:
/// Front: c0 = 3.77e-4, c1 = 5.8e-2
/// Back:  c0 = 1.3e-3,  c1 = 9.6e-2
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


/***************************************************************************
 *                         Returns PSF basis function                      *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/// This method returns the PSF basis function for a single value.
/// The PSF basis functions is defined by
///   psf(u,gamma) = (1.0 - 1.0/gamma) * (1.0 + u/gamma)^(-gamma)
/// where
///   u     is an angle^2 (typically the scaled angular distance squared)
///   gamma is a powerlaw index
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
/// This method returns the PSF basis function for a vector.
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

/*==========================================================================
 =                                                                         =
 =             GLATResponse private energy dispersion methods              =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                      Initialise energy dispersion                       *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::edisp_init(void)
{
    // Build filename
    std::string filename = "edisp_"  + m_rspname + "_" + m_rsptype + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer towards energy dispersion HDU
    GFitsHDU* hdu = file.hdu("ENERGY DISPERSION");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_EDISP, "ENERGY DISPERSION");

    // Get the data
    GVector energ_lo  = get_fits_vector(hdu, "ENERG_LO");
    GVector energ_hi  = get_fits_vector(hdu, "ENERG_HI");
    GVector ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
    GVector dnorm     = get_fits_vector(hdu, "DNORM");
    GVector ltail     = get_fits_vector(hdu, "LTAIL");
    GVector rwidth    = get_fits_vector(hdu, "RWIDTH");
    GVector nr2       = get_fits_vector(hdu, "NR2");
    GVector lt2       = get_fits_vector(hdu, "LT2");
    GVector rt2       = get_fits_vector(hdu, "RT2");

    // ...

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        GLATResponse private methods                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::init_members(void)
{
    // Initialise members
    //m_aeff_name.clear();
    //m_psf_name.clear();
    //m_edisp_name.clear();

    // Initialise PSF members
    m_psf_angle_num  = 0;
    m_psf_angle_min  = 0.0;
    m_psf_angle_max  = 0.0;
    m_psf_angle_bin  = 0.0;
    m_psf_energy_num = 0;
    m_psf_ctheta_num = 0;
    m_psf            = NULL;
    m_norm           = NULL;
    m_sigma          = NULL;

    // By default use HANDOFF response database.
    char* handoff = getenv("HANDOFF_IRF_DIR");
    if (handoff != NULL)
        m_caldb.assign(handoff);

    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::copy_members(const GLATResponse& rsp)
{
    // Copy attributes
    //m_aeff_name  = rsp.m_aeff_name;
    //m_psf_name   = rsp.m_psf_name;
    //m_edisp_name = rsp.m_edisp_name;

    // Copy other membres

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::free_members(void)
{
    // Free PSF memory
    if (m_psf   != NULL) delete [] m_psf;
    if (m_norm  != NULL) delete [] m_norm;
    if (m_sigma != NULL) delete [] m_sigma;

    // Signal that PSF memory is free
    m_psf   = NULL;
    m_norm  = NULL;
    m_sigma = NULL;

    // Return
    return;
}


/***************************************************************************
 *          Load floating point data from vector column into vector        *
 * ----------------------------------------------------------------------- *
 * All information is stored in the first table row using vector columns.  *
 ***************************************************************************/
GVector GLATResponse::get_fits_vector(const GFitsHDU* hdu, const std::string& colname, int row)
{
    // Get pointer to
    GFitsTableCol* ptr = hdu->column(colname);
    if (ptr == NULL)
        throw GException::fits_column_not_found(G_GET_FITS_VECTOR, colname);

    // Determine number of entries
    int num = ptr->number();

    // Load data into vector
    GVector data(num);
    for (int i = 0; i < num; ++i)
        data(i) = ptr->real(row,i);
    //cout << colname << ": " << num << endl;

    // Return vector
    return data;
}


/*==========================================================================
 =                                                                         =
 =                           GLATResponse friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GLATResponse                  =
 =                                                                         =
 ==========================================================================*/
