/***************************************************************************
 *        GLATPsfV3.cpp - Fermi-LAT point spread function version 3        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GLATPsfV3.cpp
 * @brief Fermi-LAT point spread function version 3 class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATPsfV3.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                                 "GLATPsfV3::read(GFitsTable*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_APPROXIMATE_PSF_INTEGRAL        //!< Use approximate PSF integral

/* __ Debug definitions __________________________________________________ */
//#define G_CHECK_PSF_NORM                       //!< Check PSF normalization

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATPsfV3::GLATPsfV3(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsfV3::GLATPsfV3(const GLATPsfV3& psf)
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
GLATPsfV3::~GLATPsfV3(void)
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
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsfV3& GLATPsfV3::operator= (const GLATPsfV3& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Copy base class members
        this->GLATPsfBase::operator=(psf);

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
void GLATPsfV3::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GLATPsfBase::free_members();

    // Initialise members
    this->GLATPsfBase::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATPsfV3* GLATPsfV3::clone(void) const
{
    return new GLATPsfV3(*this);
}


/***********************************************************************//**
 * @brief Read point spread function from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 *
 * Reads point spread function information from FITS HDU. In addition to the
 * energy and costheta binning information, 6 columns are expected:
 * NCORE, NTAIL, SCORE, STAIL, GCORE, and GTAIL.
 *
 * The method assures that NCORE is set properly for each energy and
 * cos(theta) bin so that the integral over the PSF amount to unity. This
 * normalization is done by the method normalize_psf.
 ***************************************************************************/
void GLATPsfV3::read(const GFitsTable* hdu)
{
    // Clear arrays
    m_ncore.clear();
    m_ntail.clear();
    m_score.clear();
    m_stail.clear();
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
        m_ntail.reserve(size);
        m_score.reserve(size);
        m_stail.reserve(size);
        m_gcore.reserve(size);
        m_gtail.reserve(size);

        // Get pointer to columns
        const GFitsTableCol* ncore = &(*hdu)["NCORE"];
        const GFitsTableCol* ntail = &(*hdu)["NTAIL"];
        const GFitsTableCol* score = &(*hdu)["SCORE"];
        const GFitsTableCol* stail = &(*hdu)["STAIL"];
        const GFitsTableCol* gcore = &(*hdu)["GCORE"];
        const GFitsTableCol* gtail = &(*hdu)["GTAIL"];

        // Check consistency of columns
        if (ncore->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       ncore->number(), size);
        }
        if (ntail->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       ntail->number(), size);
        }
        if (score->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       score->number(), size);
        }
        if (stail->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       stail->number(), size);
        }
        if (gcore->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       gcore->number(), size);
        }
        if (gtail->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       gtail->number(), size);
        }

        // Copy data
        for (int i = 0; i < size; ++i) {
            m_ncore.push_back(ncore->real(0,i));
            m_ntail.push_back(ntail->real(0,i));
            m_score.push_back(score->real(0,i));
            m_stail.push_back(stail->real(0,i));
            m_gcore.push_back(gcore->real(0,i));
            m_gtail.push_back(gtail->real(0,i));
        }

        // Normalize PSF for all parameters
        normalize_psf();

    } // endif: there were bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file
 *
 * @param[in] file FITS file.
 *
 * Writes the PSF into the extension "RPSF" of a FITS file. This method
 * does not check if a "RPSF" extension exists so far, it simply adds one
 * each time it is called.
 *
 * Nothing is done if the PSF size is 0.
 *
 * @todo Check if a RPSF extension exists already in FITS file
 ***************************************************************************/
void GLATPsfV3::write(GFits& file) const
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
        GFitsTableFloatCol col_ntail = GFitsTableFloatCol("NTAIL",  1, size);
        GFitsTableFloatCol col_score = GFitsTableFloatCol("SCORE",  1, size);
        GFitsTableFloatCol col_stail = GFitsTableFloatCol("STAIL",  1, size);
        GFitsTableFloatCol col_gcore = GFitsTableFloatCol("GCORE",  1, size);
        GFitsTableFloatCol col_gtail = GFitsTableFloatCol("GTAIL",  1, size);

        // Fill columns
        for (int i = 0; i < size; ++i) {
            col_ncore(0,i) = m_ncore[i];
            col_ntail(0,i) = m_ntail[i];
            col_score(0,i) = m_score[i];
            col_stail(0,i) = m_stail[i];
            col_gcore(0,i) = m_gcore[i];
            col_gtail(0,i) = m_gtail[i];
        }

        // Append columns to table
        hdu_rpsf->append_column(col_ncore);
        hdu_rpsf->append_column(col_ntail);
        hdu_rpsf->append_column(col_score);
        hdu_rpsf->append_column(col_stail);
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
 * @brief Return point spread function value
 *
 * @param[in] offset Offset angle (deg).
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * Evaluates point spread function by doing a bi-linear interpolation of
 * PSF values obtained at the 4 corners that bound the specified energy and
 * cos(theta) value.
 ***************************************************************************/
double GLATPsfV3::psf(const double& offset, const double& logE,
                      const double& ctheta)
{
    // Initialise response
    double psf = 0.0;

    // Compute point spread function
    if (ctheta >= m_min_ctheta) {

        // Set interpolation indices and weights
        m_rpsf_bins.set(logE, ctheta);

        // Recover information for interpolation
        std::vector<int>    index  = m_rpsf_bins.indices();
        std::vector<double> energy = m_rpsf_bins.energies();
        std::vector<double> weight = m_rpsf_bins.weights();

        // Compute offset angle in radians
        double offset_rad = offset * deg2rad;
        
        // Compute PSF values for the four corners
        double psf0 = eval_psf(offset_rad, energy[0], index[0]);
        double psf1 = eval_psf(offset_rad, energy[1], index[1]);
        double psf2 = eval_psf(offset_rad, energy[2], index[2]);
        double psf3 = eval_psf(offset_rad, energy[3], index[3]);

        // Perform bi-linear interpolation
        psf = weight[0] * psf0 + weight[1] * psf1 +
              weight[2] * psf2 + weight[3] * psf3;

    } // endif: cos(theta) was in valid range

    // Return point spread function
    return psf;
}


/***********************************************************************//**
 * @brief Print point spread function
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GLATPsfV3::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATPsfV3 ===");

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
void GLATPsfV3::init_members(void)
{
    // Initialise members
    m_ncore.clear();
    m_ntail.clear();
    m_score.clear();
    m_stail.clear();
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
void GLATPsfV3::copy_members(const GLATPsfV3& psf)
{
    // Copy members
    m_ncore = psf.m_ncore;
    m_ntail = psf.m_ntail;
    m_score = psf.m_score;
    m_stail = psf.m_stail;
    m_gcore = psf.m_gcore;
    m_gtail = psf.m_gtail;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsfV3::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return point spread base function value
 *
 * @param[in] u Function argument.
 * @param[in] gamma Index.
 *
 * The version 3 PSF base function is given by
 * \f[\left(1 - \frac{1}{\Gamma} \right)
 *    \left(1 + \frac{u}{\Gamma} \right)^{-\Gamma}\f]
 ***************************************************************************/
double GLATPsfV3::base_fct(const double& u, const double& gamma)
{
    // Get base function value. The special case of gamma==1 is a ugly
    // kluge because of sloppy programming in handoff response when
    // setting boundaries of fit parameters for the PSF.
    double base = (gamma == 1)
                  ? (1.0 - 1.0/1.001) * std::pow(1.0 + u/1.001, -1.001)
                  : (1.0 - 1.0/gamma) * std::pow(1.0 + u/gamma, -gamma);

    // Return base function
    return base;
}


/***********************************************************************//**
 * @brief Return approximation of point spread base function integral
 *
 * @param[in] u Function argument.
 * @param[in] gamma Index.
 *
 * The version 3 PSF base function integral is approximated by
 * \f[1 - \left(1 + \frac{u}{\Gamma} \right)^{1-\Gamma}\f]
 * which is valid for small angles \f$u\f$. For larger angles a numerical
 * integration of the base function has to be performed.
 *
 * @todo Verify that 1+u/gamma is not negative
 ***************************************************************************/
double GLATPsfV3::base_int(const double& u, const double& gamma)
{
    // Compute integral of base function
    double integral = 1.0 - std::pow(1.0 + u/gamma, 1.0 - gamma);

    // Return integral
    return integral;
}


/***********************************************************************//**
 * @brief Evaluate PSF for a specific set of parameters
 *
 * @param[in] offset Offset angle (radians).
 * @param[in] energy Energy (MeV).
 * @param[in] index Parameter array index.
 *
 * Evaluates PSF as function of offset angle and energy for a specific set
 * of PSF parameters. The parameter set that is used is specified by the
 * index parameter. The energy parameter only serves to scale the score and
 * stail parameters of the PSF.
 ***************************************************************************/
double GLATPsfV3::eval_psf(const double& offset, const double& energy,
                           const int& index)
{
    // Get energy scaling
    double scale = scale_factor(energy);

    // Get parameters
    double ncore(m_ncore[index]);
    double ntail(m_ntail[index]);
    double score(m_score[index] * scale);
    double stail(m_stail[index] * scale);
    double gcore(m_gcore[index]);
    double gtail(m_gtail[index]);

    // Compute argument
    double rc = offset / score;
    double uc = 0.5 * rc * rc;
    double rt = offset / stail;
    double ut = 0.5 * rt * rt;

    // Evaluate PSF
    double psf = ncore * (base_fct(uc, gcore) + ntail * base_fct(ut, gtail));
    
    // Return PSF
    return psf;
}


/***********************************************************************//**
 * @brief Integrates PSF for a specific set of parameters
 *
 * @param[in] energy Energy (MeV).
 * @param[in] index Parameter array index.
 *
 * Integrates PSF for a specific set of parameters.
 *
 * Compile option G_APPROXIMATE_PSF_INTEGRAL:
 * If defined, a numerical PSF integral is only performed for energies
 * < 120 MeV, while for larger energies the small angle approximation is
 * used. In not defined, a numerical PSF integral is performed for all
 * energies.
 * This option is kept for comparison with the Fermi/LAT ScienceTools who
 * select the integration method based on the true photon energy. As the
 * normalization is only performed once upon loading of the PSF, CPU time
 * is not really an issue here, and we can afford the more precise numerical
 * integration. Note that the uncertainties of the approximation at energies
 * near to 120 MeV reaches 0.1%.
 *
 * @todo Implement gcore and gtail checking
 ***************************************************************************/
double GLATPsfV3::integrate_psf(const double& energy, const int& index)
{
    // Initialise integral
    double psf = 0.0;

    // Get energy scaling
    double scale = scale_factor(energy);

    // Get parameters
    double ncore(m_ncore[index]);
    double ntail(m_ntail[index]);
    double score(m_score[index] * scale);
    double stail(m_stail[index] * scale);
    double gcore(m_gcore[index]);
    double gtail(m_gtail[index]);

    // Make sure that gcore and gtail are not negative
    //if (gcore < 0 || gtail < 0) {
    //}

    // Do we need an exact integral?
    #if defined(G_APPROXIMATE_PSF_INTEGRAL)
    if (energy < 120) {
    #endif

        // Allocate integrand
        GLATPsfV3::base_integrand integrand(ncore, ntail, score, stail, gcore, gtail);

        // Allocate integral
        GIntegral integral(&integrand);

        // Integrate radially from 0 to 90 degrees
        psf = integral.romb(0.0, pihalf) * twopi;
    
    #if defined(G_APPROXIMATE_PSF_INTEGRAL)
    } // endif: exact integral was performed

    // No, so we use the small angle approximation
    else {

        // Compute arguments
        double rc = pihalf / score;
        double uc = 0.5 * rc * rc;
        double sc = twopi * score * score;
        double rt = pihalf / stail;
        double ut = 0.5 * rt * rt;
        double st = twopi * stail * stail;

        // Evaluate PSF integral (from 0 to 90 degrees)
        psf = ncore * (base_int(uc, gcore) * sc +
                       base_int(ut, gtail) * st * ntail);
    
    }
    #endif

    // Return PSF integral
    return psf;
}


/***********************************************************************//**
 * @brief Normalize PSF for all parameters
 *
 * Makes sure that PSF is normalized for all parameters. We assure this by
 * looping over all parameter nodes, integrating the PSF for each set of
 * parameters, and dividing the NCORE parameter by the integral.
 *
 * Compile option G_CHECK_PSF_NORM:
 * If defined, checks that the PSF is normalized correctly.
 ***************************************************************************/
void GLATPsfV3::normalize_psf(void)
{
    // Loop over all energy bins
    for (int ie = 0; ie < m_rpsf_bins.nenergies(); ++ie) {

        // Extract energy value (in MeV)
        double energy = m_rpsf_bins.energy(ie);

        // Loop over all cos(theta) bins
        for (int ic = 0; ic < m_rpsf_bins.ncostheta(); ++ic) {

            // Get parameter index
            int index = m_rpsf_bins.index(ie, ic);

            // Integrate PSF
            double norm = integrate_psf(energy, index);

            // Normalize PSF
            m_ncore[index] /= norm;

            // Compile option: check PSF normalization
            #if defined(G_CHECK_PSF_NORM)
            double scale = scale_factor(energy);
            double ncore(m_ncore[index]);
            double ntail(m_ntail[index]);
            double score(m_score[index] * scale);
            double stail(m_stail[index] * scale);
            double gcore(m_gcore[index]);
            double gtail(m_gtail[index]);
            GLATPsfV3::base_integrand integrand(ncore, ntail, score, stail, gcore, gtail);
            GIntegral integral(&integrand);
            double sum = integral.romb(0.0, pihalf) * twopi;
            std::cout << "Energy=" << energy;
            std::cout << " cos(theta)=" << m_rpsf_bins.costheta_lo(ic);
            std::cout << " error=" << sum-1.0 << std::endl;
            #endif

        } // endfor: looped over cos(theta)

    } // endfor: looped over energies
    
    // Return
    return;
}
