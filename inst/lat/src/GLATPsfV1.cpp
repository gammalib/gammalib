/***************************************************************************
 *       GLATPsfV1.cpp  -  Fermi LAT point spread function version 1       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GLATPsfV1.cpp
 * @brief Fermi LAT point spread function version 1 class implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATPsfV1.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                           "GLATPsfV1::read(GFitsTable*)"

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
GLATPsfV1::GLATPsfV1(void)
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
GLATPsfV1::GLATPsfV1(const GLATPsfV1& psf)
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
GLATPsfV1::~GLATPsfV1(void)
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
GLATPsfV1& GLATPsfV1::operator= (const GLATPsfV1& psf)
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
void GLATPsfV1::clear(void)
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
GLATPsfV1* GLATPsfV1::clone(void) const
{
    return new GLATPsfV1(*this);
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
 * energy and costheta binning information, 4 columns are expected:
 * NCORE, SIGMA, GCORE, and GTAIL.
 ***************************************************************************/
void GLATPsfV1::read(const GFitsTable* hdu)
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
        const GFitsTableCol* ncore = &(*hdu)["NCORE"];
        const GFitsTableCol* sigma = &(*hdu)["SIGMA"];
        const GFitsTableCol* gcore = &(*hdu)["GCORE"];
        const GFitsTableCol* gtail = &(*hdu)["GTAIL"];

        // Check consistency of columns
        if (ncore->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       ncore->number(), size);
        }
        if (sigma->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       sigma->number(), size);
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
            m_sigma.push_back(sigma->real(0,i));
            m_gcore.push_back(gcore->real(0,i));
            m_gtail.push_back(gtail->real(0,i));
        }

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
void GLATPsfV1::write(GFits& file) const
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
 * @brief Return point spread function value
 *
 * @param[in] offset Offset angle (deg).
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * @todo Some optimisation could be done as in many cases gcore==gtail,
 *       and for this special case ntail=ncore, hence things become a little
 *       simpler.
 ***************************************************************************/
double GLATPsfV1::psf(const double& offset, const double& logE,
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
        double ntail = ncore * (base_fct(ub, gcore) / base_fct(ub, gtail));

        // Ensure that PSF integrates to unity. For small energies perform
        // a numerical integration over the solid angle, while for larger
        // energies use a small angle approximation.
        if (energy < 120.0) {
            GLATPsfV1::base_integrand integrand(ncore, ntail, sigma, gcore, gtail);
            GIntegral integral(&integrand);
            ncore /= integral.romb(0.0, pihalf) * twopi;
        }
        else {
            double rmax = pihalf / sigma;
            double umax = 0.5 * rmax * rmax;
            double norm = ncore*base_int(umax, gcore) + ntail*base_int(umax, gtail);
            ncore /= norm * twopi * sigma * sigma;
        }

        // Re-compute normalization of tail
        ntail = ncore * (base_fct(ub, gcore) / base_fct(ub, gtail));

        // Compute PSF value
        psf = ncore * base_fct(u, gcore) + ntail * base_fct(u, gtail);

        // Compile option: check PSF normalization
        #if G_CHECK_PSF_NORM
        GLATPsfV1::base_integrand integrand(ncore, ntail, sigma, gcore, gtail);
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


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATPsfV1::init_members(void)
{
    // Initialise members
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
void GLATPsfV1::copy_members(const GLATPsfV1& psf)
{
    // Copy members
    m_ncore = psf.m_ncore;
    m_sigma = psf.m_sigma;
    m_gcore = psf.m_gcore;
    m_gtail = psf.m_gtail;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsfV1::free_members(void)
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
 * The version 1 PSF base function is given by
 * \f[\left(1 - \frac{1}{\Gamma} \right)
 *    \left(1 + \frac{u}{\Gamma} \right)^{-\Gamma}\f]
 ***************************************************************************/
double GLATPsfV1::base_fct(const double& u, const double& gamma)
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
 *
 * @param[in] u Function argument.
 * @param[in] gamma Index.
 *
 * The version 1 PSF base function integral is approximated by
 * \f[1 - \left(1 + \frac{u}{\Gamma} \right)^{1-\Gamma}\f]
 * which is valid for small angles \f$u\f$. For larger angles a numerical
 * integration of the base function has to be performed.
 ***************************************************************************/
double GLATPsfV1::base_int(const double& u, const double& gamma)
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
