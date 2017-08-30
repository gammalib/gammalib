/***************************************************************************
 *         GCOMIaq.cpp - COMPTEL instrument response representation        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMIaq.cpp
 * @brief COMPTEL instrument response representation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFilename.hpp"
#include "GIntegral.hpp"
#include "GCOMIaq.hpp"
#include "GCOMSupport.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_RESPONSE_KERNEL_LIMITS
//#define G_RESPSC_NO_WARNINGS

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_IAQWEI
//#define G_DEBUG_RESPSC
//#define G_DEBUG_RESPONSE_KERNEL

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty COMPTEL instrument response representation.
 ***************************************************************************/
GCOMIaq::GCOMIaq(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] iaq COMPTEL instrument response representation.
 **************************************************************************/
GCOMIaq::GCOMIaq(const GCOMIaq& iaq)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(iaq);

    // Return
    return;
}


/***********************************************************************//**
 * @brief COMPTEL instrument response representation constructor
 *
 * @param[in] phigeo_max Maximum geometrical scatter angle (deg).
 * @param[in] phigeo_bin_size Bin size in geometrical scatter angle (deg).
 * @param[in] phibar_max Maximum Compton scatter angle (deg).
 * @param[in] phibar_bin_size Bin size in Compton scatter angle (deg).
 * @param[in] ebounds Boundaries of observed energy.
 *
 * @todo Test input argument validity.
 **************************************************************************/
GCOMIaq::GCOMIaq(const double&   phigeo_max, const double& phigeo_bin_size,
                 const double&   phibar_max, const double& phibar_bin_size,
                 const GEbounds& ebounds)
{
    // Initialise members
    init_members();

    // Compute the number of bins in phigeo and phibar
    int nphigeo = int(phigeo_max / phigeo_bin_size + 0.5);
    int nphibar = int(phibar_max / phibar_bin_size + 0.5);

    // Store the bin definition. Use the number of bins to re-derive the
    // maximum angles
    m_phigeo_bin_size = phigeo_bin_size;
    m_phibar_bin_size = phibar_bin_size;
    m_phigeo_max      = phigeo_bin_size * double(nphigeo);
    m_phibar_max      = phibar_bin_size * double(nphibar);

    // Store the energy boundaries
    m_ebounds = ebounds;

    // Allocate FITS image
    m_iaq = GFitsImageFloat(nphigeo, nphibar);

    // Set FITS image keywords
    m_iaq.card("CTYPE1", "nphig", "Number of geometrical scatter angle pixels");
    m_iaq.card("CRVAL1", 0.5*m_phigeo_bin_size, "[deg] Geometrical scatter angle for reference bin");
    m_iaq.card("CRPIX1", 1.0, "Reference bin of geometrical scatter angle");
    m_iaq.card("CDELT1", m_phigeo_bin_size, "[deg] Geometrical scatter angle bin size");
    m_iaq.card("CTYPE2", "nphib", "Number of Compton scatter angle bins");
    m_iaq.card("CRVAL2", 0.5*m_phibar_bin_size, "[deg] Compton scatter angle for reference bin");
    m_iaq.card("CRPIX2", 1.0, "Reference bin of Compton scatter angle");
    m_iaq.card("CDELT2", m_phibar_bin_size, "[deg] Compton scatter angle bin size");
    m_iaq.card("ENERGIE", 0.0, "[MeV] Source photon energy");
    m_iaq.card("EMIN", m_ebounds.emin().MeV(), "[MeV] Minimum measured photon energy");
    m_iaq.card("EMAX", m_ebounds.emax().MeV(), "[MeV] Maximum measured photon energy");
    m_iaq.card("TELESCOP", "GRO", "Name of telescope");
    m_iaq.card("INSTRUME", "COMPTEL", "Name of instrument");
    m_iaq.card("ORIGIN", "GammaLib", "Origin of FITS file");
    m_iaq.card("OBSERVER", "Unknown", "Observer that created FITS file");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of COMPTEL response object.
 ***************************************************************************/
GCOMIaq::~GCOMIaq(void)
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
 * @param[in] iaq COMPTEL instrument response representation.
 * @return COMPTEL instrument response representation.
 ***************************************************************************/
GCOMIaq& GCOMIaq::operator=(const GCOMIaq& iaq)
{
    // Execute only if object is not identical
    if (this != &iaq) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(iaq);

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
 * Clears COMPTEL instrument response representation by resetting all members
 * to an initial state. Any information that was present in the object before
 * will be lost.
 ***************************************************************************/
void GCOMIaq::clear(void)
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
 *
 * @return Pointer to deep copy of COMPTEL instrument response
 *         representation.
 ***************************************************************************/
GCOMIaq* GCOMIaq::clone(void) const
{
    return new GCOMIaq(*this);
}


/***********************************************************************//**
 * @brief Save COMPTEL instrument response representation into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves the COMPTEL instrument response representation into a FITS file. If
 * the file exists already and the @p clobber parameter is true, the method
 * will overwrite the content of the existing file. Otherwise, an exception
 * is thrown.
 ***************************************************************************/
void GCOMIaq::save(const GFilename& filename, const bool& clobber) const
{
    // Work on a copy of the image to set the extension name
    GFitsImageFloat iaq = m_iaq;

    // Set extension name
    if (filename.has_extname()) {
        iaq.extname(filename.extname());
    }

    // Initialise empty FITS file
    GFits fits;

    // Append IAQ to FITS file
    fits.append(iaq);

    // Save FITS file (without extension name which was extracted earlier
    // and set in the HDU)
    fits.saveto(filename.url(), clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL instrument response representation information
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL instrument response representation
 *         information.
 ***************************************************************************/
std::string GCOMIaq::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMIaq ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of Phigeo bins"));
        if (m_iaq.naxis() > 0) {
            result.append(gammalib::str(m_iaq.naxes(0)));
        }
        else {
            result.append("0");
        }
        result.append("\n"+gammalib::parformat("Number of Phibar bins"));
        if (m_iaq.naxis() > 1) {
            result.append(gammalib::str(m_iaq.naxes(1)));
        }
        else {
            result.append("0");
        }
        result.append("\n"+gammalib::parformat("Maximum Phigeo value"));
        result.append(gammalib::str(m_phigeo_max)+" deg");
        result.append("\n"+gammalib::parformat("Maximum Phibar value"));
        result.append(gammalib::str(m_phibar_max)+" deg");
        result.append("\n"+gammalib::parformat("Phigeo bin size"));
        result.append(gammalib::str(m_phigeo_bin_size)+" deg");
        result.append("\n"+gammalib::parformat("Phibar bin size"));
        result.append(gammalib::str(m_phibar_bin_size)+" deg");
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(m_ebounds.emin().print()+" - ");
        result.append(m_ebounds.emax().print());

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCOMIaq::init_members(void)
{
    // Initialise members
    m_iaq.clear();
    m_ebounds.clear();
    m_response_d1.clear();
    m_response_d2.clear();
    m_ict.clear();
    m_phigeo_max      =  0.0;
    m_phibar_max      =  0.0;
    m_phigeo_bin_size =  0.0;
    m_phibar_bin_size =  0.0;
    m_e1min           =  0.070; //!< Default: 70 keV
    m_e1max           = 20.0;   //!< Default: 20 MeV
    m_e2min           =  0.650; //!< Default: 650 keV
    m_e2max           = 30.0;   //!< Default: 30 MeV

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] iaq COMPTEL instrument response representation.
 ***************************************************************************/
void GCOMIaq::copy_members(const GCOMIaq& iaq)
{
    // Copy attributes
    m_iaq             = iaq.m_iaq;
    m_ebounds         = iaq.m_ebounds;
    m_response_d1     = iaq.m_response_d1;
    m_response_d2     = iaq.m_response_d2;
    m_ict             = iaq.m_ict;
    m_phigeo_max      = iaq.m_phigeo_max;
    m_phibar_max      = iaq.m_phibar_max;
    m_phigeo_bin_size = iaq.m_phigeo_bin_size;
    m_phibar_bin_size = iaq.m_phibar_bin_size;
    m_e1min           = iaq.m_e1min;
    m_e1max           = iaq.m_e1max;
    m_e2min           = iaq.m_e2min;
    m_e2max           = iaq.m_e2max;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMIaq::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Add IAQ for a given photon energy
 *
 * @param[in] energy Input photon energy.
 * @param[in] weight Spectral weight.
 *
 * The code implemented is based on the COMPASS RESPSIT2 function IAQWEI.F
 * (release 1.0, 24-FEB-93).
 ***************************************************************************/
void GCOMIaq::iaqwei(const GEnergy& energy, const double& weight)
{
    // Get response information
    GCaldb caldb("cgro","comptel");
    m_response_d1 = GCOMD1Response(caldb, "DEFAULT");
    m_response_d2 = GCOMD2Response(caldb, "DEFAULT");
    m_ict         = GCOMInstChars(caldb, "DEFAULT");

    // Get IAQ dimensions
    int n_phigeo = m_iaq.naxes(0);
    int n_phibar = m_iaq.naxes(1);

    // Set fine phibar step size to improve the accuracy
    double finstp = 0.25;
    int    nfin   = int(m_phibar_bin_size / finstp + 0.5);
    finstp        = m_phibar_bin_size / double(nfin);

    // Loop over phigeo
    for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {

        // Get geometrical scatter angle (deg)
        double phigeo = (double(i_phigeo) + 0.5) * m_phigeo_bin_size;

        // Compute the true D1 and D2 energy deposits based on the
        // geometrical scatter angle
        double etrue2 = com_energy2(energy.MeV(), phigeo);
        double etrue1 = energy.MeV() - etrue2;

        // Debug
        #if defined(G_DEBUG_IAQWEI)
        std::cout << "phigeo=" << phigeo;
        std::cout << " etrue1=" << etrue1;
        std::cout << " etrue2=" << etrue2 << std::endl;
        double sum = 0.0;
        #endif

        // Loop over phibar
        for (int i_phibar = 0; i_phibar < n_phibar; ++i_phibar) {

            // Initialise start of phibar for this layer
            double phibar = double(i_phibar) * m_phibar_bin_size + 0.5 * finstp;

            // Initialise response for this phibar layer
            double response = 0.0;

            // Loop fine phibar layers
            for (int i = 0; i < nfin; ++i, phibar += finstp) {

                // Compute response
                response += respsc(etrue1, etrue2, phibar);
//std::cout << " phibar=" << phibar << " response=" << response << std::endl;

            } // endfor: looper over fine phibar layers

            // Multiply response by size of fine phibar layers
            response *= finstp;

            // Store result
            m_iaq(i_phigeo, i_phibar) = response;

            // Debug
            #if defined(G_DEBUG_IAQWEI)
            sum += response;
            #endif

        } // endfor: looped over phibar

        // Debug
        #if defined(G_DEBUG_IAQWEI)
        std::cout << " Sum=" << sum << std::endl;
        #endif

    } // endfor: looped over phigeo

    // Return
    return;
}


/***********************************************************************//**
 * @brief Integrates response over energy interval
 *
 * @param[in] etrue1 True D1 energy deposit (MeV).
 * @param[in] etrue2 True D2 energy deposit (MeV).
 * @param[in] phibar Compton scatter angle (deg).
 * @param[in] e1min Minimum measured D1 energy deposit for integration (MeV).
 * @param[in] e1max Maximum measured D1 energy deposit for integration (MeV).
 *
 * Compute the integral
 *
 * \f[
 *    R(\bar{\varphi}|\hat{E_1},\hat{E_2}) =
 *    \int_{E_1^{\rm min}}^{E_1^{\rm max}}
 *    R(E_1,E_2|\hat{E_1},\hat{E_2}) dE_1
 * \f]
 *
 * where
 * \f$R(E_1,E_2|\hat{E_1},\hat{E_2})\f$ is computed using
 * GCOMIaq::response_kernel::eval,
 *
 * \f$E_1\f$ and \f$E_2\f$ are the measured D1 and D2 energy deposits,
 * respectively,
 * \f$\hat{E_1}\f$ and \f$\hat{E_2}\f$ are the true D1 and D2 energy
 * deposits, respectively, and \f$\bar{\varphi}\f$ the Compton scatter
 * angle.
 *
 * The code implementation is based on the COMPASS RESPSIT2 functions
 * RESPSC.F (release 2.0, 14-Mar-91) and RESINT.F (release 3.0, 21-Oct-91).
 * Since RESPSC.F only defines the threshold while RESINT.F does the real
 * job, both functions were merged.
 ***************************************************************************/
double GCOMIaq::respsc(const double& etrue1,
                       const double& etrue2,
                       const double& phibar)
{
    // Initialise response
    double response = 0.0;

    // Set integration limits
    double e1min = m_response_d1.emin(etrue1) - 0.5 * m_response_d1.ewidth(etrue1);
    double e1max = m_response_d1.emax(etrue1);

    // Set integration interval. Do not integrate more than 3 standard
    // deviations away from photo peak.
    double position = m_response_d1.position(etrue1);
    double sigma    = m_response_d1.sigma(etrue1);
    double e_low    = position - 3.0 * sigma;
    double e_high   = position + 3.0 * sigma;
    double e_min    = (e1min < e_low)  ? e_low  : e1min;
    double e_max    = (e1max > e_high) ? e_high : e1max;

    // Continue only if the integration interval is positive
    if (e_min < e_max) {

        // Setup integration kernel
        response_kernel integrand(m_response_d1,
                                  m_response_d2,
                                  etrue1,
                                  etrue2,
                                  phibar,
                                  m_ebounds.emin().MeV(),
                                  m_ebounds.emax().MeV(),
                                  m_e1min,
                                  m_e1max,
                                  m_e2min,
                                  m_e2max);

        // Setup integral
        GIntegral integral(&integrand);

        // No warnings
        #if defined(G_RESPSC_NO_WARNINGS)
        integral.silent(true);
        #endif

        // Perform integration
        response = integral.romberg(e_min, e_max);

    } // endif: integration interval was positive

    // Debug
    #if defined(G_DEBUG_RESPSC)
    std::cout << " phibar=" << phibar;
    std::cout << " e1min=" << e1min;
    std::cout << " e1max=" << e1max;
    std::cout << " e_min=" << e_min;
    std::cout << " e_max=" << e_max;
    std::cout << " response=" << response << std::endl;
    #endif

    // Return response
    return response;
}


/***********************************************************************//**
 * @brief Computes product of D1 and D2 responses at energy E1 and phibar
 *
 * @param[in] energy1 D1 energy deposit (MeV).
 *
 * The method computes
 *
 * \f[
 *    R(E_1,E_2|\hat{E_1},\hat{E_2}) =
      R_1(E_1|\hat{E_1}) \times R_2(E_2|\hat{E_2}) \times J
 * \f]
 *
 * where
 * \f$R_1(E_1|\hat{E_1})\f$ is the D1 module response function,
 * \f$R_2(E_2|\hat{E_2})\f$ is the D2 module response function,
 * \f$E_1\f$ and \f$E_2\f$ are the measured D1 and D2 energy deposits,
 * respectively,
 * \f$\hat{E_1}\f$ and \f$\hat{E_2}\f$ are the true D1 and D2 energy
 * deposits, respectively, and \f$J\f$ is the Jacobian.
 *
 * The measured D2 energy deposit is computed using
 *
 * \f[
 *    E_2 = \frac{1}{2} \times \left(
 *          \sqrt{\frac{4 m_e c^2 E_1}{1 - \cos \bar{\varphi}} + E_1^2} - E_1
 *          \right)
 * \f]
 *
 * where
 * \f$m_e\,c^2\f$ the electron rest mass, and
 * \f$\bar{\varphi}\f$ the Compton scatter angle.
 *
 * The Jacobian is given by
 *
 * \f[
 *    J = \frac{m_e c^2 \sqrt{E_1} \sin \bar{\varphi}}
 *             {(1 - \cos \bar{\varphi})^2
 *              \sqrt{\frac{4 m_e c^2}{1 - \cos \bar{\varphi}} + E_1}}
 * \f]
 *
 * The code implementation is based on the COMPASS RESPSIT2 function FUNC.F
 * (release 4.0, 21-Oct-91).
 ***************************************************************************/
double GCOMIaq::response_kernel::eval(const double& energy1)
{
    // Initialise result
    double value = 0.0;

    // Set E1 and 1-cos(phibar). If the values are too small, limit to small
    // positive numbers.
    #if defined(G_RESPONSE_KERNEL_LIMITS)
    double e1 = (energy1 < 1.0e-20) ? 1.0e-20 : energy1;
    double cp = 1.0 - m_cos_phibar;
    if (std::abs(cp) < 1.0e-20) {
        cp = 1.0e-20;
    }
    #else
    double e1 = energy1;
    double cp = 1.0 - m_cos_phibar;
    #endif

    // Compute E2 and Jacobian
    double term1 = 4.0 * gammalib::mec2 / cp;
    double e2    = 0.5 * (std::sqrt(term1 * e1 + e1*e1) - e1);
    double jc    = gammalib::mec2 * std::sqrt(e1) * m_sin_phibar /
                   (cp * cp * std::sqrt(term1 + e1));

    // Debug
    #if defined(G_DEBUG_RESPONSE_KERNEL)
    std::cout << " e1=" << e1;
    std::cout << " e2=" << e2;
    std::cout << " jc=" << jc << std::endl;
    #endif

    // Allow to break now at any time
    do {

        // Break if the D1 energy deposit is outside the energy range
        if ((e1 < m_e1min) || (e1 > m_e1max)) {
            break;
        }

        // Break if the D2 energy deposit is outside the energy range
        if ((e2 < m_e2min) || (e2 > m_e2max)) {
            break;
        }

        // Break if the total energy deposit is outside the energy range
        double etot = e1 + e2;
        if ((etot < m_etmin) || (etot > m_etmax)) {
            break;
        }

        // Break if the Jacobian precision is too poor
        #if defined(G_RESPONSE_KERNEL_LIMITS)
        if (jc < 1.0e-20) {
            break;
        }
        #endif

        // Get D1 response. Break if it is too small.
        double d1 = m_response_d1(m_etrue1, e1);
        #if defined(G_RESPONSE_KERNEL_LIMITS)
        if (d1 < 1.0e-20) {
            break;
        }
        #endif

        // Get D2 response. Break if it is too small.
        double d2 = m_response_d2(m_etrue2, e2);
        #if defined(G_RESPONSE_KERNEL_LIMITS)
        if (d2 < 1.0e-20) {
            break;
        }
        #endif

        // Compute value
        value = jc * d1 * d2;

    } while(false);

    // Return value
    return value;
}
