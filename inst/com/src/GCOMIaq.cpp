/***************************************************************************
 *         GCOMIaq.cpp - COMPTEL instrument response representation        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2019 by Juergen Knoedlseder                         *
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
#include "GEbounds.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlawPhotonFlux.hpp"
#include "GModelSpectralPlawEnergyFlux.hpp"
#include "GCOMIaq.hpp"
#include "GCOMSupport.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_RESPONSE_KERNEL_LIMITS
//#define G_APPLY_PSD_CORRECTION_IN_RESPONSE_KERNEL
#define G_COMPUTE_IAQ_BIN_NO_WARNINGS
#define G_LOCATION_SMEARING_NO_WARNINGS

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_COMPTON_KINEMATICS
//#define G_DEBUG_COMPUTE_IAQ_BIN
//#define G_DEBUG_KLEIN_NISHINA
//#define G_DEBUG_KLEIN_NISHINA_INTEGRAL
//#define G_DEBUG_RESPONSE_KERNEL
//#define G_DEBUG_WEIGHT_IAQ
//#define G_DEBUG_SET_CONTINUUM
//#define G_DEBUG_LOCATION_SMEARING
//#define G_DEBUG_LOCATION_SPREAD

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
 *
 * @todo Test input argument validity.
 **************************************************************************/
GCOMIaq::GCOMIaq(const double&   phigeo_max, const double& phigeo_bin_size,
                 const double&   phibar_max, const double& phibar_bin_size)
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

    // Allocate FITS image
    m_iaq = GFitsImageFloat(nphigeo, nphibar);

    // Set FITS image keywords
    m_iaq.card("CTYPE1", "Phigeo", "Geometrical scatter angle");
    m_iaq.card("CRVAL1", 0.5*m_phigeo_bin_size, "[deg] Geometrical scatter angle for reference bin");
    m_iaq.card("CDELT1", m_phigeo_bin_size, "[deg] Geometrical scatter angle bin size");
    m_iaq.card("CRPIX1", 1.0, "Reference bin of geometrical scatter angle");
    m_iaq.card("CTYPE2", "Phibar", "Compton scatter angle");
    m_iaq.card("CRVAL2", 0.5*m_phibar_bin_size, "[deg] Compton scatter angle for reference bin");
    m_iaq.card("CDELT2", m_phibar_bin_size, "[deg] Compton scatter angle bin size");
    m_iaq.card("CRPIX2", 1.0, "Reference bin of Compton scatter angle");
    m_iaq.card("BUNIT", "Probability per bin", "Unit of bins");
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
 * @brief Set mono-energetic IAQ
 *
 * @param[in] energy Input photon energy.
 * @param[in] ebounds Boundaries of observed energy.
 *
 * The code implemented is based on the COMPASS RESPSIT2 function SPCIAQ.F
 * (release 1.0, 18-DEC-92).
 *
 * @todo Implement geometrical smearing.
 ***************************************************************************/
void GCOMIaq::set(const GEnergy& energy, const GEbounds& ebounds)
{
    // Store the energy boundaries
    m_ebounds = ebounds;

    // Initialise COMPTEL response information and instrument characteristics
    init_response();

    // Remove any extra header cards
    remove_cards();

    // Generate IAQ using Compton kinematics
    compton_kinematics(energy.MeV());

    // Multiply IAQ with Klein-Nishina factors
    klein_nishina(energy.MeV());

    // Weight IAQ
    weight_iaq(energy.MeV());

    // Apply location smearing
    location_smearing(m_zenith);

    // Set source photon energy
    m_iaq.card("ENERGY", energy.MeV(), "[MeV] Source photon energy");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set continuum IAQ
 *
 * @param[in] spectrum Input spectrum.
 * @param[in] ebounds Boundaries of observed energy.
 *
 * Computes the continuum IAQ based on an input @p spectrum. The method
 * computes the line IAQ for a m_num_energies logarithmically spaced energies
 * within an energy range that conceivably covers the observed energy band.
 *
 * The code implementation is loosly based on the COMPASS RESPSIT2 functions
 * SPCIAQ.F, GENWEI.F and DERFAN.F (release 1.0, 18-DEC-92). The code was
 * much simplified.
 *
 * @todo Implement geometrical smearing.
 ***************************************************************************/
void GCOMIaq::set(const GModelSpectral& spectrum, const GEbounds& ebounds)
{
    // Store the energy boundaries
    m_ebounds = ebounds;

    // Initialise COMPTEL response information and instrument characteristics
    init_response();

    // Remove any extra header cards
    remove_cards();

    // Get minimum and maximum input energy for continuum IAQ
    double energy_min = 0.7 * m_ebounds.emin().MeV();
    double energy_max = (m_response_d1.emax() < m_response_d2.emax()) ?
                         m_response_d1.emax() : m_response_d2.emax();

    // Debug
    #if defined(G_DEBUG_SET_CONTINUUM)
    std::cout << "Minimum energy=" << energy_min << std::endl;
    std::cout << "Maximum energy=" << energy_max << std::endl;
    #endif

    // Setup array of energy boundaries
    GEbounds ebds(m_num_energies, GEnergy(energy_min, "MeV"),
                                  GEnergy(energy_max, "MeV"));

    // Compute flux over measured total energy interval
    double flux_total = spectrum.flux(m_ebounds.emin(), m_ebounds.emax());

    // Store empty copy of IAQ
    GFitsImageFloat iaq = m_iaq;
    for (int k = 0; k < iaq.npix(); ++k) {
        iaq(k) = 0.0;
    }

    // Loop over energy bins
    for (int i = 0; i < ebds.size(); ++i) {

        // Get log mean energy for bin
        GEnergy energy = ebds.elogmean(i);

        // Get weight for this bin
        double weight = spectrum.flux(ebds.emin(i), ebds.emax(i)) /
                        flux_total;

        // Initialise sum
        double sum = 0.0;

        // Continue only if we have weight
        if (weight > 0.0) {

            // Generate IAQ using Compton kinematics
            compton_kinematics(energy.MeV());

            // Multiply IAQ with Klein-Nishina factors
            klein_nishina(energy.MeV());

            // Weight IAQ
            weight_iaq(energy.MeV());

            // Copy weigthed IAQ values
            for (int k = 0; k < iaq.npix(); ++k) {
                double value = m_iaq(k) * weight;
                iaq(k)      += value;
                sum         += value;
            }

        } // endif: we has weight

        // Debug
        #if defined(G_DEBUG_SET_CONTINUUM)
        std::cout << "Energy=" << energy;
        std::cout << " weight=" << weight;
        std::cout << " sum=" << sum << std::endl;
        #endif

    } // endfor: looped over energy bins

    // Store IAQ matrix
    m_iaq = iaq;

    // Apply location smearing
    location_smearing(m_zenith);

    // Set parameters
    m_iaq.card("SPECTRUM", spectrum.classname(), "Source spectrum");
    const GModelSpectralPlaw*           plaw  =
          dynamic_cast<const GModelSpectralPlaw*>(&spectrum);
    const GModelSpectralPlawPhotonFlux* pplaw =
          dynamic_cast<const GModelSpectralPlawPhotonFlux*>(&spectrum);
    const GModelSpectralPlawEnergyFlux* eplaw =
          dynamic_cast<const GModelSpectralPlawEnergyFlux*>(&spectrum);
    if (plaw != NULL) {
        m_iaq.card("PLAWINX", plaw->index(), "Power law spectral index");
    }
    else if (pplaw != NULL) {
        m_iaq.card("PLAWINX", pplaw->index(), "Power law spectral index");
    }
    else if (eplaw != NULL) {
        m_iaq.card("PLAWINX", eplaw->index(), "Power law spectral index");
    }
    m_iaq.card("NENG", m_num_energies, "Number of incident energies");
    m_iaq.card("EIMIN", energy_min, "[MeV] Minimum incident energy");
    m_iaq.card("EIMAX", energy_max, "[MeV] Maximum incident energy");

    // Return
    return;
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

    // Set FITS image keywords
    iaq.card("EMIN", m_ebounds.emin().MeV(), "[MeV] Minimum measured photon energy");
    iaq.card("EMAX", m_ebounds.emax().MeV(), "[MeV] Maximum measured photon energy");
    iaq.card("E1MIN", m_e1min, "[MeV] Minimum D1 energy deposit");
    iaq.card("E1MAX", m_e1max, "[MeV] Maximum D1 energy deposit");
    iaq.card("E2MIN", m_e1min, "[MeV] Minimum D2 energy deposit");
    iaq.card("E2MAX", m_e1max, "[MeV] Maximum D2 energy deposit");
    iaq.card("ZENITH", m_zenith, "[deg] Zenith angle for location smearing");
    iaq.card("PSDCORR", m_psd_correct, "PSD correction for 0-110");
    iaq.card("PHIRES", m_phibar_resolution, "[deg] Phibar resolution in computation");

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

        // Append parameters
        result.append("\n"+gammalib::parformat("D1 energy range"));
        result.append(gammalib::str(m_e1min)+" MeV - ");
        result.append(gammalib::str(m_e1max)+" MeV");
        result.append("\n"+gammalib::parformat("D2 energy range"));
        result.append(gammalib::str(m_e2min)+" MeV - ");
        result.append(gammalib::str(m_e2max)+" MeV");
        result.append("\n"+gammalib::parformat("Zenith angle"));
        result.append(gammalib::str(m_zenith)+" deg");
        result.append("\n"+gammalib::parformat("Phibar resolution"));
        result.append(gammalib::str(m_phibar_resolution)+" deg");
        result.append("\n"+gammalib::parformat("PSD correction"));
        if (m_psd_correct) {
            result.append("yes");
        }
        else {
            result.append("no");
        }
        result.append("\n"+gammalib::parformat("Number of input energies"));
        result.append(gammalib::str(m_num_energies));

        // Append IAQ integral
        double sum = 0.0;
        for (int i = 0; i < m_iaq.npix(); ++i) {
            sum += m_iaq(i);
        }
        result.append("\n"+gammalib::parformat("Response integral"));
        result.append(gammalib::str(sum));

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

    // Initialise parameters
    m_phibar_resolution =  0.25;  //!< Default: 0.25 deg
    m_e1min             =  0.070; //!< Default: 70 keV
    m_e1max             = 20.0;   //!< Default: 20 MeV
    m_e2min             =  0.650; //!< Default: 650 keV
    m_e2max             = 30.0;   //!< Default: 30 MeV
    m_num_energies      = 50;     //!< Default: 50 input energies
    m_psd_correct       = true;   //!< Default: use PSD correction
    m_zenith            = 25.0;   //!< Default: 25 deg

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

    // Copy parameters
    m_phibar_resolution = iaq.m_phibar_resolution;
    m_e1min             = iaq.m_e1min;
    m_e1max             = iaq.m_e1max;
    m_e2min             = iaq.m_e2min;
    m_e2max             = iaq.m_e2max;
    m_num_energies      = iaq.m_num_energies;
    m_psd_correct       = iaq.m_psd_correct;
    m_zenith            = iaq.m_zenith;

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
 * @brief Initialise COMPTEL response function and instrument characteristics
 ***************************************************************************/
void GCOMIaq::init_response(void)
{
    // Initialise COMPTEL response function  and instrument characteristics
    GCaldb caldb("cgro","comptel");
    m_response_d1 = GCOMD1Response(caldb, "DEFAULT");
    m_response_d2 = GCOMD2Response(caldb, "DEFAULT");
    m_ict         = GCOMInstChars(caldb, "DEFAULT");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove any extra header cards
 *
 * Remove all header cards that may have been attached by the set() methods.
 ***************************************************************************/
void GCOMIaq::remove_cards(void)
{
    // Remove extra header cards
    if (m_iaq.has_card("ENERGY")) {
        const_cast<GFitsHeader&>(m_iaq.header()).remove("ENERGY");
    }
    if (m_iaq.has_card("SPECTRUM")) {
        const_cast<GFitsHeader&>(m_iaq.header()).remove("SPECTRUM");
    }
    if (m_iaq.has_card("PLAWINX")) {
        const_cast<GFitsHeader&>(m_iaq.header()).remove("PLAWINX");
    }
    if (m_iaq.has_card("NENG")) {
        const_cast<GFitsHeader&>(m_iaq.header()).remove("NENG");
    }
    if (m_iaq.has_card("EIMIN")) {
        const_cast<GFitsHeader&>(m_iaq.header()).remove("EIMIN");
    }
    if (m_iaq.has_card("EIMAX")) {
        const_cast<GFitsHeader&>(m_iaq.header()).remove("EIMAX");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate IAQ matrix based on Compton kinematics
 *
 * @param[in] energy Input photon energy (MeV).
 *
 * Generates an IAQ matrix based on Compton kinematics only. See the
 * compute_iaq_bin() method for the formula used to compute each bin of the
 * IAQ matrix.
 *
 * The Phibar dimension of the IAQ is sampled at a resolution that is set
 * by the m_phibar_resolution member.
 *
 * The code implemented is based on the COMPASS RESPSIT2 function IAQWEI.F
 * (release 1.0, 24-FEB-93).
 ***************************************************************************/
void GCOMIaq::compton_kinematics(const double& energy)
{
    // Get IAQ dimensions
    int n_phigeo = m_iaq.naxes(0);
    int n_phibar = m_iaq.naxes(1);

    // Set the phibar step size for the internal computation. The phibar
    // step size is set by the m_phibar_resolution member. This ensures that
    // the computation is done with a sufficient resolution in phibar.
    int n_fine = int(m_phibar_bin_size / m_phibar_resolution + 0.5);
    if (n_fine < 1) {
        n_fine = 1;
    }
    double dphibar     = m_phibar_bin_size / double(n_fine);
    double dphibar_rad = dphibar * gammalib::deg2rad;

    // Loop over phigeo
    for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {

        // Get geometrical scatter angle (deg)
        double phigeo = (double(i_phigeo) + 0.5) * m_phigeo_bin_size;

        // Compute the true D1 and D2 energy deposits based on the
        // geometrical scatter angle
        double etrue2 = gammalib::com_energy2(energy, phigeo);
        double etrue1 = energy - etrue2;

        // Debug
        #if defined(G_DEBUG_COMPTON_KINEMATICS)
        std::cout << "phigeo=" << phigeo;
        std::cout << " etrue1=" << etrue1;
        std::cout << " etrue2=" << etrue2;
        #endif

        // Initialise sum
        double sum = 0.0;

        // Loop over phibar
        for (int i_phibar = 0; i_phibar < n_phibar; ++i_phibar) {

            // Initialise start of phibar for this layer
            double phibar = double(i_phibar) * m_phibar_bin_size + 0.5 * dphibar;

            // Initialise response for this phibar layer
            double response = 0.0;

            // Loop fine phibar layers
            for (int i = 0; i < n_fine; ++i, phibar += dphibar) {

                // Compute IAQ bin
                response += compute_iaq_bin(etrue1, etrue2, phibar);

            } // endfor: looper over fine phibar layers

            // Multiply response by size of fine phibar layers in radians
            response *= dphibar_rad;

            // Store result
            m_iaq(i_phigeo, i_phibar) = response;

            // Add response to sum
            sum += response;

        } // endfor: looped over phibar

        // Debug
        #if defined(G_DEBUG_COMPTON_KINEMATICS)
        std::cout << " Prob=" << sum << std::endl;
        #endif

    } // endfor: looped over phigeo

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes the IAQ for one bin
 *
 * @param[in] etrue1 True D1 energy deposit (MeV).
 * @param[in] etrue2 True D2 energy deposit (MeV).
 * @param[in] phibar Compton scatter angle (deg).
 * @return Response for one IAQ bin
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
 *
 * @todo Study impact of integration precision.
 ***************************************************************************/
double GCOMIaq::compute_iaq_bin(const double& etrue1,
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

        // Set precision
        integral.eps(1.0e-4);

        // No warnings
        #if defined(G_COMPUTE_IAQ_BIN_NO_WARNINGS)
        integral.silent(true);
        #endif

        // Perform integration
        response = integral.romberg(e_min, e_max);

    } // endif: integration interval was positive

    // Debug
    #if defined(G_DEBUG_COMPUTE_IAQ_BIN)
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
 * @brief Multiply IAQ matrix by the Klein-Nishina formula
 *
 * @param[in] energy Input photon energy (MeV).
 *
 * The code implemented is based on the COMPASS RESPSIT2 function IAQWEI.F
 * (release 1.0, 24-FEB-93).
 ***************************************************************************/
void GCOMIaq::klein_nishina(const double& energy)
{
    // Get IAQ dimensions
    int n_phigeo = m_iaq.naxes(0);
    int n_phibar = m_iaq.naxes(1);

    // Debug
    #if defined(G_DEBUG_KLEIN_NISHINA)
    double sum = 0.0;
    #endif

    // Loop over phigeo
    for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {

        // Get geometrical scatter angle (deg)
        double phigeo = (double(i_phigeo) + 0.5) * m_phigeo_bin_size;

        // Get Klein-Nishina value for Phigeo
        double prob_kn = klein_nishina_bin(energy, phigeo);

        // Debug
        #if defined(G_DEBUG_KLEIN_NISHINA)
        std::cout << "phigeo=" << phigeo << " prob_kn=" << prob_kn << std::endl;
        sum += prob_kn;
        #endif

        // Loop over phibar
        for (int i_phibar = 0; i_phibar < n_phibar; ++i_phibar) {

            // Store result
            m_iaq(i_phigeo, i_phibar) *= prob_kn;

        } // endfor: looped over phibar

    } // endfor: looped over phigeo

    // Debug
    #if defined(G_DEBUG_KLEIN_NISHINA)
    std::cout << "Sum of probabilities = " << sum << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes Klein-Nishina probability for one bin
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometric scatter angle (deg).
 * @return Klein-Nishina probability.
 *
 * Compute the probability that a photon of @p energy is Compton scattered
 * in a given \f$\varphi_{\rm geo}\f$ bin using
 *
 * \f[
 *    P = \frac{\int_{\varphi_{\rm geo,min}}^{\varphi_{\rm geo,max}}
 *              \sigma_{\rm KN}(E, \varphi_{\rm geo}) d\varphi_{\rm geo}}
 *             {\int_{0}^{\pi}
 *              \sigma_{\rm KN}(E, \varphi_{\rm geo}) d\varphi_{\rm geo}}
 * \f]
 *
 * where \f$\sigma_{\rm KN}(E, \varphi_{\rm geo})\f$ is the Klein-Nishina
 * cross section.
 *
 * The code implementation is based on the COMPASS RESPSIT2 function
 * RESPSG.F (release 2.0, 26-Apr-90).
 ***************************************************************************/
double GCOMIaq::klein_nishina_bin(const double& energy, const double& phigeo)
{
    // Compute phigeo boundaries in radians. Make sure that they remain in
    // the range [0,pi].
    double phigeo_lo = phigeo - 0.5 * m_phigeo_bin_size;
    double phigeo_hi = phigeo + 0.5 * m_phigeo_bin_size;
    if (phigeo_lo < 0.0) {
        phigeo_lo = 0.0;
    }
    if (phigeo_hi < 0.0) {
        phigeo_hi = 0.0;
    }
    if (phigeo_lo > 180.0) {
        phigeo_lo = 180.0;
    }
    if (phigeo_hi > 180.0) {
        phigeo_hi = 180.0;
    }
    phigeo_lo *= gammalib::deg2rad;
    phigeo_hi *= gammalib::deg2rad;

    // Express input photon energy in units of electron rest mass
    double alpha = energy / gammalib::mec2;

    // Calculate the probability of phigeo bin
    double v_low    = 1.0 - 1.0 / (1.0 + (1.0 - std::cos(phigeo_lo)) * alpha);
    double v_high   = 1.0 - 1.0 / (1.0 + (1.0 - std::cos(phigeo_hi)) * alpha);
    double r_low    = klein_nishina_integral(v_low,  alpha);
    double r_high   = klein_nishina_integral(v_high, alpha);
    double prob_bin = (r_high - r_low);

    // Calculate the probability of [0,pi] range
    const double v_zero   = 0.0;
    const double v_pi     = 1.0 - 1.0 / (1.0 + 2.0 * alpha);
    const double r_zero   = klein_nishina_integral(v_zero, alpha);
    const double r_pi     = klein_nishina_integral(v_pi, alpha);
    const double prob_tot = (r_pi - r_zero);

    // Calculate the probability of phigeo bin
    double prob_phigeo = prob_bin / prob_tot;

    // Return probability of phigeo bin
    return prob_phigeo;
}


/***********************************************************************//**
 * @brief Computes Klein-Nishina integral
 *
 * @param[in] v Integration limit.
 * @param[in] a Normalized energy.
 * @return Integrated Klein-Nishina cross section.
 *
 * Computes the indefinite integral of the Klein-Nishina cross section
 * evaluated at the integration limit @p v.
 *
 * The code implementation is based on the COMPASS RESPSIT2 function
 * RESSLV.F (release 1.0, 23-Nov-88). The formula used in that code is
 * the following:
 *
 * W = V*( 1.D0/A  +  1.D0 )**2
 * W = W + ( 2.D0/(A**2) + 2.D0/A - 1.D0 )* DLOG(1.D0 - V)
 * W = W - (V**2)/2.D0
 * W = W + 1.D0/(A**2) * 1.D0/(1.D0 - V)
 *
 * It has been check that the original formula gives the same results as the
 * version that was implemented.
 ***************************************************************************/
double GCOMIaq::klein_nishina_integral(const double& v, const double& a)
{
    // Compute integral (hope this is okay :o)
    double w = v * (1.0/a + 1.0) * (1.0/a + 1.0) +
               (2.0/(a*a) + 2.0/a - 1.0) * std::log(1.0 - v) -
               v*v/2.0 +
               1.0/(a*a) * 1.0/(1.0 - v);

    // Debug
    #if defined(G_DEBUG_KLEIN_NISHINA_INTEGRAL)
    double w_test;
    w_test = v*( 1.0/a  +  1.0 )*( 1.0/a  +  1.0 );
    w_test = w_test + ( 2.0/(a*a) + 2.0/a - 1.0 )* std::log(1.0 - v);
    w_test = w_test - (v*v)/2.0;
    w_test = w_test + 1.0/(a*a) * 1.0/(1.0 - v);
    if (w != w_test) {
        std::cout << "w=" << w << " w_test=" << w_test << std::endl;
    }
    #endif

    // Return integral
    return w;
}


/***********************************************************************//**
 * @brief Weight IAQ matrix
 *
 * @param[in] energy Input photon energy (MeV).
 *
 * The code implemented is based on the COMPASS RESPSIT2 function IAQWEI.F
 * (release 1.0, 24-FEB-93).
 ***************************************************************************/
void GCOMIaq::weight_iaq(const double& energy)
{
    // Compute the total D1 interaction coefficient for source on axis
    double d1prob = m_ict.prob_D1inter(energy);

    // Compute the fraction of Compton interactions within the total D1
    // interaction probability (COM-RP-ROL-DRG-41)
    double cfract = 1.0;
    if (energy >= 2.0) {
        cfract = 1.067 - 0.0295 * energy + 3.4e-4 * energy * energy;
    }
    if (cfract > 1.0) {
        cfract = 1.0;
    }
    else if (cfract < 0.0) {
        cfract = 0.0;
    }

    // Compute the transmission of material above D1
    double ad1trans = m_ict.trans_D1(energy);

    // Compute the transmission of V1 veto dome
    //double v1trans = m_ict.trans_V1(energy);

    // Compute the probability that a photon was not self vetoed assuming
    // a source on axis
    double sveto = m_ict.prob_no_selfveto(energy, 0.0);

    // Compute the probability that an event is not a multihit
    double mlhit = m_ict.prob_no_multihit(energy);

    // Compute the overall (shape independent) transmission coefficients
    //double oalltr = ad1trans * v1trans * d1prob * cfract * mlhit * sveto;
    double oalltr = ad1trans * d1prob * cfract * mlhit * sveto;

    // Debug
    #if defined(G_DEBUG_WEIGHT_IAQ)
    std::cout << "Transmission coefficients:" << std::endl;
    std::cout << "==========================" << std::endl;
    std::cout << " Above D1 transmission ..........: " << ad1trans << std::endl;
    std::cout << " V1 veto dome transmission ......: " << v1trans << std::endl;
    std::cout << " D1 interaction probability .....: " << d1prob << std::endl;
    std::cout << " Compton scatter fraction .......: " << cfract << std::endl;
    std::cout << " Compton interaction probability : " << d1prob*cfract << std::endl;
    std::cout << " Multihit transmission ..........: " << mlhit << std::endl;
    std::cout << " Self-vetoing transmission ......: " << sveto << std::endl;
    std::cout << " Overall shape-independent trans.: " << oalltr << std::endl;
    #endif

    // Get IAQ dimensions
    int n_phigeo = m_iaq.naxes(0);
    int n_phibar = m_iaq.naxes(1);

    // Loop over phigeo
    for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {

        // Get geometrical scatter angle (deg)
        double phigeo = (double(i_phigeo) + 0.5) * m_phigeo_bin_size;

        // Compute the transmission between D1 and D2
        double ad2trans = m_ict.trans_D2(energy, phigeo);

        // Compute the transmission of V2+V3 veto domes
        //double v23trans = m_ict.trans_V23(energy, phigeo);

        // Compute the D2 interaction coefficient
        double d2prob = m_ict.prob_D2inter(energy, phigeo);

        // Compute multi-scatter transmission inside D1 module
        double mscatt = m_ict.multi_scatter(energy, phigeo);

        // Optionally compute PSD correction for the default corrected PSD
        // channel selection 0-110
        #if defined(G_APPLY_PSD_CORRECTION_IN_RESPONSE_KERNEL)
        double psdtrn = 1.0;
        #else
        double psdtrn = (m_psd_correct) ? m_ict.psd_correction(energy, phigeo) : 1.0;
        #endif

        // Compute the overall phigeo dependent correction
        //double oallpg = ad2trans * v23trans * d2prob * mscatt * psdtrn;
        double oallpg = ad2trans * d2prob * mscatt * psdtrn;

        // Compute the overall correction
        double weight = oalltr * oallpg;

        // Debug
        #if defined(G_DEBUG_WEIGHT_IAQ)
        std::cout << " Phigeo .........................: " << phigeo << std::endl;
        std::cout << "  Between D1 & D2 transmission ..: " << ad2trans << std::endl;
        std::cout << "  V2+V3 veto dome transmission ..: " << v23trans << std::endl;
        std::cout << "  D2 interaction probability ....: " << d2prob << std::endl;
        std::cout << "  D1 multi-scatter transmission .: " << mscatt << std::endl;
        std::cout << "  PSD correction ................: " << psdtrn << std::endl;
        std::cout << "  Overall shape-dependent trans. : " << oallpg << std::endl;
        std::cout << "  Overall transmission ..........: " << weight << std::endl;
        #endif

        // Apply weight to all phibar layers
        for (int i_phibar = 0; i_phibar < n_phibar; ++i_phibar) {
            m_iaq(i_phigeo, i_phibar) *= weight;
        }

    } // endfor: looped over phigeo

    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform location smearing
 *
 * @param[in] zenith Zenith angle of source (deg).
 *
 * The code implementation is based on the COMPASS RESPSIT2 function LOCSPR.F
 * (release 1.0, 05-JAN-93).
 ***************************************************************************/
void GCOMIaq::location_smearing(const double& zenith)
{
    // Compute location spread for all phigeo angles
    std::vector<double> sigmas = location_spread(zenith);

    // Get IAQ dimensions
    int n_phigeo = m_iaq.naxes(0);
    int n_phibar = m_iaq.naxes(1);

    // Setup phigeo node array for interpolated
    GNodeArray phigeos;
    for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {
        phigeos.append((double(i_phigeo) + 0.5) * m_phigeo_bin_size);
    }

    // Loop over phibar
    for (int i_phibar = 0; i_phibar < n_phibar; ++i_phibar) {

        // Copy phigeo vector
        std::vector<double> values;
        double              sum_before = 0.0;
        for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {
            values.push_back(m_iaq(i_phigeo, i_phibar));
            sum_before += m_iaq(i_phigeo, i_phibar);
        }

        // Compute convolution integral for each phigeo pixel
        std::vector<double> convolved_values;
        for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {

            // Setup integration kernel
            smearing_kernel integrand(phigeos[i_phigeo],
                                      sigmas[i_phigeo],
                                      phigeos,
                                      values);

            // Setup integral
            GIntegral integral(&integrand);

            // Set precision
            integral.eps(1.0e-4);

            // No warnings
            #if defined(G_LOCATION_SMEARING_NO_WARNINGS)
            integral.silent(true);
            #endif

            // Get integration boundaries
            double phigeo_min = phigeos[i_phigeo] - 3.0 * sigmas[i_phigeo];
            double phigeo_max = phigeos[i_phigeo] + 3.0 * sigmas[i_phigeo];

            // Perform integration
            double value = integral.romberg(phigeo_min, phigeo_max);

            // Store result
            convolved_values.push_back(value);

        } // endfor: convolution integral

        // Restore phigeo vector
        double sum_after = 0.0;
        for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {
            m_iaq(i_phigeo, i_phibar) = convolved_values[i_phigeo];
            sum_after += convolved_values[i_phigeo];
        }

        // Debug
        #if defined(G_DEBUG_LOCATION_SMEARING)
        std::cout << "phibar=" << (double(i_phibar) + 0.5) * m_phibar_bin_size;
        std::cout << " before=" << sum_before;
        std::cout << " after=" << sum_after;
        if (sum_before > 0.0) {
            std::cout << " (" << sum_after/sum_before << ")";
        }
        std::cout << std::endl;
        #endif

    } // endfor: looped over phibar

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute location spread vector
 *
 * @param[in] zenith Zenith angle of source (deg).
 * @return Vector of location spreads in Gaussian sigma
 *
 * The code implementation is based on the COMPASS RESPSIT2 function LOCSPR.F
 * (release 1.0, 05-JAN-93).
 ***************************************************************************/
std::vector<double> GCOMIaq::location_spread(const double& zenith) const
{
    // Set distance between the detector midplanes of D1 and D2 (cm)
    const double zdist = 158.0;

    // Set horizontal (x,y) location spread in D1 resp. D2 (cm)
    const double sighd1 = 2.30;
    const double sighd2 = 1.96;

    // Set vertical (z) location spread in D1 resp. D2 based on homogeneous
    // event distributions in z-direction (cm)
    const double sigvd1 = 2.45;
    const double sigvd2 = 2.17;

    // Derive constants
    const double zdist2 = zdist * zdist;
    const double sight2 = sighd1*sighd1 + sighd2*sighd2;
    const double sigvt2 = sigvd1*sigvd1 + sigvd2*sigvd2;

    // Initialise result
    std::vector<double> sigmas;

    // Compute terms
    double a1     = std::sin(zenith * gammalib::deg2rad);
    double a1sq   = a1 * a1;
    double a2     = std::cos(zenith * gammalib::deg2rad);
    double a2sq   = a2 * a2;

    // Get IAQ dimensions
    int n_phigeo = m_iaq.naxes(0);
    //int n_phibar = m_iaq.naxes(1);

    // Loop over phigeo
    for (int i_phigeo = 0; i_phigeo < n_phigeo; ++i_phigeo) {

        // Get geometrical scatter angle (deg)
        double phigeo = (double(i_phigeo) + 0.5) * m_phigeo_bin_size;

        // Compute spread due to (x,y) uncertainty
        double cphig  = std::cos(phigeo * gammalib::deg2rad);
        double sphig  = std::sin(phigeo * gammalib::deg2rad);
        double cphig2 = cphig*cphig;
        double sphig2 = sphig*sphig;
        double sigxy2 = (a2sq * cphig2 *
                         (4.0 * a1sq * sphig2 + cphig2 - a1sq) +
                          0.5 * a1sq *
                         (a2sq * cphig2 + a1sq * sphig2 * (1.0 - 0.75 * cphig2))) *
                        sight2 / zdist2;

        // Compute spread due to (z) uncertainty
        double sigz2 = (a2sq * a2sq * sphig2 * cphig2 +
                        a2sq * a2sq * (0.5 - 3.0 * cphig2 * sphig2) +
                        3.0 * a1sq * a1sq * sphig2 * cphig2 / 8.0) *
                       sigvt2 / zdist2;
    
        // Calculate sigma in degrees
        double siggeo = std::sqrt(sigxy2 + sigz2) * gammalib::rad2deg;

        // Append sigma
        sigmas.push_back(siggeo);

        // Debug
        #if defined(G_DEBUG_LOCATION_SPREAD)
        std::cout << "phigeo=" << phigeo << " sigma=" << siggeo << std::endl;
        #endif

        // We should now create a warning if siggeo+zenith is larger than
        // 90.0 deg, at least the original code does that

    } // endfor: looped over phigeo

    // Return sigmas
    return sigmas;
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
 *    R_1(E_1|\hat{E_1}) \times R_2(E_2|\hat{E_2}) \times J
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

        // Optionally apply PSD correction
        #if defined(G_APPLY_PSD_CORRECTION_IN_RESPONSE_KERNEL)
        const double a1 = 1727.9;
        const double a2 = 2.530;
        d1             *= 1.0 - (1.0 / (a1 * std::pow(e1, a2) + 1.0));
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


/***********************************************************************//**
 * @brief Computes product of D1 and D2 responses at energy E1 and phibar
 *
 * @param[in] phigeo Geometrical scatter angle.
 * @return Bla ...
 *
 ***************************************************************************/
double GCOMIaq::smearing_kernel::eval(const double& phigeo)
{
    // Get interpolated IAQ value for geometrical scatter angle. Make sure
    // that the value is positive
    double value = m_phigeos.interpolate(phigeo, m_values);
    if (value < 0.0) {
        value = 0.0;
    }

    // Multiply with Gaussian
    double arg  = (m_phigeo - phigeo) * m_wgt;
    value      *= m_norm * std::exp(-0.5 * arg * arg);

    // Return value
    return value;
}
