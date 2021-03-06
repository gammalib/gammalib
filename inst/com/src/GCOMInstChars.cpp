/***************************************************************************
 *       GCOMInstChars.cpp - COMPTEL Instrument Characteristics class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMInstChars.cpp
 * @brief COMPTEL Instrument Characteristics class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <algorithm>
#include "GMath.hpp"
#include "GEnergy.hpp"
#include "GCOMInstChars.hpp"
#include "GCOMSupport.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define USE_ICT_MEAN_FREE_PATH

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_READ_SELFVETO

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates empty COMPTEL instrument characteristics.
 ***************************************************************************/
GCOMInstChars::GCOMInstChars(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ict COMPTEL instrument characteristics.
 **************************************************************************/
GCOMInstChars::GCOMInstChars(const GCOMInstChars& ict)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(ict);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * @param[in] caldb Calibration database.
 * @param[in] ictname ICT response name.
 *
 * Create COMPTEL instrument characteristics by loading an SDA file from a
 * calibration database.
 ***************************************************************************/
GCOMInstChars::GCOMInstChars(const GCaldb& caldb, const std::string& ictname)
{
    // Initialise members
    init_members();

    // Set calibration database
    this->caldb(caldb);

    // Load instrument characteristics
    this->load(ictname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of COMPTEL instrument characteristics.
 ***************************************************************************/
GCOMInstChars::~GCOMInstChars(void)
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
 * @param[in] ict COMPTEL instrument characteristics.
 * @return COMPTEL instrument characteristics.
 ***************************************************************************/
GCOMInstChars& GCOMInstChars::operator=(const GCOMInstChars& ict)
{
    // Execute only if object is not identical
    if (this != &ict) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(ict);

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
 * Clears COMPTEL instrument characteristics object by resetting all members
 * to an initial state. Any information that was present in the object before
 * will be lost.
 ***************************************************************************/
void GCOMInstChars::clear(void)
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
 * @return Pointer to deep copy of COMPTEL instrument characteristics.
 ***************************************************************************/
GCOMInstChars* GCOMInstChars::clone(void) const
{
    return new GCOMInstChars(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL instrument characteristics.
 *
 * @param[in] ictname COMPTEL instrument characteristics.
 *
 * Loads the COMPTEL instrument characteristics with specified name
 * @p ictname. The method first searchs for an appropriate response in the
 * calibration database. If no appropriate response is found, the method
 * takes the database root path and response name to build the full path to
 * the response file, and tries to load the response from these paths.
 ***************************************************************************/
void GCOMInstChars::load(const std::string& ictname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // First attempt reading the response using the GCaldb interface
    GFilename filename = m_caldb.filename("","","ICT","","",ictname);

    // If filename is empty then build filename from CALDB root path and
    // response name
    if (filename.is_empty()) {
        filename = gammalib::filepath(m_caldb.rootdir(), ictname);
        if (!filename.exists()) {
            GFilename testname = filename + ".fits";
            if (testname.exists()) {
                filename = testname;
            }
        }
    }

    // Open FITS file
    GFits fits(filename);

    // Get ICT tables
    const GFitsTable& d1inter  = *fits.table("D1INTER");
    const GFitsTable& d1dens   = *fits.table("D1DENS");
    const GFitsTable& d1pos    = *fits.table("D1POS");
    const GFitsTable& d1rad    = *fits.table("D1RAD");
    const GFitsTable& d1thick  = *fits.table("D1THICK");
    const GFitsTable& d2inter  = *fits.table("D2INTER");
    const GFitsTable& d2dens   = *fits.table("D2DENS");
    const GFitsTable& d2pos    = *fits.table("D2POS");
    const GFitsTable& d2rad    = *fits.table("D2RAD");
    const GFitsTable& d2thick  = *fits.table("D2THICK");
    const GFitsTable& thbar    = *fits.table("THBAR");
    const GFitsTable& delz     = *fits.table("DELZ");
    const GFitsTable& alu      = *fits.table("ALU");
    const GFitsTable& aldens   = *fits.table("ALDENS");
    const GFitsTable& althick  = *fits.table("ALTHICK");
    const GFitsTable& aboved1  = *fits.table("ABOVED1");
    const GFitsTable& abdens   = *fits.table("ABDENS");
    const GFitsTable& abthick  = *fits.table("ABTHICK");
    const GFitsTable& vetodome = *fits.table("VETODOME");
    const GFitsTable& vetodens = *fits.table("VETODENS");
    const GFitsTable& v1thick  = *fits.table("V1THICK");
    const GFitsTable& vthick   = *fits.table("VTHICK");
    const GFitsTable& selfveto = *fits.table("SELFVETO");
    const GFitsTable& d1multi  = *fits.table("D1MULTI");
    const GFitsTable& d2multi  = *fits.table("D2MULTI");

    // Read information from ICT tables
    read_coeffs(d1inter, m_d1inter_energies, m_d1inter_coeffs);
    read_coeffs(d2inter, m_d2inter_energies, m_d2inter_coeffs);
    read_coeffs(alu, m_alu_energies, m_alu_coeffs);
    read_coeffs(aboved1, m_aboved1_energies, m_aboved1_coeffs);
    read_coeffs(vetodome, m_veto_energies, m_veto_coeffs);
    read_coeffs(d1multi, m_d1multi_energies, m_d1multi_coeffs);
    read_coeffs(d2multi, m_d2multi_energies, m_d2multi_coeffs);
    read_pos(d1pos, m_d1pos_x, m_d1pos_y);
    read_pos(d2pos, m_d2pos_x, m_d2pos_y);
    read_selfveto(selfveto);

    // Read constants
    m_d1dens   = d1dens["DENSITY"]->real(0);
    m_d1rad    = d1rad["RADIUS"]->real(0);
    m_d1thick  = d1thick["THICKNESS"]->real(0);
    m_d2dens   = d2dens["DENSITY"]->real(0);
    m_d2rad    = d2rad["RADIUS"]->real(0);
    m_d2thick  = d2thick["THICKNESS"]->real(0);
    m_thbar    = thbar["ANGLE"]->real(0);
    m_delz     = delz["DISTANCE"]->real(0);
    m_aldens   = aldens["DENSITY"]->real(0);
    m_althick  = althick["THICKNESS"]->real(0);
    m_abdens   = abdens["DENSITY"]->real(0);
    m_abthick  = abthick["THICKNESS"]->real(0);
    m_vetodens = vetodens["DENSITY"]->real(0);
    m_v1thick  = v1thick["THICKNESS"]->real(0);
    m_vthick   = vthick["THICKNESS"]->real(0);

    // Close ICT FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return transmission above D1
 *
 * @param[in] energy Input photon energy (MeV).
 * @return Transmission above D1.
 *
 * Computes the transmission of material above D1 as function of energy
 * using
 *
 * \f[
 *    T(E) = \exp \left(-\mu(E) l \right)
 * \f]
 *
 * where
 * \f$\mu(E)\f$ is the energy-dependent interaction coefficient for material
 * above D1 in units of \f$cm^{-1}\f$ that is interpolated using a log-log
 * interpolation of the ICT table values, and
 * \f$l\f$ is the thickness of the material above D1 in \f$cm\f$.
 ***************************************************************************/
double GCOMInstChars::trans_D1(const double& energy) const
{
    // Initialise transmission
    double transmission = 1.0;

    // Continue only if there are energies
    if (m_aboved1_energies.size() > 0) {

        // Get log of energy
        double logE = std::log(energy);

        // Get attenuation coefficient
        double logc  = m_aboved1_energies.interpolate(logE, m_aboved1_coeffs);
        double coeff = std::exp(logc) * m_abthick;

        // Compute transmission
        transmission = std::exp(-coeff);

    }

    // Return transmission
    return transmission;
}


/***********************************************************************//**
 * @brief Return V1 veto dome transmission
 *
 * @param[in] energy Input photon energy (MeV).
 * @return V1 veto dome transmission.
 *
 * Computes the V1 veto dome transmission as function of energy
 * using
 *
 * \f[
 *    T(E) = \exp \left(-\mu(E) l \right)
 * \f]
 *
 * where
 * \f$\mu(E)\f$ is the energy-dependent interaction coefficient of the V1
 * veto dome in units of \f$cm^{-1}\f$ that is interpolated using a log-log
 * interpolation of the ICT table values, and
 * \f$l\f$ is the thickness of the V1 veto dome in \f$cm\f$.
 ***************************************************************************/
double GCOMInstChars::trans_V1(const double& energy) const
{
    // Initialise transmission
    double transmission = 1.0;

    // Continue only if there are energies
    if (m_veto_energies.size() > 0) {

        // Get log of energy
        double logE = std::log(energy);

        // Get attenuation coefficient
        double logc  = m_veto_energies.interpolate(logE, m_veto_coeffs);
        double coeff = std::exp(logc) * m_v1thick;

        // Compute transmission
        transmission = std::exp(-coeff);

    }

    // Return transmission
    return transmission;
}


/***********************************************************************//**
 * @brief Return D1 interaction probability
 *
 * @param[in] energy Input photon energy (MeV).
 * @return D1 interaction probability.
 *
 * Computes the D1 interaction probability as function of energy using
 *
 * \f[
 *    P(E) = 1 - \exp \left(-\mu_m(E) \, l \right)
 * \f]
 *
 * where
 * \f$-\mu(E)\f$ is the energy-dependent D1 attenuation coefficient in units
 * of \f$1/cm\f$ that is interpolated using a log-log interpolation of the
 * ICT table values, and \f$l\f$ is the thickness of the D1 module in
 * \f$cm\f$.
 ***************************************************************************/
double GCOMInstChars::prob_D1inter(const double& energy) const
{
    // Initialise probability
    double prob = 0.0;

    // Continue only if there are energies
    if (m_d1inter_energies.size() > 0) {

        // Get log of energy
        double logE = std::log(energy);

        // Get interaction coefficient
        double logc  = m_d1inter_energies.interpolate(logE, m_d1inter_coeffs);
        double coeff = std::exp(logc) * m_d1thick;

        // Compute interaction probability
        prob = 1.0 - std::exp(-coeff);

    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return probability that no multihit occured
 *
 * @param[in] energy Input photon energy (MeV).
 * @return Probability that no multihit occured.
 *
 * Returns the probability that there is no multihit. The probability is
 * directly interpolated using a log-log interpolation from the D1 values
 * that are given in the ICT table. The D2 values are not used.
 ***************************************************************************/
double GCOMInstChars::prob_no_multihit(const double& energy) const
{
    // Initialise probability
    double prob = 1.0;

    // Continue only if there are energies
    if (m_d1multi_energies.size() > 0) {

        // Get log of energy
        double logE = std::log(energy);

        // Compute probability that there is no multihit
        prob = std::exp(m_d1multi_energies.interpolate(logE, m_d1multi_coeffs));

        // Make sure that probability does not exceed 1
        if (prob > 1.0) {
            prob = 1.0;
        }

    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return probability that photon was not self vetoed
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] zenith Zenith angle (deg).
 * @return Probability that photon was not self vetoed.
 *
 * Returns the probability that the photon was not self vetoed. The
 * probability is directly bi-linearly interpolated from the values that
 * are given in the ICT table.
 ***************************************************************************/
double GCOMInstChars::prob_no_selfveto(const double& energy, const double& zenith) const
{
    // Initialise probability
    double prob = 1.0;

    // Get number of energies and zenith angles in array
    int n_energies = m_selfveto_energies.size();
    int n_zeniths  = m_selfveto_zeniths.size();

    // Continue only if there are energies and zenith angles
    if ((n_energies > 0) && (n_zeniths > 0)) {

        // Set energy interpolation
        m_selfveto_energies.set_value(energy);

        // Set zenith angle interpolation
        m_selfveto_zeniths.set_value(zenith);

        // Set array indices for bi-linear interpolation
        int inx_zenith_left  = m_selfveto_zeniths.inx_left()  * n_energies;
        int inx_zenith_right = m_selfveto_zeniths.inx_right() * n_energies;
        int inx1 = m_selfveto_energies.inx_left()  + inx_zenith_left;
        int inx2 = m_selfveto_energies.inx_left()  + inx_zenith_right;
        int inx3 = m_selfveto_energies.inx_right() + inx_zenith_left;
        int inx4 = m_selfveto_energies.inx_right() + inx_zenith_right;

        // Set weighting factors for bi-linear interpolation
        double wgt1 = m_selfveto_energies.wgt_left()  * m_selfveto_zeniths.wgt_left();
        double wgt2 = m_selfveto_energies.wgt_left()  * m_selfveto_zeniths.wgt_right();
        double wgt3 = m_selfveto_energies.wgt_right() * m_selfveto_zeniths.wgt_left();
        double wgt4 = m_selfveto_energies.wgt_right() * m_selfveto_zeniths.wgt_right();

        // Perform bi-linear interpolation
        prob = wgt1 * m_selfveto_coeffs[inx1] +
               wgt2 * m_selfveto_coeffs[inx2] +
               wgt3 * m_selfveto_coeffs[inx3] +
               wgt4 * m_selfveto_coeffs[inx4];

        // Make sure that probability does is within [0,1]
        if (prob < 0.0) {
            prob = 0.0;
        }
        else if (prob > 1.0) {
            prob = 1.0;
        }

        // And now convert into a no self veto probability
        prob = 1.0 - prob;

    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return transmission of material between D1 and D2
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return Transmission of material between D1 and D2.
 *
 * Computes the transmission of material between D1 and D2 as function of
 * energy using
 *
 * \f[
 *    T(E) = \exp \left(\frac{-\mu(E_2) \, l \right)
 * \f]
 *
 * with
 *
 * \f[
 *    E_2 = \frac{E_{\gamma}}
 *               {(1 - \cos \varphi_{\rm geo}) \frac{E_{\gamma}}{m_e c^2} + 1}
 * \f]
 *
 * where
 * \f$E_{\gamma}\f$ is the input photon @p energy in MeV,
 * \f$\varphi_{\rm geo}\f$ is the geometrical scatter angle,
 * \f$\mu(E_2)\f$ is the energy-dependent interaction coefficient of the
 * material between D1 and D2 in units of \f$cm^{-1}\f$ that is interpolated
 * using a log-log interpolation of the ICT table values, and
 * \f$l\f$ is the thickness of the material between D1 and D2 in \f$cm\f$.
 ***************************************************************************/
double GCOMInstChars::trans_D2(const double& energy, const double& phigeo) const
{
    // Initialise transmission
    double transmission = 1.0;

    // Continue only if there are energies
    if (m_alu_energies.size() > 0) {

        // Compute log of D2 energy deposit
        double logE2 = std::log(gammalib::com_energy2(energy, phigeo));

        // Get attenuation coefficient
        double logc  = m_alu_energies.interpolate(logE2, m_alu_coeffs);
        double coeff = std::exp(logc) * m_althick;

        // Compute transmission
        transmission = std::exp(-coeff);

    }

    // Return transmission
    return transmission;
}


/***********************************************************************//**
 * @brief Return V2+V3 veto dome transmission
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return V2+V3 veto dome transmission.
 *
 * Computes the V2+V3 veto dome transmission as function of energy
 * using
 *
 * \f[
 *    T(E) = \exp \left(-\mu(E) \, l \right)
 * \f]
 *
 * \f[
 *    E_2 = \frac{E_{\gamma}}
 *               {(1 - \cos \varphi_{\rm geo}) \frac{E_{\gamma}}{m_e c^2} + 1}
 * \f]
 *
 * where
 * \f$E_{\gamma}\f$ is the input photon @p energy in MeV,
 * \f$\varphi_{\rm geo}\f$ is the geometrical scatter angle,
 * \f$\mu(E_2)\f$ is the energy-dependent interaction coefficient of the V2
 * and V3 veto domes in units of \f$cm^{-1}\f$ that is interpolated using a
 * log-log interpolation of the ICT table values, and \f$l\f$ is the thickness
 * of the V2+V3 veto domes in \f$cm\f$.
 ***************************************************************************/
double GCOMInstChars::trans_V23(const double& energy, const double& phigeo) const
{
    // Initialise transmission
    double transmission = 1.0;

    // Continue only if there are energies
    if (m_veto_energies.size() > 0) {

        // Compute log of D2 energy deposit
        double logE2 = std::log(gammalib::com_energy2(energy, phigeo));

        // Get attenuation coefficient
        double logc  = m_veto_energies.interpolate(logE2, m_veto_coeffs);
        double coeff = std::exp(logc) * m_vthick;

        // Compute transmission
        transmission = std::exp(-coeff);

    }

    // Return transmission
    return transmission;
}


/***********************************************************************//**
 * @brief Return D2 interaction probability
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return D2 interaction probability.
 *
 * Computes the D2 interaction probability as function of energy using
 *
 * \f[
 *    P(E) = 1 - \exp \left(-\mu_m(E_2) \, l \right)
 * \f]
 *
 * with
 *
 * \f[
 *    E_2 = \frac{E_{\gamma}}
 *               {(1 - \cos \varphi_{\rm geo}) \frac{E_{\gamma}}{m_e c^2} + 1}
 * \f]
 *
 * where
 * \f$E_{\gamma}\f$ is the input photon @p energy in MeV,
 * \f$\varphi_{\rm geo}\f$ is the geometrical scatter angle,
 * \f$\mu(E_2)\f$ is the energy-dependent D2 attenuation coefficient in units
 * of \f$cm^{-1}\f$ that is interpolated using a log-log interpolation of the
 * ICT table values, and \f$l\f$ is the thickness of the D2 module in
 * \f$cm\f$.
 ***************************************************************************/
double GCOMInstChars::prob_D2inter(const double& energy, const double& phigeo) const
{
    // Initialise probability
    double prob = 0.0;

    // Continue only if there are energies
    if (m_d2inter_energies.size() > 0) {

        // Compute log of D2 energy deposit
        double logE2 = std::log(gammalib::com_energy2(energy, phigeo));

        // Get interaction coefficient
        double logc  = m_d2inter_energies.interpolate(logE2, m_d2inter_coeffs);
        double coeff = std::exp(logc) * m_d2thick;

        // Compute interaction probability
        prob = 1.0 - std::exp(-coeff);

    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return multi-scatter correction
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return Correction factor due to multi-scatter.
 *
 * Returns the fraction of photons that have undergone a single scatter
 * and which leave the D1 module unattenuated (no second interaction when
 * traversing the remaining path inside the same module).
 *
 * RES already calculates the fraction of photons which undergo a single
 * scatter for VERTICALLY indicent photons, based on the Klein-Nishina
 * cross-section and the composition of NE213A. Therefore the above mentioned
 * transmission can be applied as a multiplicative correction per phigeo.
 * However, in order to calculate that correction, an integral over the
 * module must be performed, which properly takes into account the radiative
 * transfer inside the module at all R and z. We again use the assumption
 * of vertical photon incidence, which simplifies the calculation. The
 * integral is done over the location of the first interaction.
 *
 * Note that the transmission calculated is conservative : in reality it will
 * be a bit higher because of a fraction of the photons which undergo a
 * second scatter have a final scatter angle after escaping the module within
 * the phigeo-range of the PSF. These events will, however, mainly be located
 * at large phibar (large D1E deposit).
 *
 * The code implementation is based on the COMPASS RESPSIT2 function
 * MULTIS.F (release ?, 27-NOV-92).
 ***************************************************************************/
double GCOMInstChars::multi_scatter(const double& energy, const double& phigeo) const
{
    // Set integration step size
    const int nz   = 10; //!< Number of steps in z-direction
    const int nphi = 20; //!< Number of azimuthal steps

    // Compute stuff
    double alpha = energy / gammalib::mec2;
    double sth   = std::sin(phigeo * gammalib::deg2rad);
    double cth   = std::cos(phigeo * gammalib::deg2rad);

    // Get the attenuation coefficient for NE213A
    double mu = 1.0 / ne213a_mfpath(energy);

    // Compute the energy after scattering and the corresponding attenuation
    double e_scattered  = energy / (1.0 + alpha * (1.0 - cth));
    double mu_scattered = 1.0 / ne213a_mfpath(e_scattered);

    // Compute the full vertical optical depth of the module
    double tau0 = mu * m_d1thick;

    // Compute the fraction of photons that are interacting
    double total_interactions = 1.0 - std::exp(-tau0);

    // Perform integration per rho-value (symmetry!) over
    // (1) z (the geometrical depth in the module)
    // (2) phi (azimuth of scatter)
    // at geometrical scatter angle phigeo.
    double deltaz = m_d1thick / double(nz);
    double deltau = deltaz * mu;
    double dltphi = gammalib::pi / double(nphi);

    // Initialise integration loop. The binning in radius is based on
    // constant surface per cylinder
    double rho     = 1.0;
    double delrho  = 2.0 * rho;
    double rholow  = rho - 0.5 * delrho; // Should be 0
    double rhoupp  = rho + 0.5 * delrho; // Should be 2
    double surface = 2.0 * delrho * rho; // Should be 4

    // Compute the integration normalisation value
    double kappa = deltau / double(nphi);

    // Integration loop over the module radius
    double contribution_rho = 0.0;
    double total_surface    = surface;
    int    klast            = 0;   // Will be >0 if the loop should be exited
    for (int krho = 0; krho < 100; ++krho) {

        // If we have reached the end then exit the loop now
        if (klast > 0) {
            break;
        }

        // Determine radius limits. If D1 module radius is reached or
        // exceeded, limit the upper rho value to the D1 module radius,
        // re-compute the surface, and signal to exit the loop in the
        // next round.
        if (krho > 0) {
            rholow = rhoupp;
            rhoupp = std::sqrt(surface + rholow*rholow);
            rho    = 0.5 * (rhoupp + rholow);
            if (rhoupp > m_d1rad) {
                rhoupp  = m_d1rad;
                rho     = 0.5 * (rhoupp + rholow);
                surface = (rhoupp*rhoupp - rholow*rholow);
                klast   = krho;
            }
            total_surface += surface;
        }

        // Initialise azimuthal results
        std::vector<double> r(nphi, 0.0);

        // Compute remaining radius as function of azimuth angle if the
        // first interaction was at radius rho.
        for (int lphi = 0; lphi < nphi; ++lphi) {
            double phi  = dltphi * (lphi + 0.5);
            double term = rho * std::sin(phi);
            r[lphi]     = -rho * std::cos(phi) +
                                 std::sqrt(m_d1rad*m_d1rad - term*term);
        }

        // Perform integration over depth
        double contribution_z = 0.0;
        for (int kz = 0; kz < nz; ++kz) {

            // Calculate depth
            double z   = (kz + 0.5) * deltaz;
            double tau = z * mu;

            // Compute the remaining length in the module after the first
            // interaction. Test for foward/backward scattering although
            // this should not be relevant for any reasonable phigeo values.
            double length0 = 1000.0;
            if (sth < 0.99) {                   // True for all reasonable phigeo's
                if (cth < 0.0) {                // False for all reasonable phigeo's
                    length0 = z/std::abs(cth);
                }
                else {
                    length0 = (m_d1thick-z)/cth; // All reasonable phigeo's
                }
            }

            // Perform integration over azimuth
            double contribution_phi = 0.0;
            for (int kphi = 0; kphi < nphi; ++kphi) {

                // Compute the actual remaining length for a given azimuth
                // angle phi. Limit the length to the remaining radius.
                double length = length0;
                if (length * sth > r[kphi]) {
                    length = r[kphi] / sth;
                }

                // Now compute the probability that the photon is not
                // absorbed during the remaining path through the detector.
                contribution_phi += std::exp(-length * mu_scattered);

            } // endfor: looped over azimuth

            // Average contribution must later be divided by nphi for
            // phi-pixels and multiplied by deltau; this is done by
            // multiplying with kappa. Do this at end, to save computation
            contribution_z += contribution_phi * std::exp(-tau);

        } // endfor: looped over depth

        // Add the average contribution of radius rho
        contribution_rho += contribution_z * kappa * surface;

    } // endfor: looped over radius

    // Compute average
    double transmission = contribution_rho / (total_surface * total_interactions);

    // Return transmission
    return transmission;
}


/***********************************************************************//**
 * @brief Return PSD correction
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return Correction factor due to PSD selection.
 *
 * Returns the D1 energy dependent PSD correction as described in
 * COM-RP-ROL-DRG-016 and COM-RP-ROL-DRG-35. It applies to a standard PSD
 * selection of 0-110.
 *
 * The acceptance probability fit formula
 *
 * \f[
 *    P_{\rm acc} = 1 - \frac{1}{a_1 \times E_1^{a_2} + 1}
 * \f]
 *
 * where \f$a_1=1727.9\f$, \f$a_2=2.53\f$ and \f$E_1\f$ is the D1 energy
 * deposit in MeV. Coefficients are taken from Boron calibration data
 * (ZA=1.14) to remain consistent with Rob van Dijk's SIMPSF corrections.
 *
 * The code implementation is based on the COMPASS RESPSIT2 function
 * PSDACP.F (release 1.0, 11-DEC-92).
 ***************************************************************************/
double GCOMInstChars::psd_correction(const double& energy, const double& phigeo) const
{
    // Set constants
    const double a1 = 1727.9;
    const double a2 = 2.530;

    // Compute D1 energy deposit
    double e1 = gammalib::com_energy1(energy, phigeo);

    // Original COMPASS code
    double psdacp = 1.0 - (1.0 / (a1 * std::pow(e1, a2) + 1.0));

    // Return fraction
    return psdacp;
}


/***********************************************************************//**
 * @brief Print COMPTEL instrument characteristics
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL instrument characteristics.
 ***************************************************************************/
std::string GCOMInstChars::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMInstChars ===");

        // Append D1INTER information
        result.append("\n"+gammalib::parformat("D1 interaction coeffs."));
        if (m_d1inter_energies.size() > 0) {
            int last = m_d1inter_energies.size()-1;
            result.append(gammalib::str(std::exp(m_d1inter_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_d1inter_energies[last]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_d1inter_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_d1inter_coeffs))));
            result.append("] cm^2/g");
        }
        else {
            result.append("not defined");
        }

        // Append D2INTER information
        result.append("\n"+gammalib::parformat("D2 interaction coeffs."));
        if (m_d2inter_energies.size() > 0) {
            int last = m_d2inter_energies.size()-1;
            result.append(gammalib::str(std::exp(m_d2inter_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_d2inter_energies[last]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_d2inter_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_d2inter_coeffs))));
            result.append("] cm^2/g");
        }
        else {
            result.append("not defined");
        }

        // Append ALU information
        result.append("\n"+gammalib::parformat("Al interaction coeffs."));
        if (m_alu_energies.size() > 0) {
            result.append(gammalib::str(std::exp(m_alu_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_alu_energies[m_alu_energies.size()-1]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_alu_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_alu_coeffs))));
            result.append("] 1/cm");
        }
        else {
            result.append("not defined");
        }

        // Append ABOVED1 information
        result.append("\n"+gammalib::parformat("Above D1 atten. coeffs."));
        if (m_aboved1_energies.size() > 0) {
            result.append(gammalib::str(std::exp(m_aboved1_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_aboved1_energies[m_aboved1_energies.size()-1]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_aboved1_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_aboved1_coeffs))));
            result.append("] 1/cm");
        }
        else {
            result.append("not defined");
        }

        // Append VETODOME information
        result.append("\n"+gammalib::parformat("Vetodome atten. coeffs."));
        if (m_veto_energies.size() > 0) {
            result.append(gammalib::str(std::exp(m_veto_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_veto_energies[m_veto_energies.size()-1]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_veto_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_veto_coeffs))));
            result.append("] cm^2/g");
        }
        else {
            result.append("not defined");
        }

        // Append SELFVETO information
        result.append("\n"+gammalib::parformat("Selfveto probabilities"));
        if (m_selfveto_energies.size() > 0) {
            result.append(gammalib::str(m_selfveto_energies[0])+ " - ");
            result.append(gammalib::str(m_selfveto_energies[m_selfveto_energies.size()-1])+ " MeV x ");
            if (m_selfveto_zeniths.size() > 0) {
                result.append(gammalib::str(m_selfveto_zeniths[0])+ " - ");
                result.append(gammalib::str(m_selfveto_zeniths[m_selfveto_zeniths.size()-1])+ " deg");
            }
            result.append(" [");
            result.append(gammalib::str(min_coeff(m_selfveto_coeffs)));
            result.append(" - ");
            result.append(gammalib::str(max_coeff(m_selfveto_coeffs)));
            result.append("]");
        }
        else {
            result.append("not defined");
        }

        // Append D1MULTI information
        result.append("\n"+gammalib::parformat("D1 multihit probabilities"));
        if (m_d1multi_energies.size() > 0) {
            result.append(gammalib::str(std::exp(m_d1multi_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_d1multi_energies[m_d1multi_energies.size()-1]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_d1multi_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_d1multi_coeffs))));
            result.append("]");
        }
        else {
            result.append("not defined");
        }

        // Append D2MULTI information
        result.append("\n"+gammalib::parformat("D2 multihit probabilities"));
        if (m_d2multi_energies.size() > 0) {
            result.append(gammalib::str(std::exp(m_d2multi_energies[0]))+ " - ");
            result.append(gammalib::str(std::exp(m_d2multi_energies[m_d2multi_energies.size()-1]))+ " MeV");
            result.append(" [");
            result.append(gammalib::str(std::exp(min_coeff(m_d2multi_coeffs))));
            result.append(" - ");
            result.append(gammalib::str(std::exp(max_coeff(m_d2multi_coeffs))));
            result.append("]");
        }
        else {
            result.append("not defined");
        }

        // Append D1POS information
        result.append("\n"+gammalib::parformat("D1 positions"));
        result.append(gammalib::str(m_d1pos_x.size())+"x & ");
        result.append(gammalib::str(m_d1pos_y.size())+"y");

        // Append D2POS information
        result.append("\n"+gammalib::parformat("D2 positions"));
        result.append(gammalib::str(m_d2pos_x.size())+"x & ");
        result.append(gammalib::str(m_d2pos_y.size())+"y");

        // Append D1DENS information
        result.append("\n"+gammalib::parformat("D1 density"));
        result.append(gammalib::str(m_d1dens)+ " g/cm^-3");

        // Append D1RAD information
        result.append("\n"+gammalib::parformat("D1 radius"));
        result.append(gammalib::str(m_d1rad)+ " cm");

        // Append D1THICK information
        result.append("\n"+gammalib::parformat("D1 thickness"));
        result.append(gammalib::str(m_d1thick)+ " cm");

        // Append D2DENS information
        result.append("\n"+gammalib::parformat("D2 density"));
        result.append(gammalib::str(m_d2dens)+ " g/cm^-3");

        // Append D2RAD information
        result.append("\n"+gammalib::parformat("D2 radius"));
        result.append(gammalib::str(m_d2rad)+ " cm");

        // Append D2THICK information
        result.append("\n"+gammalib::parformat("D2 thickness"));
        result.append(gammalib::str(m_d2thick)+ " cm");

        // Append THBAR information
        result.append("\n"+gammalib::parformat("Average D2 incident angle"));
        result.append(gammalib::str(m_thbar)+ " deg");

        // Append DELZ information
        result.append("\n"+gammalib::parformat("Distance between D1 and D2"));
        result.append(gammalib::str(m_delz)+ " cm");

        // Append ALDENS information
        result.append("\n"+gammalib::parformat("Dens. of Al plate above D2"));
        result.append(gammalib::str(m_aldens)+ " g/cm^-3");

        // Append ALTHICK information
        result.append("\n"+gammalib::parformat("Thk. of Al plate above D2"));
        result.append(gammalib::str(m_althick)+ " cm");

        // Append ABDENS information
        result.append("\n"+gammalib::parformat("Density above D1"));
        result.append(gammalib::str(m_abdens)+ " g/cm^-3");

        // Append ABTHICK information
        result.append("\n"+gammalib::parformat("Thickness above D1"));
        result.append(gammalib::str(m_abthick)+ " cm");

        // Append VETODENS information
        result.append("\n"+gammalib::parformat("Density of veto domes"));
        result.append(gammalib::str(m_vetodens)+ " g/cm^-3");

        // Append V1THICK information
        result.append("\n"+gammalib::parformat("Thickness of V1"));
        result.append(gammalib::str(m_v1thick)+ " cm");

        // Append VTHICK information
        result.append("\n"+gammalib::parformat("Thickness of V2+V3"));
        result.append(gammalib::str(m_vthick)+ " cm");

        // Append calibration database
        result.append("\n"+m_caldb.print(chatter));

        // Append information

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
void GCOMInstChars::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_d1inter_energies.clear();
    m_d1inter_coeffs.clear();
    m_d2inter_energies.clear();
    m_d2inter_coeffs.clear();
    m_alu_energies.clear();
    m_alu_coeffs.clear();
    m_aboved1_energies.clear();
    m_aboved1_coeffs.clear();
    m_veto_energies.clear();
    m_veto_coeffs.clear();
    m_selfveto_energies.clear();
    m_selfveto_zeniths.clear();
    m_selfveto_coeffs.clear();
    m_d1multi_energies.clear();
    m_d1multi_coeffs.clear();
    m_d2multi_energies.clear();
    m_d2multi_coeffs.clear();
    m_d1pos_x.clear();
    m_d1pos_y.clear();
    m_d2pos_x.clear();
    m_d2pos_y.clear();
    m_d1dens   = 0.0;
    m_d1rad    = 0.0;
    m_d1thick  = 0.0;
    m_d2dens   = 0.0;
    m_d2rad    = 0.0;
    m_d2thick  = 0.0;
    m_thbar    = 0.0;
    m_delz     = 0.0;
    m_aldens   = 0.0;
    m_althick  = 0.0;
    m_abdens   = 0.0;
    m_abthick  = 0.0;
    m_vetodens = 0.0;
    m_v1thick  = 0.0;
    m_vthick   = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ict COMPTEL instrument characteristics.
 ***************************************************************************/
void GCOMInstChars::copy_members(const GCOMInstChars& ict)
{
    // Copy attributes
    m_caldb             = ict.m_caldb;
    m_d1inter_energies  = ict.m_d1inter_energies;
    m_d1inter_coeffs    = ict.m_d1inter_coeffs;
    m_d2inter_energies  = ict.m_d2inter_energies;
    m_d2inter_coeffs    = ict.m_d2inter_coeffs;
    m_alu_energies      = ict.m_alu_energies;
    m_alu_coeffs        = ict.m_alu_coeffs;
    m_aboved1_energies  = ict.m_aboved1_energies;
    m_aboved1_coeffs    = ict.m_aboved1_coeffs;
    m_veto_energies     = ict.m_veto_energies;
    m_veto_coeffs       = ict.m_veto_coeffs;
    m_selfveto_energies = ict.m_selfveto_energies;
    m_selfveto_zeniths  = ict.m_selfveto_zeniths;
    m_selfveto_coeffs   = ict.m_selfveto_coeffs;
    m_d1multi_energies  = ict.m_d1multi_energies;
    m_d1multi_coeffs    = ict.m_d1multi_coeffs;
    m_d2multi_energies  = ict.m_d2multi_energies;
    m_d2multi_coeffs    = ict.m_d2multi_coeffs;
    m_d1pos_x           = ict.m_d1pos_x;
    m_d1pos_y           = ict.m_d1pos_y;
    m_d2pos_x           = ict.m_d2pos_x;
    m_d2pos_y           = ict.m_d2pos_y;
    m_d1dens            = ict.m_d1dens;
    m_d1rad             = ict.m_d1rad;
    m_d1thick           = ict.m_d1thick;
    m_d2dens            = ict.m_d2dens;
    m_d2rad             = ict.m_d2rad;
    m_d2thick           = ict.m_d2thick;
    m_thbar             = ict.m_thbar;
    m_delz              = ict.m_delz;
    m_aldens            = ict.m_aldens;
    m_althick           = ict.m_althick;
    m_abdens            = ict.m_abdens;
    m_abthick           = ict.m_abthick;
    m_vetodens          = ict.m_vetodens;
    m_v1thick           = ict.m_v1thick;
    m_vthick            = ict.m_vthick;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMInstChars::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dependent coefficients.
 *
 * @param[in] table FITS table.
 * @param[out] energies Energy node array.
 * @param[out] coeffs Coefficients.
 *
 * Read energy dependent coefficients from FITS table and store their natural
 * logarithm in the @p energies and @p coeffs vectors.
 ***************************************************************************/
void GCOMInstChars::read_coeffs(const GFitsTable&    table,
                                GNodeArray&          energies,
                                std::vector<double>& coeffs)
{
    // Initialise energies and coefficients
    energies.clear();
    coeffs.clear();

    // Extract number of entries in table
    int num = table.nrows();

    // If there are entries then read them
    if (num > 0) {

        // By default use "COEFFICIENT" column for coefficients, but if
        // table contains a "PROBABILITY" column then use that column
        std::string coeff_name = "COEFFICIENT";
        if (table.contains("PROBABILITY")) {
            coeff_name = "PROBABILITY";
        }

        // Get column pointers
        const GFitsTableCol* ptr_energy = table["ENERGY"];
        const GFitsTableCol* ptr_coeff  = table[coeff_name];

        // Copy data from table into vectors. Skip any non-positive values
        // since we store the logarithms for a log-log interpolation.
        for (int i = 0; i < num; ++i) {
            double energy = ptr_energy->real(i);
            double coeff  = ptr_coeff->real(i);
            if ((energy > 0.0) && (coeff > 0.0)) {
                energies.append(std::log(energy));
                coeffs.push_back(std::log(coeff));
            }
        }

    } // endif: there were entries

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read module positions.
 *
 * @param[in] table FITS table.
 * @param[out] x X-positions of module.
 * @param[out] y Y-positions of module.
 *
 * Read module positions from FITS table.
 ***************************************************************************/
void GCOMInstChars::read_pos(const GFitsTable& table, std::vector<double>& x,
                             std::vector<double>& y)
{
    // Initialise positions
    x.clear();
    y.clear();

    // Extract number of entries in table
    int num = table.nrows();

    // If there are entries then read them
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_detnum = table["DETNUM"];
        const GFitsTableCol* ptr_x      = table["X"];
        const GFitsTableCol* ptr_y      = table["Y"];

        // Copy data from table into vectors
        for (int i = 0; i < num; ++i) {
            int detnum = ptr_detnum->integer(i) - 1;
            x.push_back(ptr_x->real(detnum));
            y.push_back(ptr_y->real(detnum));
        }

    } // endif: there were entries

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read selfveto coefficients.
 *
 * @param[in] table FITS table.
 *
 * Read selfveto coefficients from FITS table. The selfveto coefficients
 * depend on energy and zenith angle, and the input data are not necessarily
 * regularly sampled. Therefore the method will form a regular grid from
 * the provided values and compute the interpolated values from that grid.
 ***************************************************************************/
void GCOMInstChars::read_selfveto(const GFitsTable& table)
{
    // Initialise selfveto vectors
    m_selfveto_energies.clear();
    m_selfveto_zeniths.clear();
    m_selfveto_coeffs.clear();

    // Extract number of entries in table
    int num = table.nrows();

    // If there are entries then read them
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_energies = table["ENERGY"];
        const GFitsTableCol* ptr_zeniths  = table["ZENITH"];
        const GFitsTableCol* ptr_coeffs   = table["PROBABILITY"];

        // Initialise energies and zenith vectors for sorting before creating
        // node arrays
        std::vector<double> energies;
        std::vector<double> zeniths;

        // Determine all energies and zenith angles in table
        for (int i = 0; i < num; ++i) {

            // Add energy if it is not yet in array
            double energy = ptr_energies->real(i);
            bool   add    = true;
            for (int k = 0; k < energies.size(); ++k) {
                if (energies[k] == energy) {
                    add = false;
                    break;
                }
            }
            if (add) {
                energies.push_back(energy);
            }

            // Add zenith angle if it is not yet in array
            double zenith = ptr_zeniths->real(i);
            add           = true;
            for (int k = 0; k < m_selfveto_zeniths.size(); ++k) {
                if (m_selfveto_zeniths[k] == zenith) {
                    add = false;
                    break;
                }
            }
            if (add) {
                zeniths.push_back(zenith);
            }

        } // endfor: collected all energies and zenith angles

        // Sort energies and zenith angles
        std::sort(energies.begin(), energies.end());
        std::sort(zeniths.begin(), zeniths.end());

        // Set node arrays
        m_selfveto_energies.nodes(energies);
        m_selfveto_zeniths.nodes(zeniths);

        // Get array size
        int neng    = m_selfveto_energies.size();
        int nzenith = m_selfveto_zeniths.size();

        // Fill all coefficients with -1, meaning that there is no value
        m_selfveto_coeffs.assign(neng*nzenith, -1.0);

        // Fill all coefficients
        for (int i = 0; i < num; ++i) {

            // Find nearest energy
            int    ieng   = 0;
            double energy = ptr_energies->real(i);
            double delta  = std::abs(energy - m_selfveto_energies[0]);
            for (int k = 1; k < neng; ++k) {
                double d = std::abs(energy - m_selfveto_energies[k]);
                if (d < delta) {
                    delta = d;
                    ieng  = k;
                }
            }

            // Find nearest zenith angle
            int    izenith = 0;
            double zenith  = ptr_zeniths->real(i);
            delta          = std::abs(zenith - m_selfveto_zeniths[0]);
            for (int k = 1; k < nzenith; ++k) {
                double d = std::abs(zenith - m_selfveto_zeniths[k]);
                if (d < delta) {
                    delta   = d;
                    izenith = k;
                }
            }

            // Compute index
            int index = ieng + neng*izenith;

            // Set coefficient
            m_selfveto_coeffs[index] = ptr_coeffs->real(i);

        } // endfor: filled all coefficients

        // Now loop over all energies and fill the missing zenith angles
        for (int ieng = 0; ieng < neng; ++ieng) {

            // Debug: print header
            #if defined(G_DEBUG_READ_SELFVETO)
            std::cout << "Energy: " << m_selfveto_energies[ieng] << std::endl;
            std::cout << "  Zenith angles:";
            #endif

            // Set up arrays for interpolation
            GNodeArray          nodes;
            std::vector<double> values;
            for (int izenith = 0; izenith < nzenith; ++izenith) {
                int index = ieng + neng*izenith;
                if (m_selfveto_coeffs[index] > -0.5) {
                    nodes.append(m_selfveto_zeniths[izenith]);
                    values.push_back(m_selfveto_coeffs[index]);
                    #if defined(G_DEBUG_READ_SELFVETO)
                    std::cout << " " << m_selfveto_zeniths[izenith];
                    std::cout << "(" << m_selfveto_coeffs[index] << ")";
                    #endif
                }
            }

            // Debug: print header
            #if defined(G_DEBUG_READ_SELFVETO)
            std::cout << std::endl;
            std::cout << "  Interpolated zenith angles:";
            #endif

            // Interpolate coefficients for zenith angles
            for (int izenith = 0; izenith < nzenith; ++izenith) {
                int index = ieng + neng*izenith;
                if (m_selfveto_coeffs[index] < -0.5) {
                    double coeffs = nodes.interpolate(m_selfveto_zeniths[izenith], values);
                    if (coeffs < 0.0) {
                        coeffs = 0.0;
                    }
                    else if (coeffs > 1.0) {
                        coeffs = 1.0;
                    }
                    m_selfveto_coeffs[index] = coeffs;
                    #if defined(G_DEBUG_READ_SELFVETO)
                    std::cout << " " << m_selfveto_zeniths[izenith];
                    std::cout << "(" << m_selfveto_coeffs[index] << ")";
                    #endif
                }
            }

            // Debug: print line feed
            #if defined(G_DEBUG_READ_SELFVETO)
            std::cout << std::endl;
            #endif

        } // endfor: looped over all energies

    } // endif: there were entries

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return NE213A mean free path
 *
 * @param[in] energy Energy (MeV)
 * @return Mean free path (cm)
 *
 * The mean free path as a function of energy, linearly interpolated between
 * data-points. Table of mean free paths for NE213A taken from L. Kuiper
 * calculated in JOB ROL-SPDSD1-43. This includes (1) Compton scattering and
 * (2) Pair production.
 ***************************************************************************/
double GCOMInstChars::ne213a_mfpath(const double& energy) const
{
    #if defined(USE_ICT_MEAN_FREE_PATH)
    // Get log of energy
    double logE = std::log(energy);

    // Get log of interaction coefficient
    double logc = m_d1inter_energies.interpolate(logE, m_d1inter_coeffs);

    // Convert into mean free path
    double result = 1.0 / std::exp(logc);

    #else
    // NE213A mean free path values for energies from 1,2,...,30 MeV
    const double mfpath[] = {16.19, 23.21, 29.01, 34.04, 38.46,
                             42.35, 45.64, 48.81, 51.53, 54.17,
                             56.15, 58.09, 59.96, 61.74, 63.41,
                             64.62, 65.80, 66.94, 68.04, 69.09,
                             69.86, 70.61, 71.33, 72.03, 72.70,
                             73.33, 73.93, 74.50, 75.03, 75.51};

    // Normalisation constant for low-energy extrapolation
    const double norm = mfpath[1] * kn_cross_section(2.0 / gammalib::mec2);

    // Initialise mfpath
    double result = 0.0;

    // Do a Klein-Nishina extrapolation for low energies
    if (energy < 2.0) {
        result   = norm / kn_cross_section(energy / gammalib::mec2);
    }

    // ... otherwise do a linear interpolation
    else {
        int index = int(energy) - 1;
        if (index >= 29) {
            index = 28;
        }
        double slope = mfpath[index+1] - mfpath[index];
        result       = mfpath[index] + slope * (energy - double(index + 1));
    }
    #endif

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return integrated Klein-Nishina cross section
 *
 * @param[in] k \f$E/ m_e c^2\f$
 * @return Integrated Klein-Nishina cross section
 *
 * Computes
 *
 * \f[
 *    \sigma_{\rm KN}(k) = \frac{1+k}{k^2}
 *                         \left[
 *                               \frac{2(1+k)}{1+2k} - \frac{\ln(1+2k)}{k}
 *                         \right] + \frac{\ln(1+2k)}{2k} -
 *                                   \frac{1+3k}{(1+2k)^2}
 * \f]
 *
 * where \f$k = E/ m_e c^2\f$.
 ***************************************************************************/
double GCOMInstChars::kn_cross_section(const double& k) const
{
    // Compute integrated Klein-Nishina cross section
    double k1      = k + 1.0;
    double k2p1    = 1.0 + 2.0*k;
    double logk2p1 = std::log(k2p1);
    double kn      = k1/(k*k) * (2.0*k1/k2p1 - logk2p1/k) +
                     logk2p1/(2.0*k) - (1.0+3.0*k)/(k2p1*k2p1);

    // Return
    return kn;
}


/***********************************************************************//**
 * @brief Return minimum coefficient.
 *
 * @param[in] coeffs Coefficients.
 * @return Minimum coefficient.
 *
 * Returns minimum coefficient.
 ***************************************************************************/
double GCOMInstChars::min_coeff(const std::vector<double>& coeffs) const
{
    // Initialise minimum
    double min = 0.0;

    // Continue only if there are coefficients
    if (coeffs.size() > 0) {

        // Set first value as minimum
        min = coeffs[0];

        // Search for minimum
        for (int i = 1; i < coeffs.size(); ++i) {
            if (coeffs[i] < min) {
                min = coeffs[i];
            }
        }

    }

    // Return minimum
    return min;
}


/***********************************************************************//**
 * @brief Return maximum coefficient.
 *
 * @param[in] coeffs Coefficients.
 * @return Maximum coefficient.
 *
 * Returns maximum coefficient.
 ***************************************************************************/
double GCOMInstChars::max_coeff(const std::vector<double>& coeffs) const
{
    // Initialise maximum
    double max = 0.0;

    // Continue only if there are coefficients
    if (coeffs.size() > 0) {

        // Set first value as maximum
        max = coeffs[0];

        // Search for maximum
        for (int i = 1; i < coeffs.size(); ++i) {
            if (coeffs[i] > max) {
                max = coeffs[i];
            }
        }

    }

    // Return maximum
    return max;
}
