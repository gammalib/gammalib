/***************************************************************************
 *       GCOMInstChars.cpp - COMPTEL Instrument Characteristics class      *
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
 * @brief Return D1 interaction probability
 *
 * @param[in] energy Input photon energy.
 * @return D1 interaction probability.
 *
 * Returns D1 interaction probability as function of energy for a source that
 * is on-axis.
 *
 * The computation uses a log-log interpolation of the ICT table values.
 ***************************************************************************/
double GCOMInstChars::prob_D1inter(const GEnergy& energy) const
{
    // Initialise probability
    double prob = 0.0;

    // Continue only if there are energies
    if (m_d1inter_energies.size() > 0) {

        // Get log of energy
        double logE = energy.log10MeV() * gammalib::inv_loge;

        // Get interaction coefficient
        double logc   = m_d1inter_energies.interpolate(logE, m_d1inter_coeffs);
        double coeff  = std::exp(logc) * m_d1dens * m_d1thick;

        // Compute interaction probability
        prob = 1.0 - std::exp(-coeff);

    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return D2 interaction probability
 *
 * @param[in] energy Input photon energy.
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return D2 interaction probability.
 *
 * Returns D2 interaction probability as function of energy. The computation
 * uses a log-log interpolation of the ICT table values.
 *
 * @todo Verify formula.
 ***************************************************************************/
double GCOMInstChars::prob_D2inter(const GEnergy& energy, const double& phigeo) const
{
    // Initialise probability
    double prob = 0.0;

    // Continue only if there are energies
    if (m_d2inter_energies.size() > 0) {

        // Compute log of D2 energy deposit
        double logE2 = std::log(com_energy2(energy.MeV(), phigeo));

        // Get interaction coefficient
        double logc   = m_d2inter_energies.interpolate(logE2, m_d2inter_coeffs);
        double length = m_d2thick / std::cos(phigeo * gammalib::deg2rad);
        double coeff  = std::exp(logc) * m_d2dens * length;

        // Compute interaction probability
        prob = 1.0 - std::exp(-coeff);

    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return multihit probability
 *
 * @param[in] energy Input photon energy.
 * @return Multihit probability.
 *
 * @todo Implement method. Need to find out what exactly is computed.
 ***************************************************************************/
double GCOMInstChars::prob_multihit(const GEnergy& energy) const
{
    // Initialise probability
    double prob = 1.0;

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Return attenuation above D1
 *
 * @param[in] energy Input photon energy.
 * @return Attenuation above D1.
 *
 * @todo Implement method
 ***************************************************************************/
double GCOMInstChars::atten_D1(const GEnergy& energy) const
{
    // Initialise attenuation
    double attenuation = 1.0;

    // Return attenuation
    return attenuation;
}


/***********************************************************************//**
 * @brief Return attenuation between D1 and D2
 *
 * @param[in] energy Input photon energy.
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return Attenuation between D1 and D2.
 *
 * @todo Implement method
 ***************************************************************************/
double GCOMInstChars::atten_D2(const GEnergy& energy, const double& phigeo) const
{
    // Initialise attenuation
    double attenuation = 1.0;

    // Return attenuation
    return attenuation;
}


/***********************************************************************//**
 * @brief Return attenuation due to selfveto
 *
 * @param[in] energy Input photon energy.
 * @param[in] zenith Zenith angle (deg).
 * @return Attenuation due to selfveto.
 *
 * @todo Implement method
 ***************************************************************************/
double GCOMInstChars::atten_selfveto(const GEnergy& energy, const double& zenith) const
{
    // Initialise attenuation
    double attenuation = 1.0;

    // Return attenuation
    return attenuation;
}


/***********************************************************************//**
 * @brief Return multi-scatter correction
 *
 * @param[in] energy Input photon energy.
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
 * of vertical photon incidence, which simplifies the calculation.
 *
 * Note that the transmission calculated is conservative : in reality it will
 * be a bit higher because of a fraction of the photons which undergo a
 * second scatter have a final scatter angle after escaping the module within
 * the phigeo-range of the PSF. These events will,  however, mainly be located
 * at large phibar (large D1E deposit).
 *
 * @todo Implement method from the COMPASS RESPSIT2 function MULTIS.F
 * (27-NOV-92).
 ***************************************************************************/
double GCOMInstChars::multi_scatter(const GEnergy& energy, const double& phigeo) const
{
    // Initialise fraction
    double fraction = 1.0;

    // Return fraction
    return fraction;
}


/***********************************************************************//**
 * @brief Return PSD correction
 *
 * @param[in] energy Input photon energy.
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
 *    P_{\rm acc} = 1 - \frac{1}{a_1 \times E_1^a_2 + 1.0}
 * \f]
 *
 * where \f$a_1=1727.9\f$, \f$a_2=2.53\f$ and \f$E_1\f$ is the D1 energy
 * deposit in MeV. Coefficients are taken from Boron calibration data
 * (ZA=1.14) to remain consistent with Rob van Dijk's SIMPSF corrections.
 *
 * The code was implemented from the COMPASS RESPSIT2 function PSDACP.F
 * (release 1.0, 11-DEC-92).
 ***************************************************************************/
double GCOMInstChars::psd_correction(const GEnergy& energy, const double& phigeo) const
{
    // Set constants
    const double a1 = 1727.9;
    const double a2 = 2.530;

    // Compute D1 energy deposit
    double e1 = com_energy1(energy.MeV(), phigeo);

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
        result.append("\n"+gammalib::parformat("Selfveto coeffs."));
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
            result.append("] 1/cm");
        }
        else {
            result.append("not defined");
        }

        // Append D1MULTI information
        result.append("\n"+gammalib::parformat("D1 multihit atten. coeffs."));
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
        result.append("\n"+gammalib::parformat("D2 multihit atten. coeffs."));
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

        // Get column pointers
        const GFitsTableCol* ptr_energy = table["ENERGY"];
        const GFitsTableCol* ptr_coeff  = table["COEFFICIENT"];

        // Copy data from table into vectors
        for (int i = 0; i < num; ++i) {
            energies.append(std::log(ptr_energy->real(i)));
            coeffs.push_back(std::log(ptr_coeff->real(i)));
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
        const GFitsTableCol* ptr_coeffs   = table["COEFFICIENT"];

        // Initialise energies vector (for sorting before creating a node array)
        std::vector<double> energies;

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
                m_selfveto_zeniths.push_back(zenith);
            }

        } // endfor: collected all energies and zenith angles

        // Sort energies and zenith angles
        std::sort(energies.begin(), energies.end());
        std::sort(m_selfveto_zeniths.begin(), m_selfveto_zeniths.end());

        // Set node array
        m_selfveto_energies.nodes(energies);

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
