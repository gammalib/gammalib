/***************************************************************************
 *                        GEnergy.cpp - Energy class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GEnergy.cpp
 * @brief Energy value class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cfloat>
#include <cmath>
#include "GEnergy.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Constants __________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT                 "GEnergy::GEnergy(double&, std::string&)"
#define G_OPERATOR1              "GEnergy::operator()(double&, std::string&)"
#define G_OPERATOR2                       "GEnergy::operator()(std::string&)"
#define G_LOG10_GET                            "GEnergy::log10(std::string&)"
#define G_LOG10_SET                   "GEnergy::log10(double&, std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GEnergy::GEnergy(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] eng Energy.
 ***************************************************************************/
GEnergy::GEnergy(const GEnergy& eng)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(eng);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy constructor
 *
 * @param[in] eng Energy.
 * @param[in] unit Energy unit (one of erg(s), keV, MeV, GeV, TeV).
 *
 * Construct energy from an energy value and unit. The constructor interprets
 * the unit string and performs automatic conversion of the energy value. 
 ***************************************************************************/
GEnergy::GEnergy(const double& eng, const std::string& unit)
{ 
    // Initialise private members
    init_members();

    // Set energy according to unit string
    this->operator()(eng, unit);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEnergy::~GEnergy(void)
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
 * @param[in] eng Energy.
 * @return Energy.
 ***************************************************************************/
GEnergy& GEnergy::operator=(const GEnergy& eng)
{ 
    // Execute only if object is not identical
    if (this != &eng) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(eng);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Unit set operator
 *
 * @param[in] eng Energy.
 * @param[in] unit Energy unit (one of erg(s), keV, MeV, GeV, TeV, Angstrom).
 *
 * @exception GException::invalid_argument
 *            Invalid energy unit specified.
 *
 * Construct energy from an energy value and unit. The constructor interprets
 * the unit string and performs automatic conversion of the energy value. 
 ***************************************************************************/
void GEnergy::operator()(const double& eng, const std::string& unit)
{ 
    // Set energy according to unit string
    std::string eunit = gammalib::strip_whitespace(gammalib::tolower(unit));
    if (eunit == "erg" || eunit == "ergs") {
        this->erg(eng);
    }
    else if (eunit == "kev") {
        this->keV(eng);
    }
    else if (eunit == "mev") {
        this->MeV(eng);
    }
    else if (eunit == "gev") {
        this->GeV(eng);
    }
    else if (eunit == "tev") {
        this->TeV(eng);
    }
    else if (eunit == "angstrom") {
        this->MeV(0.012417281/eng);
    }
    else {
        throw GException::invalid_argument(G_OPERATOR1, unit,
              "Valid energy units are \"erg(s)\", \"keV\", \"MeV\","
              " \"GeV\", or \"TeV\" (case insensitive).");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unit access operator
 *
 * @param[in] unit Unit.
 * @return Energy in requested units.
 *
 * @exception GException::invalid_argument
 *            Invalid energy unit specified.
 *
 * Returns the energy in the requested units.
 ***************************************************************************/
double GEnergy::operator()(const std::string& unit) const
{ 
    // Initialise energy
    double energy = 0.0;

    // Set energy according to unit string
    std::string eunit = gammalib::tolower(unit);
    if (eunit == "erg" || eunit == "ergs") {
        energy = this->erg();
    }
    else if (eunit == "kev") {
        energy = this->keV();
    }
    else if (eunit == "mev") {
        energy = this->MeV();
    }
    else if (eunit == "gev") {
        energy = this->GeV();
    }
    else if (eunit == "tev") {
        energy = this->TeV();
    }
    else {
        throw GException::invalid_argument(G_OPERATOR2, unit,
              "Valid energy units are \"erg(s)\", \"keV\", \"MeV\","
              " \"GeV\", or \"TeV\" (case insensitive).");
    }
  
    // Return energy
    return energy;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GEnergy::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pointer to deep copy of energy.
 ***************************************************************************/
GEnergy* GEnergy::clone(void) const
{
    // Clone this image
    return new GEnergy(*this);
}


/***********************************************************************//**
 * @brief Return energy in erg
 *
 * @return Energy in erg.
 ***************************************************************************/
double GEnergy::erg(void) const
{
    // Compute energy
    double energy = m_energy * gammalib::MeV2erg;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Return energy in keV
 *
 * @return Energy in keV.
 ***************************************************************************/
double GEnergy::keV(void) const
{
    // Compute energy
    double energy = m_energy * 1.0e+3;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Return energy in MeV
 *
 * @return Energy in MeV.
 ***************************************************************************/
double GEnergy::MeV(void) const
{
    // Return energy
    return m_energy;
}


/***********************************************************************//**
 * @brief Return energy in GeV
 *
 * @return Energy in GeV.
 ***************************************************************************/
double GEnergy::GeV(void) const
{
    // Compute energy
    double energy = m_energy * 1.0e-3;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Return energy in TeV
 *
 * @return Energy in TeV.
 ***************************************************************************/
double GEnergy::TeV(void) const
{
    // Compute energy
    double energy = m_energy * 1.0e-6;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Return log10 of energy in erg
 *
 * @return Energy in log10 erg.
 *
 * Returns the log10 of the energy in erg.
 ***************************************************************************/
double GEnergy::log10erg(void) const
{
    // Set offset
    const double offset = std::log10(gammalib::MeV2erg);

    // Return log10 energy
    return (log10MeV()+offset); 
}


/***********************************************************************//**
 * @brief Return log10 of energy in keV
 *
 * @return Energy in log10 keV.
 *
 * Returns the log10 of the energy in keV.
 ***************************************************************************/
double GEnergy::log10keV(void) const
{
    // Return log10 energy
    return (log10MeV()+3.0); 
}


/***********************************************************************//**
 * @brief Return log10 of energy in MeV
 *
 * @return Energy in log10 MeV.
 *
 * Returns the log10 of the energy in MeV. The result is stored internally
 * and not recomputed when the method is called again with the same energy
 * value. This speeds up computation. In case that the energy is not positive
 * the method returns DBL_MIN.
 ***************************************************************************/
double GEnergy::log10MeV(void) const
{
    // If required compute log10 of energy.
    if (!m_has_log10) {
        m_elog10    = (m_energy > 0.0) ? std::log10(m_energy) : DBL_MIN;
        m_has_log10 = true;
    }
    
    // Return log10 energy
    return m_elog10; 
}


/***********************************************************************//**
 * @brief Return log10 of energy in GeV
 *
 * @return Energy in log10 GeV.
 *
 * Returns the log10 of the energy in GeV.
 ***************************************************************************/
double GEnergy::log10GeV(void) const
{
    // Return log10 energy
    return (log10MeV()-3.0); 
}


/***********************************************************************//**
 * @brief Return log10 of energy in TeV
 *
 * @return Energy in log10 TeV.
 *
 * Returns the log10 of the energy in TeV.
 ***************************************************************************/
double GEnergy::log10TeV(void) const
{
    // Return log10 energy
    return (log10MeV()-6.0); 
}


/***********************************************************************//**
 * @brief Set log10 of energy with unit specification
 *
 * @param[in] unit Unit of log10 of energy.
 * @return eng log10 of energy.
 *
 * @exception GException::invalid_argument
 *            Unit argument is not valid.
 ***************************************************************************/
double GEnergy::log10(const std::string& unit) const
{
    // Initialise result
    double logE = 0.0;

    // Set energy according to unit string
    std::string eunit = gammalib::tolower(unit);
    if (eunit == "erg" || eunit == "ergs") {
        logE = this->log10erg();
    }
    else if (eunit == "kev") {
        logE = this->log10keV();
    }
    else if (eunit == "mev") {
        logE = this->log10MeV();
    }
    else if (eunit == "gev") {
        logE = this->log10GeV();
    }
    else if (eunit == "tev") {
        logE = this->log10TeV();
    }
    else {
        throw GException::invalid_argument(G_LOG10_GET, unit,
              "Valid energy units are \"erg(s)\", \"keV\", \"MeV\","
              " \"GeV\", or \"TeV\" (case insensitive).");
    }
    
    // Return
    return logE;
}


/***********************************************************************//**
 * @brief Set energy in erg
 *
 * @param[in] eng Energy in erg.
 ***************************************************************************/
void GEnergy::erg(const double& eng)
{
    // Set energy
    m_energy = eng * gammalib::erg2MeV;

    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in keV
 *
 * @param[in] eng Energy in keV.
 ***************************************************************************/
void GEnergy::keV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e-3;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in MeV
 *
 * @param[in] eng Energy in MeV.
 ***************************************************************************/
void GEnergy::MeV(const double& eng)
{
    // Set energy
    m_energy = eng;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in GeV
 *
 * @param[in] eng Energy in GeV.
 ***************************************************************************/
void GEnergy::GeV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e+3;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in TeV
 *
 * @param[in] eng Energy in TeV.
 ***************************************************************************/
void GEnergy::TeV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e+6;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in erg
 *
 * @param[in] eng log10 of energy in erg.
 ***************************************************************************/
void GEnergy::log10erg(const double& eng)
{
    // Set offset
    const double offset = std::log10(gammalib::MeV2erg);

    // Set energy
    log10MeV(eng-offset);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in keV
 *
 * @param[in] eng log10 of energy in keV.
 ***************************************************************************/
void GEnergy::log10keV(const double& eng)
{
    // Set energy
    log10MeV(eng-3.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in MeV
 *
 * @param[in] eng log10 of energy in MeV.
 ***************************************************************************/
void GEnergy::log10MeV(const double& eng)
{
    // Set energy
    m_elog10    = eng;
    m_energy    = std::pow(10.0, eng);
    m_has_log10 = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in GeV
 *
 * @param[in] eng log10 of energy in GeV.
 ***************************************************************************/
void GEnergy::log10GeV(const double& eng)
{
    // Set energy
    log10MeV(eng+3.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in TeV
 *
 * @param[in] eng log10 of energy in TeV.
 ***************************************************************************/
void GEnergy::log10TeV(const double& eng)
{
    // Set energy
    log10MeV(eng+6.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy with unit specification
 *
 * @param[in] eng log10 of energy.
 * @param[in] unit Unit of log10 of energy.
 *
 * @exception GException::invalid_argument
 *            Unit argument is not valid.
 ***************************************************************************/
void GEnergy::log10(const double& eng, const std::string& unit)
{
    // Set energy according to unit string
    std::string eunit = gammalib::tolower(unit);
    if (eunit == "erg" || eunit == "ergs") {
        this->log10erg(eng);
    }
    else if (eunit == "kev") {
        this->log10keV(eng);
    }
    else if (eunit == "mev") {
        this->log10MeV(eng);
    }
    else if (eunit == "gev") {
        this->log10GeV(eng);
    }
    else if (eunit == "tev") {
        this->log10TeV(eng);
    }
    else {
        throw GException::invalid_argument(G_LOG10_SET, unit,
              "Valid energy units are \"erg(s)\", \"keV\", \"MeV\","
              " \"GeV\", or \"TeV\" (case insensitive).");
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print energy
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy information.
 ***************************************************************************/
std::string GEnergy::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append energy
        if (GeV() > 1000.0) {
            result.append(gammalib::str(TeV())+" TeV");
        }
        else if (MeV() > 1000.0) {
            result.append(gammalib::str(GeV())+" GeV");
        }
        else if (keV() > 1000.0) {
            result.append(gammalib::str(MeV())+" MeV");
        }
        else {
            result.append(gammalib::str(keV())+" keV");
        }

        // VERBOSE: append energy and log10 energy
        if (chatter == VERBOSE) {
            result.append(" (E="+gammalib::str(m_energy));
            if (m_has_log10) {
                result.append(", log10(E)="+gammalib::str(m_elog10)+")");
            }
            else {
                result.append(", no log10(E) value)");
            }
        }

    } // endif: chatter was not silent

    // Return
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
void GEnergy::init_members(void)
{
    // Initialise members
    m_energy    = 0.0;
    m_elog10    = DBL_MIN;
    m_has_log10 = false;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] eng Energy.
 ***************************************************************************/
void GEnergy::copy_members(const GEnergy& eng)
{
    // Copy time
    m_energy    = eng.m_energy;
    m_elog10    = eng.m_elog10;
    m_has_log10 = eng.m_has_log10;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEnergy::free_members(void)
{
    // Return
    return;
}
