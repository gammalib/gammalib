/***************************************************************************
 *                  GEnergies.cpp - Energy container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2020 by Juergen Knoedlseder                         *
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
 * @file GEnergies.cpp
 * @brief Energy container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GEnergies.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                            "GEnergies::at(int&)"
#define G_INSERT                          "GEnergies::insert(int&, GEnergy&)"
#define G_REMOVE                                    "GEnergies::remove(int&)"
#define G_SET       "GEnergies::set(int&, GEnergy&, GEnergy&, std::string&, "\
                                                                   "double&)"
#define G_SET_LIN              "GEnergies::set_lin(int&, GEnergy&, GEnergy&)"
#define G_SET_LOG              "GEnergies::set_log(int&, GEnergy&, GEnergy&)"
#define G_SET_POW     "GEnergies::set_pow(int&, GEnergy&, GEnergy&, double&)"

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
 *
 * Constructs empty energy container.
 ***************************************************************************/
GEnergies::GEnergies(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS file constructor
 *
 * @param[in] filename FITS file name.
 *
 * Constructs energy container from a FITS file.
 ***************************************************************************/
GEnergies::GEnergies(const GFilename& filename)
{
    // Initialise members
    init_members();

    // Load energies
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy boundaries constructor
 *
 * @param[in] ebounds Energy boundaries.
 *
 * Constructs energy container from energy boundaries.
 ***************************************************************************/
GEnergies::GEnergies(const GEbounds& ebounds)
{
    // Initialise members
    init_members();

    // Set energies from energy boundaries
    set(ebounds);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param energies[in] Energy container.
 *
 * Construct energy container by copying from another energy container.
 ***************************************************************************/
GEnergies::GEnergies(const GEnergies& energies)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(energies);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Interval constructor
 *
 * @param[in] num Number of energies.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] method Energy spacing method (one of "LIN", "LOG" or "POW").
 * @param[in] gamma Power law index for @p POW method.
 *
 * Constructs energy container by defining @p num energies between @p emin
 * and @p emax. The @p method parameter controls the energy spacing of the
 * energies. See the set() method for more information.
 ***************************************************************************/
GEnergies::GEnergies(const int&         num,
                     const GEnergy&     emin,
                     const GEnergy&     emax,
                     const std::string& method,
                     const double&      gamma)
{
    // Initialise members
    init_members();

    // Set intervals
    set(num, emin, emax, method, gamma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEnergies::~GEnergies(void)
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
 * @param[in] energies Energy container.
 * @return Energy container.
 ***************************************************************************/
GEnergies& GEnergies::operator=(const GEnergies& energies)
{
    // Execute only if object is not identical
    if (this != &energies) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(energies);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear energy container
 *
 * Removes all energies from the container.
 ***************************************************************************/
void GEnergies::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone energy container
 *
 * @return Pointer to deep copy of energy container
 *
 * Makes a deep copy of the energy container instance.
 ***************************************************************************/
GEnergies* GEnergies::clone(void) const
{
    return new GEnergies(*this);
}


/***********************************************************************//**
 * @brief Return reference to energy
 *
 * @param[in] index Energy index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Returns a reference to the energy with the specified @p index.
 ***************************************************************************/
GEnergy& GEnergies::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return m_energies[index];
}


/***********************************************************************//**
 * @brief Return reference to energy (const version)
 *
 * @param[in] index Energy index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Returns a reference to the energy with the specified @p index.
 ***************************************************************************/
const GEnergy& GEnergies::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return m_energies[index];
}


/***********************************************************************//**
 * @brief Append energy to container
 *
 * @param[in] energy Energy.
 * @return Reference to appended energy.
 *
 * Appends energy to the container by making a deep copy of the energy.
 ***************************************************************************/
GEnergy& GEnergies::append(const GEnergy& energy)
{
    // Append energy to list
    m_energies.push_back(energy);

    // Return reference
    return m_energies[size()-1];
}


/***********************************************************************//**
 * @brief Insert energy into container
 *
 * @param[in] index Energy index (0,...,size()-1).
 * @param[in] energy Energy.
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Inserts an @p energy into the container before the energy with the
 * specified @p index.
 ***************************************************************************/
GEnergy& GEnergies::insert(const int& index, const GEnergy& energy)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, index, 0, size()-1);
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, index, 0, size()-1);
        }
    }
    #endif

    // Inserts energy
    m_energies.insert(m_energies.begin()+index, energy);

    // Return reference
    return m_energies[index];
}


/***********************************************************************//**
 * @brief Remove energy from container
 *
 * @param[in] index Energy index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Remove energy of specified @p index from container.
 ***************************************************************************/
void GEnergies::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, index, 0, size()-1);
    }
    #endif

    // Erase energy from container
    m_energies.erase(m_energies.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append energy container
 *
 * @param[in] energies Energy container.
 *
 * Append energy container to the container.
 ***************************************************************************/
void GEnergies::extend(const GEnergies& energies)
{
    // Do nothing if energy container is empty
    if (!energies.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = energies.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_energies.push_back(energies[i]);
        }

    } // endif: energy container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energies from energy boundaries
 *
 * @param[in] ebounds Energy boundaries.
 *
 * Sets the energies from energy boundaries. Each unique minimum and maximum
 * energy boundary will be appended as energy to the container.
 ***************************************************************************/
void GEnergies::set(const GEbounds& ebounds)
{
    // Initialise members
    clear();

    // Get number of energy boundaries
    int num = ebounds.size();

    // Append energy boundaries
    for (int i = 0; i < num; ++i) {

        // Loop over minimum and maximum boundaries
        for (int j = 0; j < 2; ++j) {

            // Get energy
            GEnergy energy = (j == 0) ? ebounds.emin(i) : ebounds.emax(i);

            // Append energy if it does not yet exist in the container
            bool    not_found = true;
            for (int k = 0; k < m_energies.size(); ++k) {
                if (energy == m_energies[k]) {
                    not_found = false;
                    break;
                }
            }
            if (not_found) {
                m_energies.push_back(energy);
            }

        } // endfor: looped over minimum and maximum boundaries

    } // endfor: looped over energy boundaries

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energies
 *
 * @param[in] num Number of energies.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] method Energy spacing method (one of "LIN", "LOG" or "POW").
 * @param[in] gamma Power law index for @p POW method.
 *
 * @exception GException::invalid_argument
 *            Invalid number of energies, energy interval or energy spacing
 *            method specified.
 *
 * Sets energies by defining @p num successive energies between @p emin and
 * @p emax. The @p method parameter controls the energy spacing. See the
 * set_lin(), set_log() and set_pow() methods for more information.
 ***************************************************************************/
void GEnergies::set(const int&         num,
                    const GEnergy&     emin,
                    const GEnergy&     emax,
                    const std::string& method,
                    const double&      gamma)
{
    // Check validity of number of energies
    if (num < 0) {
        std::string msg = "Negative number of energies "+gammalib::str(num)+
                          " specified. Please specify a non-negative number "
                          "of energies.";
        throw GException::invalid_argument(G_SET, msg);
    }

    // Check validity of energies
    if (emin > emax) {
        std::string msg = "Minimum energy "+emin.print()+" is larger than "
                          "maximum energy "+emax.print()+". Please specify "
                          "a minimum energy that does not exceed the "
                          "maximum energy.";
        throw GException::invalid_argument(G_SET, msg);
    }

    // If one energy is requested then check that emin and emax are equal
    if (num == 1) {

        // Check that emin and emax are equal
        if (emin != emax) {
            std::string msg = "Single energy is requested but the minimum "
                              "energy "+emin.print()+" differs from the "
                              "maximum energy "+emax.print()+". Please "
                              "specify identical energies.";
            throw GException::invalid_argument(G_SET, msg);
        }

    }

    // Convert method to upper-case string
    std::string umethod = gammalib::toupper(method);

    // Dispatch to corresponding set methos
    if (umethod == "LIN") {
        set_lin(num, emin, emax);
    }
    else if (umethod == "LOG") {
        set_log(num, emin, emax);
    }
    else if (umethod == "POW") {
        set_pow(num, emin, emax, gamma);
    }
    else {
        std::string msg = "Invalid energy spacing method \""+umethod+"\" "
                          "specified. Please provide one of \"LIN\", \"LOG\""
                          ", or \"POW\".";
        throw GException::invalid_argument(G_SET, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energies from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the energies from FITS file.
 *
 * If no extension name is provided, the energies are loaded from the
 * `ENERGIES` extension.
 ***************************************************************************/
void GEnergies::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get energies table
    const GFitsTable& table = *fits.table(filename.extname(gammalib::extname_energies));

    // Read energies from table
    read(table);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energies into FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite an existing energies extension?
 *
 * Saves energies into a FITS file. If a file with the given @p filename
 * does not yet exist it will be created, otherwise the method opens the
 * existing file. Energies can only be appended to an existing file if the
 * @p clobber flag is set to `true` (otherwise an exception is thrown).
 *
 * The method will append a binary FITS table containing the energies to the
 * FITS file. The extension name can be specified as part of the
 * @p filename. For example the @p filename
 *
 *      myfile.fits[ENERGY VALUES]
 *
 * will save the energies in the `ENERGY VALUES` extension of the
 * `myfile.fits` file. If the extension exists already in the file it will be
 * replaced, otherwise a new extension will be created. If no extension name
 * is provided, the method will use `ENERGIES` as the default extension name
 * for energies.
 ***************************************************************************/
void GEnergies::save(const GFilename& filename, const bool& clobber) const
{
    // Open or create FITS file (without extension name since the requested
    // extension may not yet exist in the file)
    GFits fits(filename.url(), true);

    // Write energies to FITS file
    write(fits, filename.extname(gammalib::extname_energies));

    // Save to file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energies from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads the energies from a FITS table.
 ***************************************************************************/
void GEnergies::read(const GFitsTable& table)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Extract number of energy bins in FITS file
    int num = table.integer("NAXIS2");

    // Continue only if there are energy bins
    if (num > 0) {

        // Get the unit of the energies. If no TUNIT1 header keyword is
        // found then use MeV
        std::string unit = "MeV";
        if (table.has_card("TUNIT1")) {
            unit = table.string("TUNIT1");
        }

        // Get the column with the name "Energy"
        const GFitsTableCol* col_energy = table["Energy"];

        // Set energies
        for (int i = 0; i < num; ++i) {
            append(GEnergy(col_energy->real(i), unit));
        }

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energies into FITS object
 *
 * @param[in] fits FITS file.
 * @param[in] extname Energy extension name.
 *
 * Writes energies into FITS object.
 ***************************************************************************/
void GEnergies::write(GFits& fits, const std::string& extname) const
{
    // Set number of energies
    int num = m_energies.size();

    // Create Energy column
    GFitsTableDoubleCol col_energy("Energy", num);

    // Fill energy column in units of MeV
    for (int i = 0; i < num; ++i) {
        col_energy(i) = (*this)[i].MeV();
    }
    col_energy.unit("MeV");

    // Create energies table
    GFitsBinTable table(num);
    table.append(col_energy);
    table.extname(extname);

    // If the FITS object contains already an extension with the same
    // name then remove now this extension
    if (fits.contains(extname)) {
        fits.remove(extname);
    }

    // Append energies table to FITS file
    fits.append(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print energy container information
 *
 * @param[in] chatter Chattiness.
 * @return String containing energy container information.
 ***************************************************************************/
std::string GEnergies::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GEnergies ===");

        // Append energy container information
        result.append("\n"+gammalib::parformat("Number of energies"));
        result.append(gammalib::str(size()));

        // EXPLICIT: Append energies
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n");
                result.append(gammalib::parformat("Energy "+gammalib::str(i)));
                result.append(m_energies[i].print(chatter));
            }
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEnergies::init_members(void)
{
    // Initialise members
    m_energies.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] energies Energy container.
 ***************************************************************************/
void GEnergies::copy_members(const GEnergies& energies)
{
    // Copy attributes
    m_energies = energies.m_energies;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEnergies::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set linearly spaced energies
 *
 * @param[in] num Number of energies (>0).
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * Creates \f$N\f$ linearly spaced energy boundaries running from
 * \f$E_{\rm min}\f$ to \f$E_{\rm max}\f$ defined by
 *
 * \f[
 *    E_i = E_{\rm min} + \frac{E_{\rm max} - E_{\rm min}}{N-1}
 *                        \times i
 * \f]
 *
 * The method does not check the validity of the @p num, @p emin and @p emax
 * arguments.
 ***************************************************************************/
void GEnergies::set_lin(const int&     num,
                        const GEnergy& emin,
                        const GEnergy& emax)
{
    // Initialise members
    clear();

    // If a single energy is requested then just append the minimum energy
    if (num == 1) {
        append(emin);
    }

    // ... otherwise append linearly spaced energies
    else if (num > 1) {

        // Compute bin width
        GEnergy ebin = (emax - emin)/double(num-1); 

        // Append energies
        GEnergy energy = emin;
        for (int i = 0; i < num; ++i) {
            append(energy);
            energy += ebin;
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set logarithmically spaced energies
 *
 * @param[in] num Number of energies (>0).
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * @exception GException::invalid_argument
 *            Minimum or maximum energy are not positive.
 *
 * Creates \f$N\f$ logarithmically spaced energy boundaries running from
 * \f$E_{\rm min}\f$ to \f$E_{\rm max}\f$ defined by
 *
 * \f[
 *    E_i = 10^{\log_{10} E_{\rm min} +
 *              \frac{\log_{10} E_{\rm max} - \log_{10} E_{\rm min}}{N-1}
 *              \times i}
 * \f]
 *
 * The method does not check the validity of the @p num, @p emin and @p emax
 * arguments except for the positivity of @p emin and @p emax.
 ***************************************************************************/
void GEnergies::set_log(const int&     num,
                        const GEnergy& emin,
                        const GEnergy& emax)
{
    // Initialise members
    clear();

    // Throw an exception if the minimum or maximum energy is not positive
    if (emin.MeV() <= 0.0) {
        std::string msg = "Non-positive minimum energy "+emin.print()+
                          " specified. Please provide a positive minimum "
                          "energy value.";
        throw GException::invalid_argument(G_SET_LOG, msg);
    }
    if (emax.MeV() <= 0.0) {
        std::string msg = "Non-positive maximum energy "+emax.print()+
                          " specified. Please provide a positive minimum "
                          "energy value.";
        throw GException::invalid_argument(G_SET_LOG, msg);
    }

    // If a single energy is requested then just append the minimum energy
    if (num == 1) {
        append(emin);
    }

    // ... otherwise append logarithmically spaced energies
    else if (num > 1) {

        // Compute bin width
        double elogmin = std::log10(emin.MeV());
        double elogmax = std::log10(emax.MeV());
        double elogbin = (elogmax - elogmin)/double(num-1);

        // Append energies
        GEnergy energy;
        for (int i = 0; i < num; ++i) {
            energy.MeV(std::pow(10.0, double(i)*elogbin + elogmin));
            append(energy);
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set power-law spaced energy intervals
 *
 * @param[in] num Number of energy intervals (>0).
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 * @param[in] gamma Power law index.
 *
 * @exception GException::invalid_argument
 *            Minimum or maximum energy are not positive.
 *
 * Creates \f$N\f$ power-law spaced energy boundaries running from
 * \f$E_{\rm min}\f$ to \f$E_{\rm max}\f$ defined by
 *
 * \f[
 *    E_i = \left \{
 *    \begin{array}{l l}
 *      \exp \left( \frac{\log E_{\rm max} - \log E_{\rm min}}{N-1} +
 *           \log E_{i-1} \right)
 *      & \mbox{if $a = 0$} \\
 *      \\
 *      \exp \left( \frac{\log \left(\frac{E_{\rm max}^a - E_{\rm min}^a}{N-1} \right) +
 *           E_{i-1}^{a}}{a} \right)
 *      & \mbox{else}
 *  \end{array}
 *  \right .
 * \f]
 *
 * where \f$i>0\f$, \f$a=1-\gamma\f$ and \f$E_0=E_{\rm min}\f$.
 *
 * The method does not check the validity of the @p num, @p emin and @p emax
 * arguments except for the positivity of @p emin and @p emax.
 ***************************************************************************/
void GEnergies::set_pow(const int&     num,
                        const GEnergy& emin,
                        const GEnergy& emax,
                        const double&  gamma)
{
    // Initialise members
    clear();

    // Throw an exception if the minimum or maximum energy is not positive
    if (emin.MeV() <= 0.0) {
        std::string msg = "Non-positive minimum energy "+emin.print()+
                          " specified. Please provide a positive minimum "
                          "energy value.";
        throw GException::invalid_argument(G_SET_POW, msg);
    }
    if (emax.MeV() <= 0.0) {
        std::string msg = "Non-positive maximum energy "+emax.print()+
                          " specified. Please provide a positive minimum "
                          "energy value.";
        throw GException::invalid_argument(G_SET_POW, msg);
    }

    // If a single energy is requested then just append the minimum energy
    if (num == 1) {
        append(emin);
    }

    // ... otherwise append logarithmically spaced energies
    else if (num > 1) {

        // Append first energy
        append(emin);

        // Precomputation (we work in MeV from now on)
        double a = 1.0 - gamma;
        double c = (a == 0.0) ? 1.0 / (std::log(emax.MeV())   - std::log(emin.MeV()))
                              : a   / (std::pow(emax.MeV(),a) - std::pow(emin.MeV(),a));
        double b = double(c*(num-1.0));

        // Initialse first energy
        double e = emin.MeV();

        // Loop over energies
        for (int i = 0; i < num-1; ++i) {

            // Compute next energy
            double log_e_next = (a == 0.0)
                                ? 1.0/b + std::log(e)
                                : std::log(a/b + std::pow(e,a)) / a;
            double e_next = std::exp(log_e_next);

            // Append next energy
            append(GEnergy(e_next, "MeV"));

            // Set next energy as starting energy for next iteration
            e = e_next;

        } // endfor: looped over energies

    } // endif: number of energies was positive

    // Return
    return;
}
