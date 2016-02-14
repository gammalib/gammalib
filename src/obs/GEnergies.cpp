/***************************************************************************
 *                  GEnergies.cpp - Energy container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Juergen Knoedlseder                         *
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
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                            "GEnergies::at(int&)"
#define G_INSERT                          "GEnergies::insert(int&, GEnergy&)"
#define G_REMOVE                                    "GEnergies::remove(int&)"
#define G_SET_LIN              "GEnergies::set_lin(int&, GEnergy&, GEnergy&)"
#define G_SET_LOG              "GEnergies::set_log(int&, GEnergy&, GEnergy&)"

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
 * @brief Copy constructor
 *
 * @param energies Energy container.
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
 * @param[in] log Use logarithmic spacing? (defaults to true).
 *
 * Constructs energy container by defining @p num energies between @p emin
 * and @p emax. The @p log parameter controls whether the energy spacing is
 * logarihmic (default) or linear.
 ***************************************************************************/
GEnergies::GEnergies(const int& num, const GEnergy& emin, const GEnergy& emax,
                     const bool& log)
{
    // Initialise members
    init_members();

    // Set intervals
    if (log) {
        this->set_log(num, emin, emax);
    }
    else {
        this->set_lin(num, emin, emax);
    }

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
 * @brief Set linearly spaced energies
 *
 * @param[in] num Number of energies.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * @exception GException::invalid_argument
 *            Invalid number of energies or minimum and maximum energy
 *            requested.
 *
 * Creates @p num linearly spaced energies running from @p emin to @p emax.
 ***************************************************************************/
void GEnergies::set_lin(const int&     num,
                        const GEnergy& emin,
                        const GEnergy& emax)
{
    // Initialise members
    clear();

    // Check validity of number of energies
    if (num < 0) {
        std::string msg = "Negative number of energies "+gammalib::str(num)+
                          " specified. Specify a non-negative number of "
                          "energies.";
        throw GException::invalid_argument(G_SET_LIN, msg);
    }

    // Check validity of energies
    if (emin > emax) {
        std::string msg = "Minimum energy "+emin.print()+" is larger than "
                          "maximum energy "+emax.print()+". Specify "
                          "a minimum energy that does not exceed the "
                          "maximum energy.";
        throw GException::invalid_argument(G_SET_LIN, msg);
    }

    // Case A: we have 1 energy
    if (num == 1) {

        // Check that emin and emax are equal
        if (emin != emax) {
            std::string msg = "Single energy is requested but the minimum "
                              "energy "+emin.print()+" differs from the "
                              "maximum energy "+emax.print()+". Specify "
                              "identical energies.";
            throw GException::invalid_argument(G_SET_LIN, msg);
        }

        // Append energy
        append(emin);

    }

    // Case B: more than 1 energy requested
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
 * @param[in] num Number of energies.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * @exception GException::invalid_argument
 *            Invalid number of energies or minimum and maximum energy
 *            requested.
 *
 * Creates @p num logarithmically spaced energies running from @p emin to
 * @p emax.
 ***************************************************************************/
void GEnergies::set_log(const int&     num,
                        const GEnergy& emin,
                        const GEnergy& emax)
{
    // Initialise members
    clear();

    // Check validity of number of energies
    if (num < 0) {
        std::string msg = "Negative number of energies "+gammalib::str(num)+
                          " specified. Specify a non-negative number of "
                          "energies.";
        throw GException::invalid_argument(G_SET_LOG, msg);
    }

    // Check validity of energies
    if (emin > emax) {
        std::string msg = "Minimum energy "+emin.print()+" is larger than "
                          "maximum energy "+emax.print()+". Specify "
                          "a minimum energy that does not exceed the "
                          "maximum energy.";
        throw GException::invalid_argument(G_SET_LOG, msg);
    }

    // Case A: we have 1 energy
    if (num == 1) {

        // Check that emin and emax are equal
        if (emin != emax) {
            std::string msg = "Single energy is requested but the minimum "
                              "energy "+emin.print()+" differs from the "
                              "maximum energy "+emax.print()+". Specify "
                              "identical energies.";
            throw GException::invalid_argument(G_SET_LOG, msg);
        }

        // Append energy
        append(emin);

    }

    // Case B: more than 1 energy requested
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
 * @brief Load energies from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the energies from FITS file.
 *
 * If no extension name is provided, the energies are loaded from the
 * "ENERGIES" extension.
 ***************************************************************************/
void GEnergies::load(const GFilename& filename)
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get energies table
    const GFitsTable& table = *file.table(filename.extname("ENERGIES"));

    // Read energies from table
    read(table);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energies to FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite any existing energies extension?
 *
 * Saves energiesinto a FITS file. If a file with the given @p filename does
 * not yet exist it will be created, otherwise the method opens the existing
 * file. The method will create (or replace an existing) energies extension.
 * The extension name can be specified as part of the @p filename, or if no
 * extension name is given, is assumed to be "ENERGIES".
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GEnergies::save(const GFilename& filename, const bool& clobber) const
{
    // Open or create FITS file
    GFits fits(filename, true);

    // Write energies to FITS file
    write(fits, filename.extname("ENERGIES"));

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
 * @param[in] file FITS file.
 * @param[in] extname Energy extension name (default: "ENERGIES")
 *
 * Writes energies into FITS object.
 ***************************************************************************/
void GEnergies::write(GFits& file, const std::string& extname) const
{
    // Set number of energies
    int num = m_energies.size();

    // Create Energy column
    GFitsTableDoubleCol col_energy = GFitsTableDoubleCol("Energy", num);

    // Fill energy column in units of MeV
    for (int i = 0; i < num; ++i) {
        col_energy(i) = (*this)[i].MeV();
    }
    col_energy.unit("MeV");

    // Create energies table
    GFitsBinTable* table = new GFitsBinTable(num);
    table->append(col_energy);
    table->extname(extname);

    // Write to FITS file
    file.append(*table);

    // Free table
    delete table;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print energy container information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
                result.append("\n"+m_energies[i].print(chatter));
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
