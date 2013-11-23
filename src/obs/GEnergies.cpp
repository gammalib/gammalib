/***************************************************************************
 *                  GEnergies.cpp - Energy container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
GEnergies::GEnergies(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename FITS filename.
 * @param[in] extname Energies extension name (defaults to "ENERGIES")
 *
 * Constructs energy container by loading energies from FITS file.
 ***************************************************************************/
GEnergies::GEnergies(const std::string& filename, const std::string& extname)
{
    // Initialise members
    init_members();

    // Load energies
    load(filename, extname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param energies Energy container.
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
 * @brief Clear container
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
 * @brief Clone instance
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
    if (isempty()) {
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
    if (!energies.isempty()) {

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
 * @brief Load energies from FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] extname Energies extension name (defaults to "ENERGIES")
 *
 * Loads the energies from FITS file.
 ***************************************************************************/
void GEnergies::load(const std::string& filename, const std::string& extname)
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get energies HDU
    GFitsTable* hdu = file.table(extname);

    // Read energies from HDU
    read(hdu);

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
 * @param[in] extname Energies extension name (defaults to "ENERGIES")
 *
 * Saves energies into extension @p extname of a FITS file. If the file does
 * not exist it is created. If the file exists the energies are appended as
 * extension. If another energies extension exists already it is overwritten
 * if @p clobber=true.
 ***************************************************************************/
void GEnergies::save(const std::string& filename, bool clobber,
                     const std::string& extname) const
{
    // Allocate FITS file
    GFits file;

    // Write energies to FITS file
    write(&file, extname);

    // Save to file
    file.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energies from FITS table
 *
 * @param[in] hdu Pointer to FITS table.
 *
 * Reads the energies from a FITS table.
 ***************************************************************************/
void GEnergies::read(const GFitsTable* hdu)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract number of energy bins in FITS file
        int num = hdu->integer("NAXIS2");

        // Continue only if there are energy bins
        if (num > 0) {

            // Get the unit of the energies. If no TUNIT1 header keyword is
            // found then use MeV
            std::string unit = "MeV";
            if (hdu->hascard("TUNIT1")) {
                unit = hdu->string("TUNIT1");
            }

            // Get the column with the name "Energy"
            const GFitsTableCol* col_energy = (*hdu)["Energy"];

            // Set energies
            for (int i = 0; i < num; ++i) {
                append(GEnergy(col_energy->real(i), unit));
            }

        } // endif: there were energy bins

    } // endif: the HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energies into FITS object
 *
 * @param[in] file Pointer to FITS file.
 * @param[in] extname Energy extension name (defaults to "ENERGIES")
 *
 * Writes energies into FITS object.
 ***************************************************************************/
void GEnergies::write(GFits* file, const std::string& extname) const
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
    file->append(*table);

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
