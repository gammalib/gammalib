/***************************************************************************
 *         GCOMHkds.cpp - COMPTEL Housekeeping Data collection class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMHkds.hpp
 * @brief COMPTEL Housekeeping Data collection class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMHkds.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                         "GCOMHkds::operator[](std::string&)"
#define G_AT                                             "GCOMHkds::at(int&)"
#define G_SET1                                "GCOMHkds::set(int&, GCOMHkd&)"
#define G_SET2                        "GCOMHkds::set(std::string&, GCOMHkd&)"
#define G_INSERT                           "GCOMHkds::insert(int&, GCOMHkd&)"
#define G_REMOVE                                     "GCOMHkds::remove(int&)"
#define G_READ                                  "GCOMHkds::read(GFitsTable&)"

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
 * Constructs an empty Housekeeping Data collection.
 ***************************************************************************/
GCOMHkds::GCOMHkds(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Filename constructor
 *
 * @param[in] filename Housekeeping Data FITS file.
 *
 * Constructs a Housekeeping Data collection from a HKD FITS file.
 ***************************************************************************/
GCOMHkds::GCOMHkds(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load data
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] hkds Housekeeping Data collection.
 ***************************************************************************/
GCOMHkds::GCOMHkds(const GCOMHkds& hkds)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(hkds);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMHkds::~GCOMHkds(void)
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
 * @param[in] hkds Housekeeping Data collection.
 * @return Housekeeping Data collection.
 ***************************************************************************/
GCOMHkds& GCOMHkds::operator=(const GCOMHkds& hkds)
{
    // Execute only if object is not identical
    if (this != &hkds) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(hkds);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to Housekeeping Data container
 *
 * @param[in] name Housekeeping parameter name.
 * @return Reference to Housekeeping Data container.
 *
 * @exception GException::invalid_argument
 *            Housekeeping parameter with specified @p name not found in
 *            collection.
 *
 * Returns a reference to the Housekeeping Data container with the specified
 * @p name. An exception is thrown if the specified @p name is not found
 * in the collection.
 ***************************************************************************/
GCOMHkd& GCOMHkds::operator[](const std::string& name)
{
    // Get model index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Housekeeping parameter \""+name+"\" not found in "
                          "Housekeeping Data collection. Please specify a "
                          "valid housekeeping parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return m_hkds[index];
}


/***********************************************************************//**
 * @brief Return reference to Housekeeping Data container (const version)
 *
 * @param[in] name Housekeeping parameter name.
 * @return Reference to Housekeeping Data container.
 *
 * @exception GException::invalid_argument
 *            Housekeeping parameter with specified @p name not found in
 *            collection.
 *
 * Returns a reference to the Housekeeping Data container with the specified
 * @p name. An exception is thrown if the specified @p name is not found
 * in the collection.
 ***************************************************************************/
const GCOMHkd& GCOMHkds::operator[](const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Housekeeping parameter \""+name+"\" not found in "
                          "Housekeeping Data collection. Please specify a "
                          "valid housekeeping parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return m_hkds[index];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Housekeeping Data collection
 ***************************************************************************/
void GCOMHkds::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Housekeeping Data collection
 *
 * @return Pointer to deep copy of Housekeeping Data collection.
 ***************************************************************************/
GCOMHkds* GCOMHkds::clone(void) const
{
    return new GCOMHkds(*this);
}


/***********************************************************************//**
 * @brief Return reference to Housekeeping Data container
 *
 * @param[in] index Housekeeping Data index [0,...,size()-1].
 * @return Reference to Housekeeping Data container.
 *
 * @exception GException::out_of_range
 *            Housekeeping Data container @p index is out of range.
 *
 * Returns a reference to the Housekeeping Data container with the specified
 * @p index.
 ***************************************************************************/
GCOMHkd& GCOMHkds::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        std::string msg = "Housekeeping Data container index";
        throw GException::out_of_range(G_AT, msg, index, size());
    }

    // Return reference
    return m_hkds[index];
}


/***********************************************************************//**
 * @brief Return reference to Housekeeping Data container (const version)
 *
 * @param[in] index Housekeeping Data index [0,...,size()-1].
 * @return Reference to Housekeeping Data container.
 *
 * @exception GException::out_of_range
 *            Housekeeping Data container @p index is out of range.
 *
 * Returns a reference to the Housekeeping Data container with the specified
 * @p index.
 ***************************************************************************/
const GCOMHkd& GCOMHkds::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        std::string msg = "Housekeeping Data container index";
        throw GException::out_of_range(G_AT, msg, index, size());
    }

    // Return reference
    return m_hkds[index];
}


/***********************************************************************//**
 * @brief Set Housekeeping Data container in collection
 *
 * @param[in] index Housekeeping Data container index [0,...,size()[.
 * @param[in] hkd Housekeeping Data container.
 * @return Reference to deep copy of Housekeeping Data container.
 *
 * @exception GException::out_of_range
 *            Housekeeping Data container @p index is out of range.
 * @exception GException::invalid_value
 *            Housekeeping parameter exists already in collection at other
 *            index.
 *
 * Set Housekeeping Data container in the collection. Note that each
 * housekeeping parameter can only exist once in an collection, hence an
 * exception will be thrown if the same housekeeping parameter exists already
 * at another index. A deep copy of the Housekeeping Data container will be
 * made.
 ***************************************************************************/
GCOMHkd& GCOMHkds::set(const int& index, const GCOMHkd& hkd)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        std::string msg = "Housekeeping Data container index";
        throw GException::out_of_range(G_SET1, msg, index, size());
    }
    #endif

    // Check if a container with specified name does not yet exist
    int inx = get_index(hkd.name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set Housekeeping Data container with name \""+
            hkd.name()+"\" in collection at index "+gammalib::str(index)+
            ", but a Housekeeping Data container with the same name exists "
            "already at index "+gammalib::str(inx)+" in the collection. "
            "Every Housekeeping Data container in the collection needs a "
            "unique name.";
        throw GException::invalid_value(G_SET1, msg);
    }

    // Assign Housekeeping Data container
    m_hkds[index] = hkd;

    // Return reference to Housekeeping Data container
    return m_hkds[index];
}


/***********************************************************************//**
 * @brief Set Housekeeping Data container in collection
 *
 * @param[in] name Housekeeping parameter name.
 * @param[in] hkd Housekeeping Data container.
 * @return Reference to deep copy of Housekeeping Data container.
 *
 * @exception GException::invalid_argument
 *            Housekeeping parameter with specified @p name not found in
 *            collection.
 * @exception GException::invalid_value
 *            Name of Housekeeping parameter exists already in collection.
 *
 * Set Housekeeping Data container in the collection replacing the container
 * with the specified Housekeeping parameter @p name. Note that each
 * housekeeping parameter can only exist once in an collection, hence an
 * exception will be thrown if the same housekeeping parameter exists already
 * at another index. A deep copy of the Housekeeping Data container will be
 * made.
 ***************************************************************************/
GCOMHkd& GCOMHkds::set(const std::string& name, const GCOMHkd& hkd)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        std::string msg = "Housekeeping parameter \""+name+"\" not found in "
                          "Housekeeping Data collection. Please specify a "
                          "Housekeeping parameter name that exists in the "
                          "collection.";
        throw GException::invalid_argument(G_SET2, msg);
    }

    // Check if a container with specified name does not yet exist
    int inx = get_index(hkd.name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set Housekeeping Data container with name \""+
            hkd.name()+"\" in collection at index "+gammalib::str(index)+
            ", but a Housekeeping Data container with the same name exists "
            "already at index "+gammalib::str(inx)+" in the collection. "
            "Every Housekeeping Data container in the collection needs a "
            "unique name.";
        throw GException::invalid_value(G_SET2, msg);
    }

    // Assign Housekeeping Data container
    m_hkds[index] = hkd;

    // Return reference to Housekeeping Data container
    return m_hkds[index];
}


/***********************************************************************//**
 * @brief Append Housekeeping Data container to collection
 *
 * @param[in] hkd Housekeeping Data container.
 * @return Reference to appended Housekeeping Data container.
 *
 * Appends Housekeeping Data container to the collection by making a deep
 * copy of specified container.
 ***************************************************************************/
GCOMHkd& GCOMHkds::append(const GCOMHkd& hkd)
{
    // Append housekeeping data container to list
    m_hkds.push_back(hkd);

    // Return reference
    return m_hkds[size()-1];
}


/***********************************************************************//**
 * @brief Insert Housekeeping Data container into collection
 *
 * @param[in] index Housekeeping Data container index (0,...,size()-1).
 * @param[in] hkd Housekeeping Data container.
 * @return Reference to inserted Housekeeping Data container.
 *
 * @exception GException::out_of_range
 *            Housekeeping Data container index is out of range.
 *
 * Inserts an @p Housekeeping Data container into the collection before the
 * Housekeeping Data container with the specified @p index.
 ***************************************************************************/
GCOMHkd& GCOMHkds::insert(const int& index, const GCOMHkd& hkd)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            std::string msg = "Housekeeping Data container index";
            throw GException::out_of_range(G_INSERT, msg, index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            std::string msg = "Housekeeping Data container index";
            throw GException::out_of_range(G_INSERT, msg, index, size());
        }
    }
    #endif

    // Inserts Housekeeping Data container into collection
    m_hkds.insert(m_hkds.begin()+index, hkd);

    // Return reference
    return m_hkds[index];
}


/***********************************************************************//**
 * @brief Remove Housekeeping Data container from collection
 *
 * @param[in] index Housekeeping Data container index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Housekeeping Data container index is out of range.
 *
 * Remove Housekeeping Data container of specified @p index from collection.
 ***************************************************************************/
void GCOMHkds::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        std::string msg = "Housekeeping Data container index";
        throw GException::out_of_range(G_REMOVE, msg, index, size());
    }
    #endif

    // Erase Housekeeping Data container from collection
    m_hkds.erase(m_hkds.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if Housekeeping parameter exists in collection
 *
 * @param[in] name Housekeeping parameter name.
 * @return True if Housekeeping parameter with specified @p name exists.
 *
 * Searches the collection for a given Housekeeping parameter @p name. The
 * method returns @c true if the specified @p name was found.
 ***************************************************************************/
bool GCOMHkds::contains(const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Extend Housekeeping Data collection
 *
 * @param[in] hkds Housekeeping Data collection.
 *
 * Extend Housekeeping Data collection by extending the container for each
 * Housekeeping parameter.
 ***************************************************************************/
void GCOMHkds::extend(const GCOMHkds& hkds)
{
    // Continue only if specified Housekeeping Data collection is not empty
    if (!hkds.is_empty()) {

        // Loop over all housekeeping parameters
        for (int i = 0; i < hkds.size(); ++i) {

            // If parameter exists already in collection then extend the
            // corresponding container
            if (contains(hkds[i].name())) {
                m_hkds[i].extend(hkds[i]);
            }

            // ... otherwise append container to collection
            else {
                append(hkds[i]);
            }

        } // endfor: looped over housekeeping parameters

    } // endif: Specified Housekeeping Data collection was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Housekeeping Data collection from FITS file
 *
 * @param[in] filename COMPTEL HKD FITS file name.
 *
 * Load a Housekeeping Data collection from a HKD FITS file into the
 * collection.
 ***************************************************************************/
void GCOMHkds::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get HDU (pointer is always valid)
    const GFitsTable& hdu = *fits.table(1);

    // Read Housekeeping Data FITS table
    read(hdu);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Housekeeping Data collection from FITS table
 *
 * @param[in] table COMPTEL HKD FITS table.
 *
 * Reads a Housekeeping Data collection from a FITS table into the
 * collection.
 ***************************************************************************/
void GCOMHkds::read(const GFitsTable& table)
{
    // Clear
    clear();

    // Extract number of records in FITS table
    int num = table.nrows();

    // If there are records in FITS table then extract housekeeping data
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_tjd       = table["TJD"];       // days
        const GFitsTableCol* ptr_tics      = table["TICS"];      // ticks
        const GFitsTableCol* ptr_parameter = table["PARAMETER"];
        const GFitsTableCol* ptr_value     = table["VALUE"];

        // Loop over HKD records
        for (int i = 0; i < num; ++i) {

            // Get reference to Housekeeping Data container by either
            // returning an existing reference or by appending a new
            // housekeeping parameter to the collection
            GCOMHkd* hkd;
            int index = get_index(ptr_parameter->string(i));
            if (index == -1) {
                hkd = &(append(GCOMHkd(ptr_parameter->string(i))));
            }
            else {
                hkd = &((*this)[index]);
            }

            // Append time and value to Housekeeping Data container
            hkd->append(gammalib::com_time(ptr_tjd->integer(i), ptr_tics->integer(i)),
                        ptr_value->real(i));

        } // endfor: looped over HKD records

    } // endif: there were records in FITS table

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Housekeeping Data collection
 *
 * @param[in] chatter Chattiness.
 * @return String containing Housekeeping Data collection information.
 ***************************************************************************/
std::string GCOMHkds::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMHkds ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Parameters"));
        result.append(gammalib::str(size()));

        // Append information on housekeeping parameters
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+gammalib::parformat((*this)[i].name()));
            result.append(gammalib::str((*this)[i].size()));
        }

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
void GCOMHkds::init_members(void)
{
    // Initialise members
    m_hkds.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] hkds Housekeeping Data collection.
 ***************************************************************************/
void GCOMHkds::copy_members(const GCOMHkds& hkds)
{
    // Copy members
    m_hkds = hkds.m_hkds;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMHkds::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return Housekeeping Data container index by parameter name
 *
 * @param[in] name Housekeeping parameter name.
 * @return Housekeeping parameter index (-1 if Housekeeping parameter name
 *         has not been found)
 *
 * Returns Housekeeping Data container index based on the specified
 * housekeeping parameter @p name. If no Housekeeping Data container with the
 * specified @p name is found the method returns -1.
 ***************************************************************************/
int GCOMHkds::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search Housekeeping Data with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_hkds[i].name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
