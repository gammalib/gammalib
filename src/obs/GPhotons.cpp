/***************************************************************************
 *                   GPhotons.cpp - Photon container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GPhotons.cpp
 * @brief Photon container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GPhotons.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS                              "GPhotons::operator[](int&)"
#define G_INSERT                           "GPhotons::insert(int&, GPhoton&)"
#define G_REMOVE                                     "GPhotons::remove(int&)"

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
GPhotons::GPhotons(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param photons Photon container.
 ***************************************************************************/
GPhotons::GPhotons(const GPhotons& photons)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(photons);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPhotons::~GPhotons(void)
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
 * @param[in] photons Photon container.
 * @return Photon container.
 ***************************************************************************/
GPhotons& GPhotons::operator=(const GPhotons& photons)
{
    // Execute only if object is not identical
    if (this != &photons) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(photons);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to photon
 *
 * @param[in] index Photon index [0,...,size()[.
 * @return Photon.
 *
 * @exception GException::out_of_range
 *            Photon index is out of range.
 ***************************************************************************/
GPhoton& GPhotons::operator[](const int& index)
{
    // If index is outside boundary then throw an error
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OP_ACCESS, "Photon index",
                                       index, size());
    }
    #endif
    
    // Return reference
    return m_photons[index];
}


/***********************************************************************//**
 * @brief Return reference to photon (const version)
 *
 * @param[in] index Photon index [0,...,size()[.
 * @return Photon.
 *
 * @exception GException::out_of_range
 *            Photon index is out of range.
 ***************************************************************************/
const GPhoton& GPhotons::operator[](const int& index) const
{
    // If index is outside boundary then throw an error
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OP_ACCESS, "Photon index",
                                       index, size());
    }
    #endif

    // Return reference
    return m_photons[index];
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear container
 ***************************************************************************/
void GPhotons::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pointer to deep copy of photon container.
 ***************************************************************************/
GPhotons* GPhotons::clone(void) const
{
    // Clone this image
    return new GPhotons(*this);
}


/***********************************************************************//**
 * @brief Append photon to container
 *
 * @param[in] photon Observation.
 *
 * This method appends a photon to the container by copying it.
 ***************************************************************************/
void GPhotons::append(const GPhoton& photon)
{
    // Append photon to list
    m_photons.push_back(photon);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert photon into container
 *
 * @param[in] index Photon index [0,...,size()[.
 * @param[in] photon Photon.
 *
 * @exception GException::out_of_range
 *            Photon index is out of range.
 *
 * Inserts a @p photon into the container before the photon with the
 * specified @p index.
 ***************************************************************************/
void GPhotons::insert(const int& index, const GPhoton& photon)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Photon index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Photon index",
                                           index, size());
        }
    }
    #endif

    // Inserts photon
    m_photons.insert(m_photons.begin()+index, photon);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove photon from container
 *
 * @param[in] index Photon index [0,...,size()[.
 *
 * @exception GException::out_of_range
 *            Photon index is out of range.
 *
 * Remove photon of specified @p index from container.
 ***************************************************************************/
void GPhotons::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Photon index",
                                       index, size());
    }
    #endif

    // Erase photon component from container
    m_photons.erase(m_photons.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append photon container
 *
 * @param[in] photons Photon container.
 *
 * Append photon container to the container.
 ***************************************************************************/
void GPhotons::extend(const GPhotons& photons)
{
    // Do nothing if photon container is empty
    if (!photons.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = photons.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_photons.push_back(photons[i]);
        }

    } // endif: photon container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve memory for photons in container
 *
 * @param[in] num Number of photons.
 *
 * This method reserves memory for @p num photons in the container.
 ***************************************************************************/
void GPhotons::reserve(const int& num)
{
    // Reserve memory
    m_photons.reserve(num);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print photon container information
 *
 * @param[in] chatter Chattiness.
 * @return String containing photon container information.
 ***************************************************************************/
std::string GPhotons::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GPhotons ===");

        // Append photon container information
        result.append("\n"+gammalib::parformat("Number of photons"));
        result.append(gammalib::str(size()));

        // EXPLICIT: Append photons
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n");
                result.append(m_photons[i].print());
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
void GPhotons::init_members(void)
{
    // Initialise members
    m_photons.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] photons Photon container.
 ***************************************************************************/
void GPhotons::copy_members(const GPhotons& photons)
{
    // Copy attributes
    m_photons = photons.m_photons;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPhotons::free_members(void)
{
    // Return
    return;
}
