/***************************************************************************
 *        GCOMHkds.hpp - COMPTEL Housekeeping Data collection class        *
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
 * @brief COMPTEL Housekeeping Data collection class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMHKDS_HPP
#define GCOMHKDS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCOMHkd.hpp"
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFitsTable;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMHkds
 *
 * @brief COMPTEL Housekeeping Data collection class
 *
 * The COMPTEL Housekeeping Data collection class holds Housekeeping
 * Data containers for various parameters.
 ***************************************************************************/
class GCOMHkds : public GContainer {

public:
    // Constructors and destructors
    GCOMHkds(void);
    explicit GCOMHkds(const GFilename& filename);
    GCOMHkds(const GCOMHkds& hkds);
    virtual ~GCOMHkds(void);

    // Operators
    GCOMHkds&      operator=(const GCOMHkds& hkds);
    GCOMHkd&       operator[](const int& index);
    const GCOMHkd& operator[](const int& index) const;
    GCOMHkd&       operator[](const std::string& name);
    const GCOMHkd& operator[](const std::string& name) const;

    // Methods
    void           clear(void);
    GCOMHkds*      clone(void) const;
    std::string    classname(void) const;
    GCOMHkd&       at(const int& index);
    const GCOMHkd& at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GCOMHkd&       set(const int& index, const GCOMHkd& hkd);
    GCOMHkd&       set(const std::string& name, const GCOMHkd& hkd);
    GCOMHkd&       append(const GCOMHkd& hkd);
    GCOMHkd&       insert(const int& index, const GCOMHkd& hkd);
    void           remove(const int& index);
    void           reserve(const int& num);
    bool           contains(const std::string& name) const;
    void           extend(const GCOMHkds& hkds);
    void           load(const GFilename& filename);
    void           read(const GFitsTable& table);
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMHkds& hkds);
    void free_members(void);
    int  get_index(const std::string& name) const;

    // Protected data members
    std::vector<GCOMHkd> m_hkds; //!< Housekeeping Data containers
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMHkds").
 ***************************************************************************/
inline
std::string GCOMHkds::classname(void) const
{
    return ("GCOMHkds");
}


/***********************************************************************//**
 * @brief Return reference to Housekeeping Data container
 *
 * @param[in] index Housekeeping Data container index [0,...,size()-1].
 * @return Reference to Housekeeping Data container.
 *
 * Returns a reference to the Housekeeping Data container with the specified
 * @p index.
 ***************************************************************************/
inline
GCOMHkd& GCOMHkds::operator[](const int& index)
{
    return (m_hkds[index]);
}


/***********************************************************************//**
 * @brief Return reference to Housekeeping Data (const version)
 *
 * @param[in] index Housekeeping Data container index [0,...,size()-1].
 * @return Reference to Housekeeping Data container.
 *
 * Returns a reference to the Housekeeping Data container with the specified
 * @p index.
 ***************************************************************************/
inline
const GCOMHkd& GCOMHkds::operator[](const int& index) const
{
    return (m_hkds[index]);
}


/***********************************************************************//**
 * @brief Return number of Housekeeping parameters in collection
 *
 * @return Number of Housekeeping parameters in collection.
 *
 * Returns the number of Housekeeping parameters in the collection.
 ***************************************************************************/
inline
int GCOMHkds::size(void) const
{
    return (int)m_hkds.size();
}


/***********************************************************************//**
 * @brief Signals if there are no Housekeeping Data containers in collection
 *
 * @return True if collection is empty, false otherwise.
 *
 * Signals if the collection does not contain any Housekeeping Data
 * containers.
 ***************************************************************************/
inline
bool GCOMHkds::is_empty(void) const
{
    return (m_hkds.empty());
}


/***********************************************************************//**
 * @brief Reserves space for Housekeeping Data containers in collection
 *
 * @param[in] num Number of Housekeeping Data containers in collection.
 *
 * Reserves space for @p num Housekeeping Data containers in the collection.
 ***************************************************************************/
inline
void GCOMHkds::reserve(const int& num)
{
    m_hkds.reserve(num);
    return;
}

#endif /* GCOMHKDS_HPP */
