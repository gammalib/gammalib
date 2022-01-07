/***************************************************************************
 *   GCOMBvcs.hpp - COMPTEL Solar System Barycentre Data container class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvcs.hpp
 * @brief COMPTEL Solar System Barycentre Data container class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMBVCS_HPP
#define GCOMBVCS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCOMBvc.hpp"
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFitsTable;
class GCOMOad;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMBvcs
 *
 * @brief COMPTEL Solar System Barycentre Data container class
 *
 * The COMPTEL Solar System Barycentre Data container class holds records
 * of Solar System Barycentre data that were extracted from one COMPTEL BVC
 * FITS file.
 ***************************************************************************/
class GCOMBvcs : public GContainer {

public:
    // Constructors and destructors
    GCOMBvcs(void);
    explicit GCOMBvcs(const GFilename& filename);
    GCOMBvcs(const GCOMBvcs& bvcs);
    virtual ~GCOMBvcs(void);

    // Operators
    GCOMBvcs&      operator=(const GCOMBvcs& bvcs);
    GCOMBvc&       operator[](const int& index);
    const GCOMBvc& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GCOMBvcs*      clone(void) const;
    std::string    classname(void) const;
    GCOMBvc&       at(const int& index);
    const GCOMBvc& at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GCOMBvc&       append(const GCOMBvc& bvc);
    GCOMBvc&       insert(const int& index, const GCOMBvc& bvc);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GCOMBvcs& bvcs);
    void           load(const GFilename& filename);
    void           read(const GFitsTable& table);
    const GCOMBvc* find(const GCOMOad& oad) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMBvcs& bvcs);
    void free_members(void);

    // Protected data members
    std::vector<GCOMBvc> m_bvcs; //!< Solar System Barycentre Data records
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMBvcs").
 ***************************************************************************/
inline
std::string GCOMBvcs::classname(void) const
{
    return ("GCOMBvcs");
}


/***********************************************************************//**
 * @brief Return reference to Solar System Barycentre Data
 *
 * @param[in] index Solar System Barycentre Data index [0,...,size()-1].
 *
 * Returns a reference to the Solar System Barycentre Data with the
 * specified @p index.
 ***************************************************************************/
inline
GCOMBvc& GCOMBvcs::operator[](const int& index)
{
    return (m_bvcs[index]);
}


/***********************************************************************//**
 * @brief Return reference to Solar System Barycentre Data (const version)
 *
 * @param[in] index Solar System Barycentre Data index [0,...,size()-1].
 *
 * Returns a reference to the Solar System Barycentre Data with the specified
 * @p index.
 ***************************************************************************/
inline
const GCOMBvc& GCOMBvcs::operator[](const int& index) const
{
    return (m_bvcs[index]);
}


/***********************************************************************//**
 * @brief Return number of Solar System Barycentre Data in container
 *
 * @return Number of Solar System Barycentre Data in container.
 *
 * Returns the number of Solar System Barycentre Data in the container.
 ***************************************************************************/
inline
int GCOMBvcs::size(void) const
{
    return (int)m_bvcs.size();
}


/***********************************************************************//**
 * @brief Signals if there are no Solar System Barycentre Data in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the Solar System Barycentre Data container does not contain
 * any Solar System Barycentre Data.
 ***************************************************************************/
inline
bool GCOMBvcs::is_empty(void) const
{
    return (m_bvcs.empty());
}


/***********************************************************************//**
 * @brief Reserves space for Solar System Barycentre Data in container
 *
 * @param[in] num Number of Solar System Barycentre Data.
 *
 * Reserves space for @p num Solar System Barycentre Data in the container.
 ***************************************************************************/
inline
void GCOMBvcs::reserve(const int& num)
{
    m_bvcs.reserve(num);
    return;
}

#endif /* GCOMBVCS_HPP */
