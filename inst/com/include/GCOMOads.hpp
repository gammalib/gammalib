/***************************************************************************
 *         GCOMOads.hpp - COMPTEL Orbit Aspect Data container class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knodlseder                          *
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
 * @file GCOMOads.hpp
 * @brief COMPTEL Orbit Aspect Data container class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMOADS_HPP
#define GCOMOADS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCOMOad.hpp"
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFitsTable;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMOads
 *
 * @brief COMPTEL Orbit Aspect Data container class
 *
 * The COMPTEL Orbit Aspect Data container class holds records of Orbit
 * Aspect data that were extracted from one COMPTEL OAD FITS file.
 ***************************************************************************/
class GCOMOads : public GContainer {

public:
    // Constructors and destructors
    GCOMOads(void);
    explicit GCOMOads(const GFilename& filename);
    GCOMOads(const GCOMOads& oads);
    virtual ~GCOMOads(void);

    // Operators
    GCOMOads&      operator=(const GCOMOads& oads);
    GCOMOad&       operator[](const int& index);
    const GCOMOad& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GCOMOads*      clone(void) const;
    std::string    classname(void) const;
    GCOMOad&       at(const int& index);
    const GCOMOad& at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GCOMOad&       append(const GCOMOad& oad);
    GCOMOad&       insert(const int& index, const GCOMOad& oad);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GCOMOads& oads);
    void           load(const GFilename& filename);
    void           read(const GFitsTable& table);
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMOads& oads);
    void free_members(void);

    // Protected data members
    std::vector<GCOMOad> m_oads; //!< Orbit Aspect Data records
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMOads").
 ***************************************************************************/
inline
std::string GCOMOads::classname(void) const
{
    return ("GCOMOads");
}


/***********************************************************************//**
 * @brief Return reference to Orbit Aspect Data
 *
 * @param[in] index Orbit Aspect Data index [0,...,size()-1].
 *
 * Returns a reference to the Orbit Aspect Data with the specified @p index.
 ***************************************************************************/
inline
GCOMOad& GCOMOads::operator[](const int& index)
{
    return (m_oads[index]);
}


/***********************************************************************//**
 * @brief Return reference to Orbit Aspect Data (const version)
 *
 * @param[in] index Orbit Aspect Data index [0,...,size()-1].
 *
 * Returns a reference to the Orbit Aspect Data with the specified @p index.
 ***************************************************************************/
inline
const GCOMOad& GCOMOads::operator[](const int& index) const
{
    return (m_oads[index]);
}


/***********************************************************************//**
 * @brief Return number of Orbit Aspect Data in container
 *
 * @return Number of Orbit Aspect Data in container.
 *
 * Returns the number of Orbit Aspect Data in the container.
 ***************************************************************************/
inline
int GCOMOads::size(void) const
{
    return (int)m_oads.size();
}


/***********************************************************************//**
 * @brief Signals if there are no Orbit Aspect Data in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the Orbit Aspect Data container does not contain any Orbit
 * Aspect Data.
 ***************************************************************************/
inline
bool GCOMOads::is_empty(void) const
{
    return (m_oads.empty());
}


/***********************************************************************//**
 * @brief Reserves space for Orbit Aspect Data in container
 *
 * @param[in] num Number of Orbit Aspect Data.
 *
 * Reserves space for @p num Orbit Aspect Data in the container.
 ***************************************************************************/
inline
void GCOMOads::reserve(const int& num)
{
    m_oads.reserve(num);
    return;
}

#endif /* GCOMOADS_HPP */
