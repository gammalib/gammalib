/***************************************************************************
 *               GSkyDirs.hpp - Sky directions container class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSkyDirs.hpp
 * @brief Sky directions container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSKYDIRS_HPP
#define GSKYDIRS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GSkyDir.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSkyDirs
 *
 * @brief Sky directions container class
 *
 * This class is a container for sky directions. Sky directions are
 * implemented by the GSkyDir class.
 ***************************************************************************/
class GSkyDirs : public GContainer {

public:
    // Constructors and destructors
    GSkyDirs(void);
    explicit  GSkyDirs(const GSkyDir& dir);
    GSkyDirs(const GSkyDirs& dirs);
    virtual ~GSkyDirs(void);
 
    // Operators
    GSkyDirs&      operator=(const GSkyDirs& dirs);
    GSkyDir&       operator[](const int& index);
    const GSkyDir& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GSkyDirs*      clone(void) const;
    std::string    classname(void) const;
    GSkyDir&       at(const int& index);
    const GSkyDir& at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GSkyDir&       append(const GSkyDir& dir);
    GSkyDir&       insert(const int& index, const GSkyDir& dir);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GSkyDirs& dirs);
    std::string    print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSkyDirs& dirs);
    void free_members(void);

    // Protected data members
    std::vector<GSkyDir> m_dirs;  //!< List of sky directions
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyDirs").
 ***************************************************************************/
inline
std::string GSkyDirs::classname(void) const
{
    return ("GSkyDirs");
}


/***********************************************************************//**
 * @brief Return reference to sky direction
 *
 * @param[in] index Sky direction index [0,...,size()-1].
 *
 * Returns a reference to the sky direction with the specified @p index.
 ***************************************************************************/
inline
GSkyDir& GSkyDirs::operator[](const int& index)
{
    return (m_dirs[index]);
}


/***********************************************************************//**
 * @brief Return reference to sky direction (const version)
 *
 * @param[in] index Sky direction index [0,...,size()-1].
 *
 * Returns a const reference to the sky direction with the specified
 * @p index.
 ***************************************************************************/
inline
const GSkyDir& GSkyDirs::operator[](const int& index) const
{
    return (m_dirs[index]);
}


/***********************************************************************//**
 * @brief Return reference to sky direction
 *
 * @param[in] index Sky direction index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Sky direction index is out of range.
 *
 * Returns a reference to the sky direction with the specified @p index.
 ***************************************************************************/
inline
GSkyDir& GSkyDirs::at(const int& index)
{
    // Return reference to sky direction using const method
    return const_cast<GSkyDir&>(static_cast<const GSkyDirs*>(this)->at(index));
}


/***********************************************************************//**
 * @brief Return number of sky directions in container
 *
 * @return Number of sky directions in container.
 *
 * Returns the number of sky directions in the sky directions container.
 ***************************************************************************/
inline
int GSkyDirs::size(void) const
{
    return ((int)m_dirs.size());
}


/***********************************************************************//**
 * @brief Signals if there are no sky directions in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the sky directions container does not contain any sky
 * direction.
 ***************************************************************************/
inline
bool GSkyDirs::is_empty(void) const
{
    return (m_dirs.empty());
}


/***********************************************************************//**
 * @brief Reserves space for sky directions in container
 *
 * @param[in] num Number of sky directions.
 *
 * Reserves space for @p num sky directions in the container.
 ***************************************************************************/
inline
void GSkyDirs::reserve(const int& num)
{
    m_dirs.reserve(num);
    return;
}

#endif /* GSKYDIRS_HPP */
