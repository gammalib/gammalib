/***************************************************************************
 *                     GTimes.hpp - Time container class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GTimes.hpp
 * @brief Time container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GTIMES_HPP
#define GTIMES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GTimes
 *
 * @brief Container class for times.
 ***************************************************************************/
class GTimes : public GBase {

public:
    // Constructors and destructors
    GTimes(void);
    GTimes(const GTimes& times);
    virtual ~GTimes(void);
 
    // Operators
    GTimes&      operator=(const GTimes& times);
    GTime&       operator[](const int& index);
    const GTime& operator[](const int& index) const;

    // Methods
    void        clear(void);
    GTimes*     clone(void) const;
    int         size(void) const { return m_times.size(); }
    void        append(const GTime& time);
    void        reserve(const int& number);
    std::string print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTimes& times);
    void free_members(void);

    // Protected data members
    std::vector<GTime> m_times;      //!< Time list
};
#endif /* GTIMES_HPP */
