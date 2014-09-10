/***************************************************************************
 *                     GTimes.hpp - Time container class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
#include "GContainer.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GTimes
 *
 * @brief Time container class.
 *
 * This class is a container for times. Times are implemented by the GTime
 * class which stores time in a system independent way.
 ***************************************************************************/
class GTimes : public GContainer {

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
    std::string classname(void) const;
    int         size(void) const;
    bool        is_empty(void) const;
    void        append(const GTime& time);
    void        insert(const int& index, const GTime& time);
    void        remove(const int& index);
    void        reserve(const int& num);
    void        extend(const GTimes& times);
    std::string print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTimes& times);
    void free_members(void);

    // Protected data members
    std::vector<GTime> m_times;  //!< List of times
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GTimes").
 ***************************************************************************/
inline
std::string GTimes::classname(void) const
{
    return ("GTimes");
}


/***********************************************************************//**
 * @brief Return number of times
 *
 * @return Number of times.
 ***************************************************************************/
inline
int GTimes::size(void) const
{
    return m_times.size();
}


/***********************************************************************//**
 * @brief Signal if there are no times
 *
 * @return True if there are no times.
 ***************************************************************************/
inline
bool GTimes::is_empty(void) const
{
    return m_times.empty();
}

#endif /* GTIMES_HPP */
