/***************************************************************************
 *                  GPhases.hpp - Phase intervals class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GPhases.hpp
 * @brief Phase intervals class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GPHASES_HPP
#define GPHASES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GPhases
 *
 * @brief Phase Intervals class
 ***************************************************************************/
class GPhases : public GContainer {

public:
    // Constructors and destructors
    GPhases(void);
    GPhases(const GPhases& phases);
    GPhases(const double& pmin, const double& pmax);
    virtual ~GPhases(void);

    // Operators
    GPhases& operator=(const GPhases& phases);

    // Methods
    void        clear(void);
    GPhases*    clone(void) const;
    std::string classname(void) const;
    int         size(void) const;
    bool        is_empty(void) const;
    bool        contains(const double& phase) const;
    void        append(const double& pmin, const double& pmax);
    void        remove(const int& index);
    void        reserve(const int& num);
    void        extend(const GPhases& phases);
    double      pmin(const int& index) const;
    double      pmax(const int& index) const;
    std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GPhases& phases);
    void  free_members(void);
    void  insert_interval(const int& index, const double& pmin,
                                            const double& pmax);

    // Protected data area
    std::vector<double> m_pmin; // Lower phase interval boundaries
    std::vector<double> m_pmax; // Upper phase interval boundaries
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPhases").
 ***************************************************************************/
inline
std::string GPhases::classname(void) const
{
    return ("GPhases");
}


/***********************************************************************//**
 * @brief Return number of phase intervals
 *
 * @return Number of phase intervals.
 ***************************************************************************/
inline
int GPhases::size(void) const
{
    return (int)m_pmin.size();
}


/***********************************************************************//**
 * @brief Signal if there are no phase intervals
 *
 * @return True if there are no phase intervals.
 ***************************************************************************/
inline
bool GPhases::is_empty(void) const
{
    return (m_pmin.size() == 0);
}

#endif /* GPHASES_HPP */
