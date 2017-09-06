/***************************************************************************
 *             GCOMStatus.hpp - COMPTEL instrument status class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMStatus.hpp
 * @brief COMPTEL instrument status class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMSTATUS_HPP
#define GCOMSTATUS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMStatus
 *
 * @brief COMPTEL instrument status class
 *
 * The class implements a database of the COMPTEL instrument status. For the
 * time being the database contains the status of the D1 and D2 modules for
 * all days that are covered by the CGRO mission.
 ***************************************************************************/
class GCOMStatus : public GBase {

public:
    // Constructors and destructors
    GCOMStatus(void);
    GCOMStatus(const GCOMStatus& status);
    virtual ~GCOMStatus(void);

    // Operators
    GCOMStatus& operator=(const GCOMStatus& status);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMStatus* clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void load(void) const;
    int  d1status(const int& tjd, const int& module) const;
    int  d2status(const int& tjd, const int& module) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMStatus& status);
    void free_members(void);
    void update_cache(const int& tjd) const;
    void load_status(void) const;
    
    // Protected members
    mutable std::vector<int>               m_tjds;     //!< TJD for status
    mutable std::vector<std::vector<int> > m_d1status; //!< D1 module status
    mutable std::vector<std::vector<int> > m_d2status; //!< D2 module status

    // Cache for last module status
    mutable int              m_last_tjd;      //!< Last TJD
    mutable std::vector<int> m_last_d1status; //!< Last D1 module status
    mutable std::vector<int> m_last_d2status; //!< Last D2 module status
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMStatus").
 ***************************************************************************/
inline
std::string GCOMStatus::classname(void) const
{
    return ("GCOMStatus");
}

#endif /* GCOMSTATUS_HPP */
