/***************************************************************************
 *                GCaldb.hpp  -  Calibration database class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Juergen Knoedlseder                              *
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
 * @file GCaldb.hpp
 * @brief Calibration database class interface definition
 * @author J. Knoedlseder
 */

#ifndef GCALDB_HPP
#define GCALDB_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <string>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GCaldb
 *
 * @brief Interface for the Calibration database class
 *
 * This class holds the definition of a calibration database. If the void
 * constructor is invoked it tries setting the calibration database access
 * path using the CALDB environment variable. Alternatively, the calibration
 * database path can be passed to the constructor to override the environment
 * variable.
 ***************************************************************************/
class GCaldb {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GCaldb& caldb);
    friend GLog&         operator<<(GLog& log,        const GCaldb& caldb);

public:
    // Constructors and destructors
    GCaldb(void);
    GCaldb(const GCaldb& caldb);
    explicit GCaldb(const std::string& pathname);
    virtual ~GCaldb(void);

    // Operators
    GCaldb& operator= (const GCaldb& caldb);

    // Methods
    void        clear(void);
    GCaldb*     clone(void) const;
    std::string dir(void) const { return m_caldb; }
    void        dir(const std::string& pathname);
    std::string print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCaldb& caldb);
    void free_members(void);
    void set_database(const std::string& pathname);

    // Protected data area
    std::string m_caldb;    //!< CALDB root directory
};

#endif /* GCALDB_HPP */
