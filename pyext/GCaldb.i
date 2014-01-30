/***************************************************************************
 *                  GCaldb.i - Calibration database class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GCaldb.i
 * @brief Calibration database class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCaldb.hpp"
%}


/***********************************************************************//**
 * @class GCaldb
 *
 * @brief Interface for the Calibration database class
 ***************************************************************************/
class GCaldb : public GBase {
public:
    // Constructors and destructors
    GCaldb(void);
    GCaldb(const GCaldb& caldb);
    explicit GCaldb(const std::string& pathname);
    GCaldb(const std::string& mission, const std::string& instrument);
    virtual ~GCaldb(void);

    // Methods
    void               clear(void);
    GCaldb*            clone(void) const;
    int                size(void) const;
    const std::string& dir(void) const;
    void               dir(const std::string& pathname);
    void               open(const std::string& mission,
                            const std::string& instrument = "");
    void               close(void);
    std::string        filename(const std::string& detector,
                                const std::string& filter,
                                const std::string& codename,
                                const std::string& date,
                                const std::string& time,
                                const std::string& expr);
};


/***********************************************************************//**
 * @brief GCaldb class extension
 ***************************************************************************/
%extend GCaldb {
    GCaldb copy() {
        return (*self);
    }
};
