/***************************************************************************
 *              GCOMStatus.i - COMPTEL instrument status class             *
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
 * @file GCOMStatus.i
 * @brief COMPTEL instrument status class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMStatus.hpp"
%}


/***********************************************************************//**
 * @class GCOMStatus
 *
 * @brief COMPTEL instrument status class
 ***************************************************************************/
class GCOMStatus : public GBase {

public:
    // Constructors and destructors
    GCOMStatus(void);
    GCOMStatus(const GCOMStatus& status);
    virtual ~GCOMStatus(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMStatus* clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    void load(void) const;
    int  d1status(const int& tjd, const int& module) const;
    int  d2status(const int& tjd, const int& module) const;
};


/***********************************************************************//**
 * @brief GCOMStatus class extension
 ***************************************************************************/
%extend GCOMStatus {
    GCOMStatus copy() {
        return (*self);
    }
};
