/***************************************************************************
 *           GCOMD1Response.i - COMPTEL D1 module response class           *
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
 * @file GCOMD1Response.i
 * @brief COMPTEL D1 module response class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMD1Response.hpp"
%}


/***********************************************************************//**
 * @class GCOMD1Response
 *
 * @brief Interface for the COMPTEL D1 module response class
 ***************************************************************************/
class GCOMD1Response : public GBase {

public:
    // Constructors and destructors
    GCOMD1Response(void);
    GCOMD1Response(const GCOMD1Response& rsp);
    GCOMD1Response(const GCaldb& caldb, const std::string& sdaname);
    ~GCOMD1Response(void);

    // Operators
    double operator()(const GEnergy& etrue, const GEnergy& ereco) const;

    // Methods
    void            clear(void);
    GCOMD1Response* clone(void) const;
    std::string     classname(void) const;
    void            caldb(const GCaldb& caldb);
    const GCaldb&   caldb(void) const;
    void            load(const std::string& sdaname);
    void            read(const GFitsTable& hdu);
};


/***********************************************************************//**
 * @brief GCOMD1Response class extension
 ***************************************************************************/
%extend GCOMD1Response {
    GCOMD1Response copy() {
        return (*self);
    }
};
