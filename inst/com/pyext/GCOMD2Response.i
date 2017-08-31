/***************************************************************************
 *           GCOMD2Response.i - COMPTEL D2 module response class           *
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
 * @file GCOMD2Response.i
 * @brief COMPTEL D2 module response class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMD2Response.hpp"
%}


/***********************************************************************//**
 * @class GCOMD2Response
 *
 * @brief Interface for the COMPTEL D2 module response class
 ***************************************************************************/
class GCOMD2Response : public GBase {

public:
    // Constructors and destructors
    GCOMD2Response(void);
    GCOMD2Response(const GCOMD2Response& rsp);
    GCOMD2Response(const GCaldb& caldb, const std::string& sdbname);
    ~GCOMD2Response(void);

    // Operators
    double operator()(const double& etrue, const double& ereco) const;

    // Methods
    void            clear(void);
    GCOMD2Response* clone(void) const;
    std::string     classname(void) const;
    void            caldb(const GCaldb& caldb);
    const GCaldb&   caldb(void) const;
    void            load(const std::string& sdaname);
    void            read(const GFitsTable& hdu);
    double          position(const double& etrue) const;
    double          sigma(const double& etrue) const;
    double          amplitude(const double& etrue) const;
    double          escape1(const double& etrue) const;
    double          escape2(const double& etrue) const;
    double          comptontail(const double& etrue) const;
    double          background(const double& etrue) const;
    double          emin(const double& etrue) const;
    double          ewidth(const double& etrue) const;
    double          emax(const double& etrue) const;
    double          emin(void) const;
    double          emax(void) const;
};


/***********************************************************************//**
 * @brief GCOMD2Response class extension
 ***************************************************************************/
%extend GCOMD2Response {
    GCOMD2Response copy() {
        return (*self);
    }
};
