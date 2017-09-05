/***************************************************************************
 *               GCOMOad.i - COMPTEL Orbit Aspect Data class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMOad.i
 * @brief COMPTEL Orbit Aspect Data class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMOad.hpp"
%}


/***********************************************************************//**
 * @class GCOMOad
 *
 * @brief COMPTEL Orbit Aspect Data class
 ***************************************************************************/
class GCOMOad : public GBase {

public:
    // Constructors and destructors
    GCOMOad(void);
    GCOMOad(const GCOMOad& oad);
    virtual ~GCOMOad(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMOad*    clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    const GTime& tstart(void) const;
    void         tstart(const GTime& tstart);
    const GTime& tstop(void) const;
    void         tstop(const GTime& tstop);
    const int&   tjd(void) const;
    void         tjd(const int& tjd);
    const int&   tics(void) const;
    void         tics(const int& tics);
};


/***********************************************************************//**
 * @brief GCOMOad class extension
 ***************************************************************************/
%extend GCOMOad {
    GCOMOad copy() {
        return (*self);
    }
};
