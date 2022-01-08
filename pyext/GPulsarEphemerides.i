/***************************************************************************
 *             GPulsarEphemerides.i - Pulsar ephemerides class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GPulsarEphemerides.i
 * @brief Pulsar ephemerides class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPulsarEphemerides.hpp"
%}


/***********************************************************************//**
 * @class GPulsarEphemerides
 *
 * @brief Pulsar ephemerides class
 ***************************************************************************/
class GPulsarEphemerides : public GBase {

public:
    // Constructors and destructors
    GPulsarEphemerides(void);
    GPulsarEphemerides(const GPulsarEphemerides& ephemerides);
    virtual ~GPulsarEphemerides(void);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GPulsarEphemerides* clone(void) const;
    virtual std::string         classname(void) const;

    // Other methods
    const GTime& tstart(void) const;
    const GTime& tstop(void) const;
    GTime        t0(void) const;
    void         t0(const GTime& t0);
    double       f0(void) const;
    void         f0(const double& f0);
    double       f1(void) const;
    void         f1(const double& f1);
    double       f2(void) const;
    void         f2(const double& f2);
    double       phase(const GTime& time) const;
};


/***********************************************************************//**
 * @brief GPulsarEphemerides class extension
 ***************************************************************************/
%extend GPulsarEphemerides {
    GPulsarEphemerides copy() {
        return (*self);
    }
};
