/***************************************************************************
 *               GPulsarEphemeris.i - Pulsar ephemeris class               *
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
 * @file GPulsarEphemeris.i
 * @brief Pulsar ephemeris class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPulsarEphemeris.hpp"
%}


/***********************************************************************//**
 * @class GPulsarEphemeris
 *
 * @brief Pulsar ephemeris class
 ***************************************************************************/
class GPulsarEphemeris : public GBase {

public:
    // Constructors and destructors
    GPulsarEphemeris(void);
    GPulsarEphemeris(const GPulsarEphemeris& ephemeris);
    virtual ~GPulsarEphemeris(void);

    // Implemented pure virtual base class methods
    virtual void              clear(void);
    virtual GPulsarEphemeris* clone(void) const;
    virtual std::string       classname(void) const;

    // Other methods
    const std::string& name(void) const;
    void               name(const std::string& name);
    const GSkyDir&     dir(void) const;
    void               dir(const GSkyDir& dir);
    const GTime&       tstart(void) const;
    void               tstart(const GTime& tstart);
    const GTime&       tstop(void) const;
    void               tstop(const GTime& tstop);
    GTime              t0(void) const;
    void               t0(const GTime& t0);
    double             f0(void) const;
    void               f0(const double& f0);
    double             f1(void) const;
    void               f1(const double& f1);
    double             f2(void) const;
    void               f2(const double& f2);
    double             phase(const GTime& time) const;
};


/***********************************************************************//**
 * @brief GPulsarEphemeris class extension
 ***************************************************************************/
%extend GPulsarEphemeris {
    GPulsarEphemeris copy() {
        return (*self);
    }
};
