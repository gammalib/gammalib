/***************************************************************************
 *                    GEphemerides.i - Ephemerides class                   *
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
 * @file GEphemerides.i
 * @brief Ephemerides class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEphemerides.hpp"
%}


/***********************************************************************//**
 * @class GEphemerides
 *
 * @brief Ephemerides class
 ***************************************************************************/
class GEphemerides : public GBase {

public:
    // Constructors and destructors
    GEphemerides(void);
    GEphemerides(const GEphemerides& ephemerides);
    virtual ~GEphemerides(void);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GEphemerides* clone(void) const;
    virtual std::string   classname(void) const;

    // Other methods
    int                size(void) const;
    bool               is_empty(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    void               load(const GFilename& filename);
    void               ephemeris(const GTime& time,
                                 GVector*     rce,
                                 GVector*     rcs,
                                 GVector*     vce,
                                 double*      etut) const;
    double             geo2ssb(const GTime&   time,
                               const GSkyDir& srcdir) const;
    double             utc2tt(const GTime& time) const;
};


/***********************************************************************//**
 * @brief GEphemerides class extension
 ***************************************************************************/
%extend GEphemerides {
    GEphemerides copy() {
        return (*self);
    }
};
