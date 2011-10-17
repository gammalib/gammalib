/***************************************************************************
 *                     GTime.i  -  Time class python I/F                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GTime.i
 * @brief GTime class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTime.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GTime
 *
 * @brief Class that handles times in a system independent way.
 *
 * The GTime class stores a time value in MJD and implements methods that
 * provide automatic conversion of the time values in other systems. This
 * makes instrument specific implementations more robust and reduces the
 * risk of unit errors.
 ***************************************************************************/
class GTime {
public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    ~GTime(void);
 
    // Methods
    void   clear(void) { m_time = 0.0; }
    double jd(void) const;
    double mjd(void) const;
    double met(void) const;
    void   jd(const double& time);
    void   mjd(const double& time);
    void   met(const double& time);
};


/***********************************************************************//**
 * @brief GTime class extension
 ***************************************************************************/
%extend GTime {
    char *__str__() {
        return tochar(self->print());
    }
    GTime __add__(const GTime& time) const {
        return ((*self) + time);
    }
    GTime __sub__(const GTime& time) const {
        return ((*self) - time);
    }
    GTime __mul__(const double& factor) const {
        return ((*self) * factor);
    }
    GTime __div__(const double& factor) const {
        return ((*self) / factor);
    }
    bool __eq__(const GTime& time) const {
        return ((*self) == time);
    }
    bool __ne__(const GTime& time) const {
        return ((*self) != time);
    }
    bool __lt__(const GTime& time) const {
        return ((*self) < time);
    }
    bool __gt__(const GTime& time) const {
        return ((*self) > time);
    }
    bool __lte__(const GTime& time) const {
        return ((*self) <= time);
    }
    bool __gte__(const GTime& time) const {
        return ((*self) >= time);
    }
    GTime copy() {
        return (*self);
    }
};


/***************************************************************************
 *                                 Typedefs                                *
 ***************************************************************************/
typedef std::vector<GTime> GTimes;
%template(GTimes) std::vector<GTime>;
