/***************************************************************************
 *                           GTime.i - Time class                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @brief Time class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTime.hpp"
%}


/***********************************************************************//**
 * @class GTime
 *
 * @brief Time class
 ***************************************************************************/
class GTime : public GBase {

public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    GTime(const double& time, const std::string& unit = "sec");
    GTime(const double& time, const GTimeReference& ref);
    virtual ~GTime(void);

    // Methods
    void           clear(void);
    GTime*         clone(void) const;
    std::string    classname(void) const;
    double         jd(void) const;
    double         mjd(void) const;
    const double&  secs(void) const;
    double         days(void) const;
    std::string    utc(void) const;
    double         convert(const GTimeReference& ref) const;
    void           jd(const double& time);
    void           mjd(const double& time);
    void           secs(const double& seconds);
    void           days(const double& days);
    void           utc(const std::string& time);
    void           set(const double& time, const GTimeReference& ref);
    GTimeReference reference(void) const;
};


/***********************************************************************//**
 * @brief GTime class extension
 ***************************************************************************/
%extend GTime {
    GTime __add__(const GTime& time) const {
        return ((*self) + time);
    }
    GTime __sub__(const GTime& time) const {
        return ((*self) - time);
    }
    GTime __mul__(const double& factor) const {
        return ((*self) * factor);
    }
    // Python 2.x
    GTime __div__(const double& factor) const {
        return ((*self) / factor);
    }
    // Python 3.x
    GTime __truediv__(const double& factor) const {
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
