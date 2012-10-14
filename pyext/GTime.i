/***************************************************************************
 *                          GTime.i  -  Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @brief Time class python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTime.hpp"
#include "GTools.hpp"
%}
// Required for template vector
%include stl.i


/***********************************************************************//**
 * @class GTime
 *
 * @brief Class that handles times in a system independent way
 *
 * The GTime class stores a time value in seconds and its MJD reference in
 * days. The time can be retrieved in any MJD reference.
 ***************************************************************************/
class GTime {
public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    explicit GTime(const double&      time,
                   const double&      mrdref,
                   const std::string& timeunit,
                   const std::string& timesys = "TT",
                   const std::string& timeref = "local");
    explicit GTime(const double&      time,
                   const int&         mjdrefi,
                   const double&      mrdreff,
                   const std::string& timeunit,
                   const std::string& timesys = "TT",
                   const std::string& timeref = "local");
    virtual ~GTime(void);
 
    // Methods
    void   clear(void);
    GTime* clone(void) const;
    double jd(void) const;
    double mjd(void) const;
    double met(void) const;
    void   jd(const double& time);
    void   mjd(const double& time);
    void   met(const double& time);
    void   time(const double&      time,
                const double&      mrdref,
                const std::string& timeunit,
                const std::string& timesys = "TT",
                const std::string& timeref = "local");
    void   time(const double&      time,
                const int&         mjdrefi,
                const double&      mrdreff,
                const std::string& timeunit,
                const std::string& timesys = "TT",
                const std::string& timeref = "local");
    double time(void) const;
    double time(const double& mrdref) const;
    double time(const int& mjdrefi, const double& mrdreff) const;
    double mjdref(void) const;
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
