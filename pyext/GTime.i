/***************************************************************************
 *                           GTime.i - Time class                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2019 by Juergen Knoedlseder                         *
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
    GTime(const std::string& time, const GTimeReference& ref);
    explicit GTime(const std::string& time);
    virtual ~GTime(void);

    // Operators
    GTime& operator+=(const double& seconds);
    GTime& operator-=(const double& seconds);

    // Methods
    void           clear(void);
    GTime*         clone(void) const;
    std::string    classname(void) const;
    double         jd(void) const;
    double         jd(const std::string& timesys) const;
    double         mjd(void) const;
    double         mjd(const std::string& timesys) const;
    const double&  secs(void) const;
    double         secs(const std::string& timesys) const;
    double         days(void) const;
    double         days(const std::string& timesys) const;
    std::string    utc(const int& precision = 0) const;
    double         gmst(void) const;
    double         gast(void) const;
    double         lmst(const double& geolon) const;
    double         last(const double& geolon) const;
    double         convert(const GTimeReference& ref) const;
    void           jd(const double& time);
    void           jd(const double& time, const std::string& timesys);
    void           mjd(const double& time);
    void           mjd(const double& time, const std::string& timesys);
    void           secs(const double& seconds);
    void           secs(const double& seconds, const std::string& timesys);
    void           days(const double& days);
    void           days(const double& days, const std::string& timesys);
    void           utc(const std::string& time);
    void           set(const double& time, const GTimeReference& ref);
    void           set(const std::string& time, const GTimeReference& ref);
    void           set(const std::string& time);
    void           now(void);
    GTimeReference reference(void) const;
};


/***********************************************************************//**
 * @brief GTime class extension
 ***************************************************************************/
%extend GTime {
    GTime __add__(const double& seconds) const {
        return ((*self) + seconds);
    }
    GTime __radd__(const double& seconds) const {
        return (seconds + (*self));
    }
    GTime __sub__(const double& seconds) const {
        return ((*self) - seconds);
    }
    double __sub__(const GTime& time) const {
        return ((*self) - time);
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
    bool __le__(const GTime& time) const {
        return ((*self) <= time);
    }
    bool __ge__(const GTime& time) const {
        return ((*self) >= time);
    }
    GTime copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = self.secs(),
        return state
    def __setstate__(self, state):
        self.__init__()
        self.secs(state[0])
}

%pythoncode %{
    def datetime(self, *args):
        """
        Convert the GTime data into a datetime.datetime object.

        Usage 1:
        import gammalib
        t = gammalib.GTime()
        d = t.datetime()   # returns a datetime.datetime object

        Usage 2:
        import gammalib, datetime
        d = datetime.datetime.now()
        t = gammalib.GTime()
        t.datetime(d)      # set the GTime to the datetime's time.

        Parameters
        ----------
        args[0] : datetime.datetime object, if present, sets gtime to 
                   this time.

        Returns
        -------
        datetime : datetime.datetime object if no input argument,
                   otherwise returns nothing
        """
        # Import modules
        import time, datetime, calendar

        # If no datetime argument then return a datetime object
        if len(args) == 0:
            #f = '%Y-%m-%dT%H:%M:%S.%f %Z'
            #d = datetime.datetime.strptime(self.utc(6) + ' UTC', f)
            (t, msecs) = self.utc(6).strip().split('.')
            d = datetime.datetime.strptime(t + ' UTC', '%Y-%m-%dT%H:%M:%S %Z')
            d = d.replace(microsecond=int(msecs))
            return d

        # ... otherwise, if an argument is given, set the gtime to the datetime
        # argument
        elif len(args) == 1:
            dt = args[0]
            if type(dt) is datetime.datetime:
                #s = dt.strftime('%Y-%m-%dT%H:%M:%S.%f')
                s = dt.strftime('%Y-%m-%dT%H:%M:%S')+'.%6.6d' % (dt.microsecond)
                self.utc(s)
            else:
                msg  = 'Argument must be a datetime.datetime object, is '
                msg += 'currently ' + str(dt.__class__)
                raise TypeError(msg)

        # ... otherwise raise an exception
        else:
            msg  = 'GTime.datetime() needs 0 or 1 arguments.  It was given '
            msg += '%d.' % len(args)
            raise ValueError(msg)

        # Return
        return
%}
};
