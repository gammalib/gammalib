/***************************************************************************
 *          GCOMHkd.i - COMPTEL Housekeeping Data container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMHkd.i
 * @brief COMPTEL Housekeeping Data container class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMHkd.hpp"
%}


/***********************************************************************//**
 * @class GCOMHkd
 *
 * @brief COMPTEL Housekeeping Data container class
 ***************************************************************************/
class GCOMHkd : public GContainer {

public:
    // Constructors and destructors
    GCOMHkd(void);
    explicit GCOMHkd(const std::string& name);
    GCOMHkd(const GCOMHkd& hkd);
    virtual ~GCOMHkd(void);

    // Methods
    void               clear(void);
    GCOMHkd*           clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    bool               is_empty(void) const;
    void               append(const GTime& time, const double& value);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               extend(const GCOMHkd& hkd);
    const std::string& name(void) const;
    void               name(const std::string& name);
    const GTime&       time(const int& index) const;
    void               time(const int& index, const GTime& time);
    const double&      value(const int& index) const;
    void               value(const int& index, const double& value);
};


/***********************************************************************//**
 * @brief GCOMHkd class extension
 ***************************************************************************/
%extend GCOMHkd {
    GCOMHkd copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = {'name'  :  self.name(),
                 'times' : [self.time(i)  for i in range(self.size())],
                 'values': [self.value(i) for i in range(self.size())]}
        return state
    def __setstate__(self, state):
        self.__init__(state['name'])
        size = len(state['times'])
        self.reserve(size)
        for i in range(size):
            self.append(state['times'][i], state['values'][i])
}
};
