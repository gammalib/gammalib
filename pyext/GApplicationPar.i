/***************************************************************************
 *                GApplicationPar.i - Application parameter class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GApplicationPar.i
 * @brief Application parameter class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GApplicationPar.hpp"
%}


/***********************************************************************//**
 * @class GApplicationPar
 *
 * @brief Application parameter class
 ***************************************************************************/
class GApplicationPar : public GBase {
public:
    // Constructors and destructors
    GApplicationPar(void);
    GApplicationPar(const std::string& name, const std::string& type,
                    const std::string& mode, const std::string& value,
                    const std::string& min, const std::string& max, 
                    const std::string& prompt);
    GApplicationPar(const GApplicationPar& par);
    virtual ~GApplicationPar(void);
 
    // Methods
    void                     clear(void);
    GApplicationPar*         clone(void) const;
    std::string              classname(void) const;
    void                     type(const std::string& type);
    void                     mode(const std::string& mode);
    void                     value(const std::string& value);
    void                     string(const std::string& value);
    void                     filename(const GFilename& value);
    void                     time(const GTime& value);
    void                     boolean(const bool& value);
    void                     integer(const int& value);
    void                     real(const double& value);
    const std::string&       name(void) const;
    const std::string&       type(void) const;
    const std::string&       mode(void) const;
    void                     query(void);
    std::string              value(void);
    std::string              string(void);
    GFilename                filename(void);
    GTime                    time(void);
    GTime                    time(const GTimeReference& ref);
    bool                     boolean(void);
    int                      integer(void);
    double                   real(void);
    const std::string&       current_value(void) const;
    const std::string&       min(void) const;
    const std::string&       max(void) const;
    const std::string&       prompt(void) const;
    bool                     is_learn(void) const;
    bool                     is_query(void) const;
    bool                     is_filename(void) const;
    bool                     is_valid(void);
    bool                     is_undefined(void);
    bool                     is_notanumber(void);
    bool                     was_queried(void) const;
    void                     pickle(const std::vector<std::string>& string);
    std::vector<std::string> pickle(void) const;
};


/***********************************************************************//**
 * @brief GApplicationPar class extension
 ***************************************************************************/
%extend GApplicationPar {
    GApplicationPar copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.pickle(),)
        return state
    def __setstate__(self, state):
        self.__init__()
        self.pickle(state[0])
}
}
