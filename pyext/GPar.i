/***************************************************************************
 *                   GPar.i - Application parameter class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GPar.i
 * @brief Application parameter class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPar.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPar
 *
 * @brief Application parameter class
 ***************************************************************************/
class GPar : public GBase {
public:
    // Constructors and destructors
    GPar(void);
    explicit GPar(const std::string& name, const std::string& type,
                  const std::string& mode, const std::string& value,
                  const std::string& min, const std::string& max, 
                  const std::string& prompt);
    GPar(const GPar& par);
    virtual ~GPar(void);
 
    // Methods
    void               clear(void);
    GPar*              clone(void) const;
    void               type(const std::string& type);
    void               mode(const std::string& mode);
    void               value(const std::string& value);
    void               string(const std::string& value);
    void               filename(const std::string& value);
    void               boolean(const bool& value);
    void               integer(const int& value);
    void               real(const double& value);
    const std::string& name(void) const;
    const std::string& type(void) const;
    const std::string& mode(void) const;
    std::string        value(void);
    std::string        string(void);
    std::string        filename(void);
    bool               boolean(void);
    int                integer(void);
    double             real(void);
    const std::string& min(void) const;
    const std::string& max(void) const;
    const std::string& prompt(void) const;
    bool               is_learn(void) const;
    bool               is_query(void) const;
    bool               is_filename(void) const;
    bool               is_valid(void);
    bool               is_undefined(void);
    bool               is_notanumber(void);
};


/***********************************************************************//**
 * @brief GPar class extension
 ***************************************************************************/
%extend GPar {
    GPar copy() {
        return (*self);
    }
}
