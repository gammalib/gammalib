/***************************************************************************
 *            GPars.i - Application parameters Python interface            *
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
 * @file GPars.i
 * @brief Application parameter container class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPars.hpp"
#include "GTools.hpp"
#include "GException.hpp"
%}


/***********************************************************************//**
 * @class GPars
 *
 * @brief Application parameter container class
 ***************************************************************************/
class GPars : public GBase {

public:
    // Constructors and destructors
    GPars(void);
    explicit GPars(const std::string& filename);
    explicit GPars(const std::string& filename, const std::vector<std::string>& args);
    GPars(const GPars& pars);
    virtual ~GPars(void);
 
    // Methods
    void   clear(void);
    GPars* clone(void) const;
    int    size(void) const;
    void   append(const GPar& par);
    void   append_standard(void);
    void   load(const std::string& filename);
    void   load(const std::string& filename, const std::vector<std::string>& args);
    void   save(const std::string& filename);
    bool   haspar(const std::string& name) const;
};


/***********************************************************************//**
 * @brief GPars class extension
 ***************************************************************************/
%extend GPars {
    char *__str__() {
        return tochar(self->print());
    }
    GPar& __getitem__(const int& index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    GPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GPar& par) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = par;
            return;
        }
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    void __setitem__(const std::string& name, const GPar& par) {
        (*self)[name] = par;
        return;
    }
    GPars copy() {
        return (*self);
    }
}
