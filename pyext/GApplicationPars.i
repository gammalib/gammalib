/***************************************************************************
 *           GApplicationPars.i - Application parameter container          *
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
 * @file GApplicationPars.i
 * @brief Application parameter container class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GApplicationPars.hpp"
#include "GException.hpp"
%}


/***********************************************************************//**
 * @class GApplicationPars
 *
 * @brief Application parameter container class
 ***************************************************************************/
class GApplicationPars : public GBase {

public:
    // Constructors and destructors
    GApplicationPars(void);
    explicit GApplicationPars(const std::string& filename);
    GApplicationPars(const std::string& filename,
                     const std::vector<std::string>& args);
    GApplicationPars(const GApplicationPars& pars);
    virtual ~GApplicationPars(void);
 
    // Methods
    void              clear(void);
    GApplicationPars* clone(void) const;
    GApplicationPar&  at(const int& index);
    int               size(void) const;
    bool              is_empty(void) const;
    GApplicationPar&  append(const GApplicationPar& par);
    void              append_standard(void);
    GApplicationPar&  insert(const int& index, const GApplicationPar& par);
    GApplicationPar&  insert(const std::string& name,
                             const GApplicationPar& par);
    void              remove(const int& index);
    void              remove(const std::string& name);
    void              reserve(const int& num);
    void              extend(const GApplicationPars& pars);
    bool              contains(const std::string& name) const;
    void              load(const std::string& filename);
    void              load(const std::string& filename,
                           const std::vector<std::string>& args);
    void              save(const std::string& filename);
};


/***********************************************************************//**
 * @brief GApplicationPars class extension
 ***************************************************************************/
%extend GApplicationPars {
    GApplicationPar& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    GApplicationPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GApplicationPar& par) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = par;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GApplicationPar& par) {
        (*self)[name] = par;
        return;
    }
    GApplicationPars copy() {
        return (*self);
    }
}
