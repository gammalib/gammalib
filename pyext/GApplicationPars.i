/***************************************************************************
 *           GApplicationPars.i - Application parameter container          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
#include "GFilename.hpp"
%}


/***********************************************************************//**
 * @class GApplicationPars
 *
 * @brief Application parameter container class
 ***************************************************************************/
class GApplicationPars : public GContainer {

public:
    // Constructors and destructors
    GApplicationPars(void);
    explicit GApplicationPars(const GFilename& filename);
    GApplicationPars(const GFilename& filename,
                     const std::vector<std::string>& args);
    GApplicationPars(const GApplicationPars& pars);
    virtual ~GApplicationPars(void);
 
    // Methods
    void              clear(void);
    GApplicationPars* clone(void) const;
    std::string       classname(void) const;
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
    void              load(const GFilename& filename);
    void              load(const GFilename& filename,
                           const std::vector<std::string>& args);
    void              save(const GFilename& filename);
};


/***********************************************************************//**
 * @brief GApplicationPars class extension
 ***************************************************************************/
%extend GApplicationPars {
    GApplicationPar& __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    GApplicationPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    GApplicationPars* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PySlice_GetIndices((PySliceObject*)param, len, &start, &stop, &step) == 0) {
                GApplicationPars* pars = new GApplicationPars;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        pars->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        pars->append((*self)[i]);
                    }
                }
                return pars;
            }
            else {
                throw GException::invalid_argument("__getitem__(PyObject)",
                                                   "Invalid slice indices");
            }
        }
        else {
            throw GException::invalid_argument("__getitem__(PyObject)","");
        }
    }
    void __setitem__(const int& index, const GApplicationPar& par) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = par;
            return;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = par;
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
