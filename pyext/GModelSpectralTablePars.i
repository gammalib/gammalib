/***************************************************************************
 *  GModelSpectralTablePars.i - Spectral table model par container class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2021 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTablePars.i
 * @brief Spectral table model parameter container class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralTablePars.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralTablePars
 *
 * @brief Spectral table model parameter container class
 ***************************************************************************/
class GModelSpectralTablePars : public GContainer {

public:
    // Constructors and destructors
    GModelSpectralTablePars(void);
    GModelSpectralTablePars(const GModelSpectralTablePars& pars);
    virtual ~GModelSpectralTablePars(void);

    // Implemented pure virtual base class methods
    void                     clear(void);
    GModelSpectralTablePars* clone(void) const;
    std::string              classname(void) const;
    int                      size(void) const;
    bool                     is_empty(void) const;
    GModelSpectralTablePar*  set(const int&                    index,
                                 const GModelSpectralTablePar& par);
    GModelSpectralTablePar*  set(const std::string&            name,
                                 const GModelSpectralTablePar& par);
    GModelSpectralTablePar*  append(const GModelSpectralTablePar& par);
    GModelSpectralTablePar*  insert(const int&                    index,
                                    const GModelSpectralTablePar& par);
    GModelSpectralTablePar*  insert(const std::string&            name,
                                    const GModelSpectralTablePar& par);
    void                     remove(const int& index);
    void                     remove(const std::string& name);
    void                     reserve(const int& num);
    void                     extend(const GModelSpectralTablePars& pars);
    bool                     contains(const std::string& name) const;
};


/***********************************************************************//**
 * @brief GModelSpectralTablePars class extension
 ***************************************************************************/
%extend GModelSpectralTablePars {
    GModelSpectralTablePar* __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Parameter index",
                                           index, self->size());
        }
    }
    GModelSpectralTablePar* __getitem__(const std::string& name) {
        return (*self)[name];
    }
    GModelSpectralTablePars* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GModelSpectralTablePars* pars = new GModelSpectralTablePars;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        pars->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        pars->append(*(*self)[i]);
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
    void __setitem__(const int& index, const GModelSpectralTablePar& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            self->set(index, val);
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            self->set(self->size()+index, val);
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Parameter index",
                                           index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GModelSpectralTablePar& val) {
        self->set(name, val);
        return;
    }
    GModelSpectralTablePars copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = tuple([x for x in self]),
        return state
    def __setstate__(self, state):
        self.__init__()
        size = len(state[0])
        self.reserve(size)
        for x in state[0]:
            self.append(x)
}
};
