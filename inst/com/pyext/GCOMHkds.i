/***************************************************************************
 *         GCOMHkds.i - COMPTEL Housekeeping Data collection class         *
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
 * @file GCOMHkds.i
 * @brief COMPTEL Housekeeping Data collection class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMHkds.hpp"
%}


/***********************************************************************//**
 * @class GCOMHkds
 *
 * @brief COMPTEL Housekeeping Data collection class
 ***************************************************************************/
class GCOMHkds : public GContainer {

public:
    // Constructors and destructors
    GCOMHkds(void);
    explicit GCOMHkds(const GFilename& filename);
    GCOMHkds(const GCOMHkds& hkds);
    virtual ~GCOMHkds(void);

    // Methods
    void           clear(void);
    GCOMHkds*      clone(void) const;
    std::string    classname(void) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GCOMHkd&       set(const int& index, const GCOMHkd& hkd);
    GCOMHkd&       set(const std::string& name, const GCOMHkd& hkd);
    GCOMHkd&       append(const GCOMHkd& hkd);
    GCOMHkd&       insert(const int& index, const GCOMHkd& hkd);
    void           remove(const int& index);
    void           reserve(const int& num);
    bool           contains(const std::string& name) const;
    void           extend(const GCOMHkds& hkds);
    void           load(const GFilename& filename);
    void           read(const GFitsTable& table);
};


/***********************************************************************//**
 * @brief GCOMHkds class extension
 ***************************************************************************/
%extend GCOMHkds {
    GCOMHkd& __getitem__(const int& index) {
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
                                           "Housekeeping Data container index",
                                           index, self->size());
        }
    }
    GCOMHkd& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    GCOMHkds* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GCOMHkds* hkds = new GCOMHkds;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        hkds->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        hkds->append((*self)[i]);
                    }
                }
                return hkds;
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
    void __setitem__(const int& index, const GCOMHkd& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Housekeeping Data container index",
                                           index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GCOMHkd& val) {
        self->set(name, val);
        return;
    }
    GCOMHkds copy() {
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
