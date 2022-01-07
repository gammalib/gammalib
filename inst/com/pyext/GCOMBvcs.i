/***************************************************************************
 *    GCOMBvcs.i - COMPTEL Solar System Barycentre Data container class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvcs.i
 * @brief COMPTEL Solar System Barycentre Data container class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMBvcs.hpp"
%}


/***********************************************************************//**
 * @class GCOMBvcs
 *
 * @brief COMPTEL Solar System Barycentre Data container class
 ***************************************************************************/
class GCOMBvcs : public GContainer {

public:
    // Constructors and destructors
    GCOMBvcs(void);
    explicit GCOMBvcs(const GFilename& filename);
    GCOMBvcs(const GCOMBvcs& bvcs);
    virtual ~GCOMBvcs(void);

    // Methods
    void           clear(void);
    GCOMBvcs*      clone(void) const;
    std::string    classname(void) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GCOMBvc&       append(const GCOMBvc& bvc);
    GCOMBvc&       insert(const int& index, const GCOMBvc& bvc);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GCOMBvcs& bvcs);
    void           load(const GFilename& filename);
    void           read(const GFitsTable& table);
    const GCOMBvc* find(const GCOMOad& oad) const;
    double         tdelta(const GSkyDir& dir, const GTime& time) const;
};


/***********************************************************************//**
 * @brief GCOMBvcs class extension
 ***************************************************************************/
%extend GCOMBvcs {
    GCOMBvc& __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "BVC index",
                                           index, self->size());
        }
    }
    GCOMBvcs* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GCOMBvcs* bvcs = new GCOMBvcs;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        bvcs->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        bvcs->append((*self)[i]);
                    }
                }
                return bvcs;
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
    void __setitem__(const int& index, const GCOMBvc& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "BVC index",
                                           index, self->size());
        }
    }
    GCOMBvcs copy() {
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
