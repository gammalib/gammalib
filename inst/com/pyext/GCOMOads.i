/***************************************************************************
 *          GCOMOads.i - COMPTEL Orbit Aspect Data container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2018 by Juergen Knodlseder                          *
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
 * @file GCOMOads.i
 * @brief COMPTEL Orbit Aspect Data container class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMOads.hpp"
%}


/***********************************************************************//**
 * @class GCOMOads
 *
 * @brief COMPTEL Orbit Aspect Data container class
 ***************************************************************************/
class GCOMOads : public GContainer {

public:
    // Constructors and destructors
    GCOMOads(void);
    explicit GCOMOads(const GFilename& filename);
    GCOMOads(const GCOMOads& oads);
    virtual ~GCOMOads(void);

    // Methods
    void        clear(void);
    GCOMOads*   clone(void) const;
    std::string classname(void) const;
    int         size(void) const;
    bool        is_empty(void) const;
    GCOMOad&    append(const GCOMOad& oad);
    GCOMOad&    insert(const int& index, const GCOMOad& oad);
    void        remove(const int& index);
    void        reserve(const int& num);
    void        extend(const GCOMOads& oads);
    void        load(const GFilename& filename);
    void        read(const GFitsTable& table);
};


/***********************************************************************//**
 * @brief GCOMOads class extension
 ***************************************************************************/
%extend GCOMOads {
    GCOMOad& __getitem__(const int& index) {
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
    GCOMOads* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GCOMOads* oads = new GCOMOads;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        oads->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        oads->append((*self)[i]);
                    }
                }
                return oads;
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
    void __setitem__(const int& index, const GCOMOad& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GCOMOads copy() {
        return (*self);
    }
};
