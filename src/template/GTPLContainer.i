/***************************************************************************
 *                 GTPLContainer.i - [WHAT] container class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GTPLContainer.i
 * @brief [WHAT] container class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTPLContainer.hpp"
%}


/***********************************************************************//**
 * @class GTPLContainer
 *
 * @brief [WHAT] container class
 ***************************************************************************/
class GTPLContainer : public GContainer {

public:
    // Constructors and destructors
    GTPLContainer(void);
    GTPLContainer(const GTPLContainer& TPL_CONTAINER);
    virtual ~GTPLContainer(void);

    // Methods
    void           clear(void);
    GTPLContainer* clone(void) const;
    std::string    classname(void) const;
    GTPLBase&      at(const int& index);
    int            size(void) const;
    bool           is_empty(void) const;
    GTPLBase&      append(const GTPLBase& TPL_OBJECT);
    GTPLBase&      insert(const int& index, const GTPLBase& TPL_OBJECT);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GTPLContainer& TPL_CONTAINER);
};


/***********************************************************************//**
 * @brief GTPLContainer class extension
 ***************************************************************************/
%extend GTPLContainer {
    GTPLBase& __getitem__(const int& index) {
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
    GTPLContainer* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GTPLContainer* TPL_CONTAINER = new GTPLContainer;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        TPL_CONTAINER->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        TPL_CONTAINER->append((*self)[i]);
                    }
                }
                return TPL_CONTAINER;
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
    void __setitem__(const int& index, const GTPLBase& val) {
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
    GTPLContainer copy() {
        return (*self);
    }
};
