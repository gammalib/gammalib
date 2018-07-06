/***************************************************************************
 *              GXXXEventList.i - [INSTRUMENT] event list class            *
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
 * @file GXXXEventList.i
 * @brief [INSTRUMENT] event list class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXEventList.hpp"
%}


/***********************************************************************//**
 * @class GXXXEventList
 *
 * @brief [INSTRUMENT] event list class
 ***************************************************************************/
class GXXXEventList : public GEventList {

public:
    // Constructors and destructors
    GXXXEventList(void);
    explicit GXXXEventList(const GFilename& filename);
    GXXXEventList(const GXXXEventList& list);
    virtual ~GXXXEventList(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GXXXEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GXXXRoi& roi(void) const;

    // Other methods
    void append(const GXXXEventAtom& event);
    void reserve(const int& number);
    void remove(const int& index, const int& number = 1);
};


/***********************************************************************//**
 * @brief GXXXEventList class extension
 ***************************************************************************/
%extend GXXXEventList {
    GXXXEventAtom* __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Event index",
                                           index, self->size());
        }
    }
    GXXXEventList* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GXXXEventList* list = new GXXXEventList;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        list->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        list->append(*(*self)[i]);
                    }
                }
                return list;
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
    void __setitem__(const int& index, const GXXXEventAtom& event) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            *(*self)[index] = event;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            *(*self)[self->size()+index] = event;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "Event index",
                                           index, self->size());
        }
    }
    GXXXEventList copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = ()
        return state
    def __setstate__(self, state):
        self.__init__()
}
};
