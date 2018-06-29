/***************************************************************************
 *                GCOMEventList.i - COMPTEL event list class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2018 by Juergen Knoedlseder                         *
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
 * @file GCOMEventList.i
 * @brief COMPTEL event list class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMEventList.hpp"
%}


/***********************************************************************//**
 * @class GCOMEventList
 *
 * @brief COMPTEL event list class
 ***************************************************************************/
class GCOMEventList : public GEventList {

public:
    // Constructors and destructors
    GCOMEventList(void);
    explicit GCOMEventList(const GFilename& filename);
    GCOMEventList(const GCOMEventList& list);
    virtual ~GCOMEventList(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GCOMRoi& roi(void) const;

    // Other methods
    void append(const GCOMEventAtom& event);
    void reserve(const int& number);
    void remove(const int& index, const int& number = 1);
};


/***********************************************************************//**
 * @brief GCOMEventList class extension
 ***************************************************************************/
%extend GCOMEventList {
    GCOMEventAtom* __getitem__(const int& index) {
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
    GCOMEventList* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GCOMEventList* list = new GCOMEventList;
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
    void __setitem__(const int& index, const GCOMEventAtom& event) {
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
    GCOMEventList copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        fits = gammalib.GFits()
        self.write(fits)
        state = (fits,)
        return state
    def __setstate__(self, state):
        self.__init__()
        if not state[0].is_empty():
            self.read(state[0])
}
};
