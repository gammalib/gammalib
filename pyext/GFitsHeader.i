/***************************************************************************
 *            GFitsHeader.i - FITS header cards container class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GFitsHeader.i
 * @brief FITS header cards container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHeader.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GFitsHeader
 *
 * @brief Interface for FITS header class
 ***************************************************************************/
class GFitsHeader : public GContainer {
public:
    // Constructors and destructors
    GFitsHeader(void);
    GFitsHeader(const GFitsHeader& header);
    virtual ~GFitsHeader(void);

    // Methods
    void             clear(void);
    GFitsHeader*     clone(void) const;
    std::string      classname(void) const;
    std::string      string(const int& cardno) const;
    std::string      string(const std::string& keyname) const;
    double           real(const int& cardno) const;
    double           real(const std::string& keyname) const;
    int              integer(const int& cardno) const;
    int              integer(const std::string& keyname) const;
    int              size(void) const;
    bool             is_empty(void) const;
    GFitsHeaderCard& append(const GFitsHeaderCard& card);
    GFitsHeaderCard& insert(const int& cardno, const GFitsHeaderCard& card);
    GFitsHeaderCard& insert(const std::string& keyname, const GFitsHeaderCard& card);
    void             remove(const int& cardno);
    void             remove(const std::string& keyname);
    void             reserve(const int& num);
    void             extend(const GFitsHeader& header);
    bool             contains(const int& cardno) const;
    bool             contains(const std::string& keyname) const;
    void             load(void* vptr);
    void             save(void* vptr) const;
};


/***********************************************************************//**
 * @brief GFitsHeader class extension
 ***************************************************************************/
%extend GFitsHeader {
    GFitsHeaderCard& __getitem__(const int& cardno) {
        // Counting from start, e.g. [2]
        if (cardno >= 0 && cardno < self->size()) {
            return (*self)[cardno];
        }
        // Counting from end, e.g. [-1]
        else if (cardno < 0 && self->size()+cardno >= 0) {
            return (*self)[self->size()+cardno];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Header card number",
                                           cardno, self->size());
        }
    }
    GFitsHeaderCard& __getitem__(const std::string& keyname) {
        return (*self)[keyname];
    }
    GFitsHeader* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GFitsHeader* header = new GFitsHeader;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        header->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        header->append((*self)[i]);
                    }
                }
                return header;
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
    void __setitem__(const int& cardno, const GFitsHeaderCard& card) {
        // Counting from start, e.g. [2]
        if (cardno >= 0 && cardno < self->size()) {
            (*self)[cardno] = card;
        }
        // Counting from end, e.g. [-1]
        else if (cardno < 0 && self->size()+cardno >= 0) {
            (*self)[self->size()+cardno] = card;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Header card number",
                                           cardno, self->size());
        }
    }
    void __setitem__(const std::string& keyname, const GFitsHeaderCard& card) {
        (*self)[keyname] = card;
        return;
    }
    GFitsHeader copy() {
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
}
