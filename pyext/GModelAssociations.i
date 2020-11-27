/***************************************************************************
 *        GModelAssociations.i - Model associations container class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GModelAssociations.i
 * @brief Model associations container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelAssociations.hpp"
%}


/***********************************************************************//**
 * @class GModelAssociations
 *
 * @brief Model associations container class
 ***************************************************************************/
class GModelAssociations : public GContainer {

public:
    // Constructors and destructors
    GModelAssociations(void);
    explicit GModelAssociations(const GXmlElement& xml);
    GModelAssociations(const GModelAssociations& associations);
    virtual ~GModelAssociations(void);

    // Methods
    void                clear(void);
    GModelAssociations* clone(void) const;
    std::string         classname(void) const;
    int                 size(void) const;
    bool                is_empty(void) const;
    GModelAssociation&  append(const GModelAssociation& association);
    GModelAssociation&  insert(const int&               index,
                               const GModelAssociation& association);
    GModelAssociation&  insert(const std::string&       name,
                               const GModelAssociation& association);
    void                remove(const int& index);
    void                remove(const std::string& name);
    void                reserve(const int& num);
    void                extend(const GModelAssociations& associations);
    bool                contains(const std::string& name) const;
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelAssociations class extension
 ***************************************************************************/
%extend GModelAssociations {
    GModelAssociation& __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Association index",
                                           index, self->size());
        }
    }
    GModelAssociation& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    GModelAssociations* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GModelAssociations* associations = new GModelAssociations;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        associations->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        associations->append((*self)[i]);
                    }
                }
                return associations;
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
    void __setitem__(const int& index, const GModelAssociation& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "Association index",
                                           index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GModelAssociation& val) {
        (*self)[name] = val;
        return;
    }
    GModelAssociations copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = {'xml': xml}
        return state
    def __setstate__(self, state):
        self.__init__(state['xml'])
}
};
