/***************************************************************************
 *             GXmlElement.i - XML element node class definition           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GXmlElement.i
 * @brief XML element node class interface definition
  * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlElement.hpp"
%}


/***********************************************************************//**
 * @class GXmlElement
 *
 * @brief XML element node class
 ***************************************************************************/
class GXmlElement : public GXmlNode {
public:
    // Constructors and destructors
    GXmlElement(void);
    GXmlElement(const GXmlElement& node);
    explicit GXmlElement(const std::string& segment);
    virtual ~GXmlElement(void);

    // Methods
    virtual void         clear(void);
    virtual GXmlElement* clone(void) const;
    virtual std::string  classname(void) const;
    const std::string&   name(void) const;
    void                 name(const std::string& name);
    int                  attributes(void) const;
    const GXmlAttribute* attribute(const int& index) const;
    std::string          attribute(const std::string& name) const;
    void                 attribute(const std::string& name, const std::string& value);
    bool                 has_attribute(const std::string& name) const;
    void                 remove_attribute(const std::string& name);
    virtual void         write(GUrl& url, const int& indent = 0) const;
    virtual NodeType     type(void) const;
};


/***********************************************************************//**
 * @brief GXmlElement class extension
 ***************************************************************************/
%extend GXmlElement {
    GXmlNode* __getitem__(const int& index) {
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
    void __setitem__(const int& index, const GXmlNode& node) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            self->set(index, node);
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            self->set(self->size()+index, node);
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GXmlElement* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GXmlElement* element = new GXmlElement;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        element->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        element->append(*(*self)[i]);
                    }
                }
                return element;
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
    GXmlElement copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (gammalib.GXmlNode.__getstate__(self), self.name(),
                 tuple([self.attribute(i) for i in range(self.attributes())]))
        return state
    def __setstate__(self, state):
        gammalib.GXmlNode.__setstate__(self, state[0])
        self.__init__()
        self.name(state[1])
        for i in range(len(state[2])):
            name  = state[2][i].name()
            value = state[2][i].value()
            self.attribute(name, value)
}
};
