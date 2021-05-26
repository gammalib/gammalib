/***************************************************************************
 *                           GXml.i - XML class                            *
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
 * @file GXml.i
 * @brief XML class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXml.hpp"
#include "GTools.hpp"
#include "GFilename.hpp"
%}


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class
 ***************************************************************************/
class GXml : public GContainer {
public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    explicit GXml(const std::string& filename);
    explicit GXml(const GXmlDocument& root);
    virtual ~GXml(void);

    // Methods
    void                clear(void);
    GXml*               clone(void) const;
    std::string         classname(void) const;
    int                 size(void) const;
    bool                is_empty(void) const;
    GXmlNode*           set(const int& index, const GXmlNode& node);
    GXmlNode*           append(const GXmlNode& node);
    GXmlElement*        append(const std::string& segment);
    GXmlNode*           insert(const int& index, const GXmlNode& node);
    void                remove(const int& index);
    void                reserve(const int& num);
    void                extend(const GXmlNode& node);
    int                 elements(void) const;
    int                 elements(const std::string& name) const;
    GXmlElement*        element(const int& index);
    GXmlElement*        element(const std::string& name);
    GXmlElement*        element(const std::string& name, const int& index);
    const GXmlDocument& root(void) const;
    void                root(const GXmlDocument& root);
    void                load(const GFilename& filename);
    void                save(const GFilename& filename) const;
    void                read(const GUrl& url);
    void                write(GUrl& url, const int& indent = 0) const;
};


/***********************************************************************//**
 * @brief GXml class extension
 ***************************************************************************/
%extend GXml {
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
            throw GException::out_of_range("__getitem__(int)",
                                           "Node index",
                                           index, self->size());
        }
    }
    GXml* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GXml* xml = new GXml;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        xml->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        xml->append(*(*self)[i]);
                    }
                }
                return xml;
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
            throw GException::out_of_range("__setitem__(int)",
                                           "Node index",
                                           index, self->size());
        }
    }
    GXml copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.root(),)
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
