/***************************************************************************
 *                 GXmlNode.i - Abstract XML node base class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GXmlNode.i
 * @brief Abstract XML node base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlNode.hpp"
#include "GTools.hpp"
#include "GException.hpp"
%}


/***********************************************************************//**
 * @class GXmlNode
 *
 * @brief Abstract XML node base class
 ***************************************************************************/
class GXmlNode : public GContainer {
public:
    // Constructors and destructors
    GXmlNode(void);
    GXmlNode(const GXmlNode& node);
    virtual ~GXmlNode(void);

    // Public enumerators
    enum NodeType {
        NT_DOCUMENT,
        NT_ELEMENT,
        NT_COMMENT,
        NT_UNKNOWN,
        NT_TEXT,
        NT_DECLARATION,
        NT_PI,
        NT_TYPECOUNT
    };

    // Methods
    virtual void         clear(void) = 0;
    virtual GXmlNode*    clone(void) const = 0;
    virtual std::string  classname(void) const = 0;
    virtual int          size(void) const;
    virtual bool         is_empty(void) const;
    virtual GXmlNode*    set(const int& index, const GXmlNode& node);
    virtual GXmlNode*    append(const GXmlNode& node);
    virtual GXmlElement* append(const std::string& segment);
    virtual GXmlNode*    insert(const int& index, const GXmlNode& node);
    virtual void         remove(const int& index);
    virtual void         reserve(const int& num);
    virtual void         extend(const GXmlNode& node);
    GXmlNode*            parent(void) const;
    void                 parent(GXmlNode* parent);
    GFilename            filename(void) const;
    virtual int          elements(void) const;
    virtual int          elements(const std::string& name) const;
    virtual GXmlElement* element(const int& index);
    virtual GXmlElement* element(const std::string& name);
    virtual GXmlElement* element(const std::string& name, const int& index);
    virtual void         write(GUrl& url, const int& indent) const = 0;
    virtual NodeType     type(void) const = 0;
};


/***********************************************************************//**
 * @brief GXmlNode class extension
 ***************************************************************************/
%extend GXmlNode {
    char *__str__() {
        return gammalib::tochar(self->print(NORMAL, 0));
    }
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
%pythoncode {
    def __getstate__(self):
        state = (self.parent(), tuple([x for x in self]))
        return state
    def __setstate__(self, state):
        self.__init__()
        self.reserve(len(state[1]))
        if state[0] != None:
            self.parent(state[0])
        for x in state[1]:
            self.append(x)
}
};
