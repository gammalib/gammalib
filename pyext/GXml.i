/***************************************************************************
 *                           GXml.i - XML class                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
    virtual ~GXml(void);

    // Methods
    void         clear(void);
    GXml*        clone(void) const;
    std::string  classname(void) const;
    int          size(void) const;
    bool         is_empty(void) const;
    GXmlNode*    set(const int& index, const GXmlNode& node);
    GXmlNode*    append(const GXmlNode& node);
    GXmlElement* append(const std::string& segment);
    GXmlNode*    insert(const int& index, const GXmlNode& node);
    void         remove(const int& index);
    void         reserve(const int& num);
    void         extend(const GXmlNode& node);
    int          elements(void) const;
    int          elements(const std::string& name) const;
    GXmlElement* element(const int& index);
    GXmlElement* element(const std::string& name);
    GXmlElement* element(const std::string& name, const int& index);
    void         load(const GFilename& filename);
    void         save(const GFilename& filename) const;
    void         read(const GUrl& url);
    void         write(GUrl& url, const int& indent = 0) const;
};


/***********************************************************************//**
 * @brief GXml class extension
 ***************************************************************************/
%extend GXml {
    GXmlNode* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index,
                                           0, self->size()-1);
        }
    }
    void __setitem__(const int& index, const GXmlNode& node) {
        if (index >= 0 && index < self->size()) {
            self->set(index, node);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index,
                                           0, self->size()-1);
        }
    }
    int __len__() {
        return (self->size());
    }
    GXml copy() {
        return (*self);
    }
};
