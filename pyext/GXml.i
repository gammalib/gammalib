/***************************************************************************
 *                        GXml.i - XML class definition                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @brief XML class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXml.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class
 ***************************************************************************/
class GXml : public GBase {
public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    explicit GXml(const std::string& filename);
    virtual ~GXml(void);

    // Methods
    void         clear(void);
    GXml*        clone(void) const;
    void         append(GXmlNode* node);
    void         load(const std::string& filename);
    void         save(const std::string& filename);
    void         read(GUrl& url);
    void         write(GUrl& url, const int& indent = 0) const;
    int          children(void) const;
    GXmlNode*    child(int index) const;
    int          elements(void) const;
    int          elements(const std::string& name) const;
    GXmlElement* element(int index) const;
    GXmlElement* element(const std::string& name, int index) const;
};


/***********************************************************************//**
 * @brief GXml class extension
 ***************************************************************************/
%extend GXml {
    char *__str__() {
        return tochar(self->print(0));
    }
    GXml copy() {
        return (*self);
    }
};
