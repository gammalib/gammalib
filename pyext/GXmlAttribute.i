/***************************************************************************
 *                   GXmlAttribute.i - XML attribute class                 *
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
 * @file GXmlAttribute.i
 * @brief XML arrtibute class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlAttribute.hpp"
%}


/***********************************************************************//**
 * @class GXmlAttribute
 *
 * @brief XML attribute class
 ***************************************************************************/
class GXmlAttribute : public GBase {
public:
    // Constructors and destructors
    GXmlAttribute(void);
    GXmlAttribute(const GXmlAttribute& attr);
    GXmlAttribute(const std::string& name, const std::string& value);
    virtual ~GXmlAttribute(void);

    // Methods
    void               clear(void);
    GXmlAttribute*     clone(void) const;
    std::string        classname(void) const;
    void               write(GUrl& url) const;
    const std::string& name(void) const;
    std::string        value(void) const;
    void               name(const std::string& name);
    void               value(std::string value);
};


/***********************************************************************//**
 * @brief GXmlAttribute class extension
 ***************************************************************************/
%extend GXmlAttribute {
    GXmlAttribute copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.name(), self.value())
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1])
}
};
