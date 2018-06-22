/***************************************************************************
 *                    GXmlText.i - XML text node class                     *
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
 * @file GXmlText.i
 * @brief XML text node class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlText.hpp"
%}


/***********************************************************************//**
 * @class GXmlText
 *
 * @brief XML text node class
 ***************************************************************************/
class GXmlText : public GXmlNode {
public:
    // Constructors and destructors
    GXmlText(void);
    GXmlText(const GXmlText& node);
    explicit GXmlText(const std::string& text);
    virtual ~GXmlText(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GXmlText*   clone(void) const;
    virtual std::string classname(void) const;
    virtual void        write(GUrl& url, const int& indent = 0) const;
    virtual NodeType    type(void) const;

    // Other methods
    const std::string&  text(void) const;
    void                text(const std::string& text);
};


/***********************************************************************//**
 * @brief GXmlText class extension
 ***************************************************************************/
%extend GXmlText {
    GXmlText copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.text(),)
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
