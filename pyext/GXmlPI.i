/***************************************************************************
 *                 GXmlPI.i - XML PI node class definition                 *
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
 * @file GXmlPI.i
 * @brief XML PI node class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlPI.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GXmlPI
 *
 * @brief XML Processing Instruction node class
 ***************************************************************************/
class GXmlPI : public GXmlNode {
public:
    // Constructors and destructors
    GXmlPI(void);
    GXmlPI(const GXmlPI& node);
    explicit GXmlPI(const std::string& segment);
    virtual ~GXmlPI(void);

    // Implemented pure virtual base class methods
    virtual void       clear(void);
    virtual GXmlPI*    clone(void) const;
    virtual void       write(GUrl& url, const int& indent = 0) const;
    virtual NodeType   type(void) const;

    // Other methods
    const std::string& pi(void) const;
    void               pi(const std::string& pi);
};


/***********************************************************************//**
 * @brief GXmlPI class extension
 ***************************************************************************/
%extend GXmlPI {
    char *__str__() {
        return tochar(self->print());
    }
    GXmlPI copy() {
        return (*self);
    }
};
