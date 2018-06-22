/***************************************************************************
 *                  GXmlComment.i - XML comment node class                 *
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
 * @file GXmlComment.i
 * @brief XML comment class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlComment.hpp"
%}


/***********************************************************************//**
 * @class GXmlComment
 *
 * @brief XML comment node class
 ***************************************************************************/
class GXmlComment : public GXmlNode {
public:
    // Constructors and destructors
    GXmlComment(void);
    GXmlComment(const GXmlComment& node);
    explicit GXmlComment(const std::string& segment);
    virtual ~GXmlComment(void);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GXmlComment* clone(void) const;
    virtual std::string  classname(void) const;
    virtual void         write(GUrl& url, const int& indent = 0) const;
    virtual NodeType     type(void) const;

    // Other methods
    const std::string&   comment(void) const;
    void                 comment(const std::string& comment);
};


/***********************************************************************//**
 * @brief GXmlComment class extension
 ***************************************************************************/
%extend GXmlComment {
    GXmlComment copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.comment(),)
        return state
    def __setstate__(self, state):
        self.__init__()
        self.comment(state[0])
}
};
