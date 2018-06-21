/***************************************************************************
 *                 GXmlDocument.i - XML document node class                *
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
 * @file GXmlDocument.i
 * @brief XML document class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlDocument.hpp"
%}


/***********************************************************************//**
 * @class GXmlDocument
 *
 * @brief XML document node class
 ***************************************************************************/
class GXmlDocument : public GXmlNode {
public:
    // Constructors and destructors
    GXmlDocument(void);
    GXmlDocument(const GFilename&   filename,
                 const std::string& version,
                 const std::string& encoding,
                 const std::string& standalone);
    GXmlDocument(const GXmlDocument& node);
    virtual ~GXmlDocument(void);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GXmlDocument* clone(void) const;
    virtual std::string   classname(void) const;
    virtual void          write(GUrl& url, const int& indent = 0) const;
    virtual NodeType      type(void) const;

    // Methods
    const GFilename& filename(void) const;
    std::string      version(void) const;
    std::string      encoding(void) const;
    std::string      standalone(void) const;
    void             filename(const GFilename& filename);
    void             version(const std::string& version);
    void             encoding(const std::string& encoding);
    void             standalone(const std::string& standalone);
};


/***********************************************************************//**
 * @brief GXmlDocument class extension
 ***************************************************************************/
%extend GXmlDocument {
    GXmlDocument copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.filename(), self.version(), self.encoding(), self.standalone())
        return state
    def __setstate__(self, state):
        self.__init__()
        self.filename(state[0])
        self.version(state[1])
        self.encoding(state[2])
        self.standalone(state[3])
}
};
