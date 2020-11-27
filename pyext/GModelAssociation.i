/***************************************************************************
 *              GModelAssociation.i - Model association class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GModelAssociation.i
 * @brief Model association class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelAssociation.hpp"
%}


/***********************************************************************//**
 * @class GModelAssociation
 *
 * @brief Model association class
 *
 * The GModelAssociation class stores association information for a model.
 ***************************************************************************/
class GModelAssociation : public GBase {

public:
    // Constructors and destructors
    GModelAssociation(void);
    explicit GModelAssociation(const std::string& name);
    explicit GModelAssociation(const GXmlElement& xml);
    GModelAssociation(const GModelAssociation& association);
    virtual ~GModelAssociation(void);

    // Methods
    void               clear(void);
    GModelAssociation* clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    bool               is_empty(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    const std::string& value(const std::string& name) const;
    const std::string& error(const std::string& name) const;
    void               property(const std::string& name,
                                const std::string& value,
                                const std::string& error = "");
    void               read(const GXmlElement& xml);
    void               write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelAssociation class extension
 ***************************************************************************/
%extend GModelAssociation {
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = {'xml': xml}
        return state
    def __setstate__(self, state):
        self.__init__(state['xml'])
}
};
