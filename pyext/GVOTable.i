/***************************************************************************
 *                       GVOTable.i - VOTable class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Thierry Louge                               *
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
 * @file GVOTable.i
 * @brief VO table class definition
 * @author Thierry Louge
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GVOTable.hpp"
%}


/***********************************************************************//**
 * @class GVOTable
 *
 * @brief VOTable class
 ***************************************************************************/
class GVOTable : public GBase {

public:
    // Constructors and destructors
    GVOTable(void);
    explicit GVOTable(const GFitsTable& table);
    GVOTable(const GVOTable& votable);
    virtual ~GVOTable(void);

    // Operators
    GVOTable& operator=(const GVOTable& votable);

    // Methods
    void               clear(void);
    GVOTable*          clone(void) const;
    std::string        classname(void) const;
    void               read(const GFitsTable& table);
    const GXml&        xml(void) const;
    const std::string& name(void) const;
};


/***********************************************************************//**
 * @brief GVOTable class extension
 ***************************************************************************/
%extend GVOTable {
    GVOTable copy() {
        return (*self);
    }
};
