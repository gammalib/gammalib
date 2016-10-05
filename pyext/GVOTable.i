/***************************************************************************
 *                       GVOTable.i - VOTable class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Thierry Louge                               *
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
 * @brief GVOTable class interface definition
 * @author Thierry Louge
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GVOTable.hpp"
#include "GXml.hpp"
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
    GVOTable(const GVOTable& votable);
    GVOTable(const std::string& filename);
    virtual ~GVOTable(void);

    // Methods
    void        	clear(void);
    GVOTable*     	clone(void) const;
    std::string 	classname(void) const;
    std::string		print(const GChatter& chatter = NORMAL) const;
    void 		open_votable(void);
    void 		fill_tabledata(const std::string& data);
    void 		fill_fields(const std::string& name, const std::string& ucd,
				const std::string& id, const std::string& datatype, 
				const std::string& width,const std::string& precision,
				const std::string& unit,const std::string& description);
    void 		close_votable(void);
    void 		init_tabledata(void);
    void 		close_tabledata(void);
};


/***********************************************************************//**
 * @brief GVOTable class extension
 ***************************************************************************/
%extend GVOTable {
    GVOTable copy() {
        return (*self);
    }
};
