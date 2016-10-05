/***************************************************************************
 *                      GVOTable.hpp - VOTable class                     *
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
 * @file GVOTable.hpp
 * @brief VOTable class definition
 * @author Thierry Louge
 */

#ifndef GVOTable_HPP
#define GVOTable_HPP

/* __ Definitions ________________________________________________________ */

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <sys/socket.h>
#include "GBase.hpp"
#include "GXml.hpp"
#include "GXmlNode.hpp"
#include "GWcs.hpp"
/***********************************************************************//**
 * @class GVOTable
 *
 * @brief VOTable class
 *
 * This class implements a VOTable for exchanges through VO-compatible
 * applications.
 ***************************************************************************/
class GVOTable : public GBase {

public:
    // Constructors and destructors
    GVOTable(void);
    GVOTable(const GVOTable& votable);
    GVOTable(const std::string& filename);
    virtual ~GVOTable(void);

    // Operators
    GVOTable& operator=(const GVOTable& votable);

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
    // Members
    GXml* 		m_tablexml;	      //!< GXml object containing the VOTable
    std::string		m_sharedtablename;    //!< name of the VOTable copied on temp dir

protected:
    // Protected methods
    void                init_members(void);
    void                copy_members(const GVOTable& client);
    void                free_members(void);
    
private:
    //Protected members
    std::string		m_header;     //!< header of the VOTable
    std::string		m_wcs;        //!< representation of wcs in the VOTable
    std::string		m_pixels;     //!< pixel system and values in the VOTable
    std::string		m_fields;     //!< <FIELD> of the VOTable
    std::string		m_data;       //!< data part of VOTable
    std::string		m_footer;     //!< footer of the VOTable
    int 		m_nx;	      //!< Number of pixels in x axis of skymap
    int 		m_ny;	      //!< Number of pixels in y axis of skymap
    
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GVOTable").
 ***************************************************************************/
inline
std::string GVOTable::classname(void) const
{
    return ("GVOTable");
}

#endif /* GVOTable_HPP */
