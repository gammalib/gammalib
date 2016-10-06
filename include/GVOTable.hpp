/***************************************************************************
 *                      GVOTable.hpp - VO table class                      *
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
 * @file GVOTable.hpp
 * @brief VO table class definition
 * @author Thierry Louge
 */

#ifndef GVOTable_HPP
#define GVOTable_HPP

/* __ Definitions ________________________________________________________ */

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GXml.hpp"

/* __ Forward declarations _______________________________________________ */
class GFitsTable;
class GFitsTableCol;
class GXmlElement;


/***********************************************************************//**
 * @class GVOTable
 *
 * @brief VOTable class
 *
 * This class implements a VOTable for exchanges through VO-compatible
 * applications. The class implements IVOA standard Recommendation 2013-09-20
 * VOTable1.3.
 *
 * See
 * http://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.pdf
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
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GVOTable& table);
    void        free_members(void);
    GXmlElement field_from_fits_column(const GFitsTableCol& column) const;
    GXmlElement data_from_fits_table(const GFitsTable& table) const;

    // Protected members
    GXml        m_xml;         //!< VO table
    std::string m_name;        //!< VO table name
    std::string m_resource;    //!< VO resource name
    std::string m_description; //!< VO table description
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


/***********************************************************************//**
 * @brief Return VO table XML file
 *
 * @return VO table XML file.
 ***************************************************************************/
inline
const GXml& GVOTable::xml(void) const
{
    return (m_xml);
}


/***********************************************************************//**
 * @brief Return VO table name
 *
 * @return VO table name string.
 ***************************************************************************/
inline
const std::string& GVOTable::name(void) const
{
    return (m_name);
}

#endif /* GVOTable_HPP */
