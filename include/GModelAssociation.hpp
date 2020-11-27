/***************************************************************************
 *              GModelAssociation.hpp - Model association class            *
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
 * @file GModelAssociation.hpp
 * @brief Model association class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELASSOCIATION_HPP
#define GMODELASSOCIATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GXmlElement;


/***********************************************************************//**
 * @class GModelAssociation
 *
 * @brief Model association class
 *
 * The GModelAssociation class stores association information for a model.
 * Each association has a name and an arbitrary number properties.
 ***************************************************************************/
class GModelAssociation : public GBase {

public:
    // Constructors and destructors
    GModelAssociation(void);
    explicit GModelAssociation(const std::string& name);
    explicit GModelAssociation(const GXmlElement& xml);
    GModelAssociation(const GModelAssociation& association);
    virtual ~GModelAssociation(void);
 
    // Operators
    GModelAssociation& operator=(const GModelAssociation& association);

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
    std::string        print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelAssociation& association);
    void         free_members(void);
    int          get_index(const std::string& name) const;
    GXmlElement* get_property_xml(GXmlElement&       xml,
                                  const std::string& name) const;

    // Protected data members
    std::string              m_name;    //!< Association name
    std::vector<std::string> m_names;   //!< Property names
    std::vector<std::string> m_values;  //!< Property values
    std::vector<std::string> m_errors;  //!< Property errors
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelAssociation").
 ***************************************************************************/
inline
std::string GModelAssociation::classname(void) const
{
    return ("GModelAssociation");
}


/***********************************************************************//**
 * @brief Return number of association properties
 *
 * @return Number of association properties.
 *
 * Returns the number of association properties.
 ***************************************************************************/
inline
int GModelAssociation::size(void) const
{
    return (int)m_names.size();
}


/***********************************************************************//**
 * @brief Signals if there are no association properties
 *
 * @return True if the association has no properties, false otherwise.
 *
 * Signals if the association has no properties.
 ***************************************************************************/
inline
bool GModelAssociation::is_empty(void) const
{
    return (m_names.empty());
}


/***********************************************************************//**
 * @brief Return association name
 *
 * @return Association name.
 ***************************************************************************/
inline
const std::string& GModelAssociation::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set association name
 *
 * @param[in] name Association name.
 ***************************************************************************/
inline
void GModelAssociation::name(const std::string& name)
{
    m_name = name;
    return;
}
#endif /* GMODELASSOCIATION_HPP */
