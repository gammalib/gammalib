/***************************************************************************
 *       GModelAssociations.hpp - Model associations container class       *
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
 * @file GModelAssociations.hpp
 * @brief Model associations container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELASSOCIATIONS_HPP
#define GMODELASSOCIATIONS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */
class GModelAssociation;
class GXmlElement;


/***********************************************************************//**
 * @class GModelAssociations
 *
 * @brief Model associations container class
 ***************************************************************************/
class GModelAssociations : public GContainer {

public:
    // Constructors and destructors
    GModelAssociations(void);
    explicit GModelAssociations(const GXmlElement& xml);
    GModelAssociations(const GModelAssociations& associations);
    virtual ~GModelAssociations(void);

    // Operators
    GModelAssociations&      operator=(const GModelAssociations& associations);
    GModelAssociation&       operator[](const int& index);
    const GModelAssociation& operator[](const int& index) const;
    GModelAssociation&       operator[](const std::string& name);
    const GModelAssociation& operator[](const std::string& name) const;

    // Methods
    void                     clear(void);
    GModelAssociations*      clone(void) const;
    std::string              classname(void) const;
    GModelAssociation&       at(const int& index);
    const GModelAssociation& at(const int& index) const;
    int                      size(void) const;
    bool                     is_empty(void) const;
    GModelAssociation&       append(const GModelAssociation& association);
    GModelAssociation&       insert(const int&               index,
                                    const GModelAssociation& association);
    GModelAssociation&       insert(const std::string&       name,
                                    const GModelAssociation& association);
    void                     remove(const int& index);
    void                     remove(const std::string& name);
    void                     reserve(const int& num);
    void                     extend(const GModelAssociations& associations);
    bool                     contains(const std::string& name) const;
    void                     read(const GXmlElement& xml);
    void                     write(GXmlElement& xml) const;
    std::string              print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelAssociations& associations);
    void         free_members(void);
    int          get_index(const std::string& name) const;
    GXmlElement* get_association_xml(GXmlElement&       xml,
                                     const std::string& name) const;

    // Proteced members
    std::vector<GModelAssociation> m_associations;  //!< List of associations
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelAssociations").
 ***************************************************************************/
inline
std::string GModelAssociations::classname(void) const
{
    return ("GModelAssociations");
}


/***********************************************************************//**
 * @brief Return reference to association
 *
 * @param[in] index Association index [0,...,size()-1].
 *
 * Returns a reference to the association with the specified @p index.
 ***************************************************************************/
inline
GModelAssociation& GModelAssociations::operator[](const int& index)
{
    return (m_associations[index]);
}


/***********************************************************************//**
 * @brief Return reference to association (const version)
 *
 * @param[in] index Association index [0,...,size()-1].
 *
 * Returns a const reference to the association with the specified @p index.
 ***************************************************************************/
inline
const GModelAssociation& GModelAssociations::operator[](const int& index) const
{
    return (m_associations[index]);
}


/***********************************************************************//**
 * @brief Return number of associations in container
 *
 * @return Number of associations in container.
 *
 * Returns the number of associations in the container.
 ***************************************************************************/
inline
int GModelAssociations::size(void) const
{
    return (int)m_associations.size();
}


/***********************************************************************//**
 * @brief Signals if there are no associations in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the association container does not contain any association.
 ***************************************************************************/
inline
bool GModelAssociations::is_empty(void) const
{
    return (m_associations.empty());
}


/***********************************************************************//**
 * @brief Reserves space for associations in container
 *
 * @param[in] num Number of associations
 *
 * Reserves space for @p num associations in the container.
 ***************************************************************************/
inline
void GModelAssociations::reserve(const int& num)
{
    m_associations.reserve(num);
    return;
}

#endif /* GMODELASSOCIATIONS_HPP */
