/***************************************************************************
 *         GSkyRegion.hpp - Abstract virtual sky region base class         *
 * ----------------------------------------------------------------------- *
 * copyright (C) 2013-2017 by Pierrick Martin                              *
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
 * @file GSkyRegion.hpp
 * @brief Abstract sky region base class interface definition
 * @author Pierrick Martin
 */

#ifndef GSKYREGION_HPP
#define GSKYREGION_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GSkyPixel;
class GSkyDir;


/***********************************************************************//**
 * @class GSkyRegion
 *
 * @brief Abstract interface for the sky region class
 *
 * This class provides an abstract interface for a sky region. The sky region 
 * is defined by an array of parameters names and values specific to the 
 * derived class where the region type or shape is implemented. 
 * 
 * Accessible information elements are:
 * - Type of the region (circle, rectangle,...)
 * - Solid angle subtended by the region
 * 
 * The object can be initialised from a string in the DS9 region file format,
 * and return its description as a string in the DS9 region file format. The 
 * input/ouput to/from files is handled in the container class GSkyRegions.
 ***************************************************************************/
class GSkyRegion : public GBase {

public:
    // Constructors and destructors
    GSkyRegion(void);
    GSkyRegion(const GSkyRegion& region); 
    virtual ~GSkyRegion(void);

    // Operators
    virtual GSkyRegion& operator=(const GSkyRegion& region);
	
    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GSkyRegion* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual void        read(const std::string& regstring) = 0;
    virtual std::string write(void) const = 0;
    virtual bool        contains(const GSkyDir& dir) const = 0;
    virtual bool        overlaps(const GSkyRegion& reg) const = 0;
    virtual bool        contains(const GSkyRegion& reg) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;	
	
    // Implemented methods
    const std::string& type(void) const;
    const std::string& name(void) const;
    const double&      solidangle(void) const;
    void               type(const std::string& type);
    void               name(const std::string& name);
    void               solidangle(const double& solidangle);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSkyRegion& region);
    void free_members(void);

    // Protected members
    std::string m_type;    //!< Type of the region (circle, rectangle,...)
    std::string m_name;    //!< Name of the region
    double      m_solid;   //!< Solid angle subtended by the region (sr)
};


/***********************************************************************//**
 * @brief Return region name
 *
 * @return Region name
 *
 * Returns the region name.
 ***************************************************************************/
inline
const std::string& GSkyRegion::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set region name
 *
 * @param[in] name Region name.
 *
 * Sets the region name.
 ***************************************************************************/
inline
void GSkyRegion::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return region type
 *
 * @return Region type
 *
 * Returns the region type.
 ***************************************************************************/
inline
const std::string& GSkyRegion::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Set region type
 *
 * @param[in] type Region type.
 *
 * Sets the region type.
 ***************************************************************************/
inline
void GSkyRegion::type(const std::string& type)
{
    m_type = type;
    return;
}


/***********************************************************************//**
 * @brief Return solid angle of region
 *
 * @return Solid angle of region.
 *
 * Returns the solid angle subtended by the region (in steradians).
 ***************************************************************************/
inline
const double& GSkyRegion::solidangle(void) const
{
    return (m_solid);
}


/***********************************************************************//**
 * @brief Set solid angle of region
 *
 * @param[in] solidangle Solid angle of region (in steradians).
 *
 * Sets the solid angle subtended by the region (in steradians).
 ***************************************************************************/
inline
void GSkyRegion::solidangle(const double& solidangle)
{
    m_solid = solidangle;
    return;
}

#endif /* GSKYREGION_HPP */
