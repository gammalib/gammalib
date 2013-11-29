/***************************************************************************
 *             GLATRoi.hpp - Fermi/LAT region of interest class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GLATRoi.hpp
 * @brief Fermi/LAT region of interest class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATROI_HPP
#define GLATROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"
#include "GLATInstDir.hpp"


/***********************************************************************//**
 * @class GLATRoi
 *
 * @brief Interface for the LAT region of interest class.
 *
 * The Fermi-LAT region of interest class defines the region of photon
 * arrival directions that is used for unbinned data analysis. A circular
 * Fermi-LAT has been implemented.
 ***************************************************************************/
class GLATRoi : public GRoi {

public:
    // Constructors and destructors
    GLATRoi(void);
    GLATRoi(const GLATInstDir& centre, const double& radius);
    GLATRoi(const GLATRoi& roi);
    virtual ~GLATRoi(void);

    // Operators
    GLATRoi& operator=(const GLATRoi& roi);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GLATRoi*    clone(void) const;
    virtual bool        contains(const GEvent& event) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GLATInstDir&  centre(void) const;
    const double&       radius(void) const;
    void                centre(const GLATInstDir& centre);
    void                radius(const double& radius);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATRoi& roi);
    void free_members(void);
    
    // Protected members
    GLATInstDir m_centre;   //!< Centre of ROI in instrument coordinates
    double      m_radius;   //!< Radius of ROI in degrees
};


/***********************************************************************//**
 * @brief Returns region of interest centre
 *
 * @return Region of interest centre sky direction.
 *
 * Returns the sky direction of the region of interest centre.
 ***************************************************************************/
inline
const GLATInstDir& GLATRoi::centre(void) const
{
    return (m_centre);
}


/***********************************************************************//**
 * @brief Returns radius of region of interest in degrees
 *
 * @return Region of interest radius (degrees).
 *
 * Returns the radius of the region of interest in degrees.
 ***************************************************************************/
inline
const double& GLATRoi::radius(void) const
{
    return (m_radius);
}


/***********************************************************************//**
 * @brief Set region of interest centre
 *
 * @param[in] centre Region of interest centre sky direction.
 *
 * Set the sky direction of the region of interest centre.
 ***************************************************************************/
inline
void GLATRoi::centre(const GLATInstDir& centre)
{
    m_centre = centre;
    return;
}


/***********************************************************************//**
 * @brief Set radius of region of interest
 *
 * @param[in] radius Region of interest radius (degrees).
 *
 * Set the radius of the region of interest.
 ***************************************************************************/
inline
void GLATRoi::radius(const double& radius)
{
    m_radius = radius;
    return;
}

#endif /* GLATROI_HPP */
