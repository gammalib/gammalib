/***************************************************************************
 *                GCTARoi.hpp - CTA region of interest class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GCTARoi.hpp
 * @brief CTA region of interest class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAROI_HPP
#define GCTAROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"
#include "GCTAInstDir.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GXmlElement;


/***********************************************************************//**
 * @class GCTARoi
 *
 * @brief Interface for the CTA region of interest class.
 *
 * The CTA region of interest class defines the region of photon arrival
 * directions that is used for unbinned data analysis. A circular ROI has
 * been implemented.
 ***************************************************************************/
class GCTARoi : public GRoi {

public:
    // Constructors and destructors
    GCTARoi(void);
    explicit GCTARoi(const GXmlElement& xml);
    GCTARoi(const GCTAInstDir& centre, const double& radius);
    GCTARoi(const GCTARoi& roi);
    virtual ~GCTARoi(void);

    // Operators
    GCTARoi& operator=(const GCTARoi& roi);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCTARoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    bool                is_valid(void) const;
    bool                contains(const GCTAInstDir& dir) const;
    const GCTAInstDir&  centre(void) const;
    const double&       radius(void) const;
    void                centre(const GCTAInstDir& centre);
    void                radius(const double& radius);
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTARoi& roi);
    void free_members(void);
    
    // Protected members
    GCTAInstDir m_centre;   //!< Centre of ROI in instrument coordinates
    double      m_radius;   //!< Radius of ROI in degrees
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTARoi").
 ***************************************************************************/
inline
std::string GCTARoi::classname(void) const
{
    return ("GCTARoi");
}


/***********************************************************************//**
 * @brief Returns region of interest centre
 *
 * @return Region of interest centre sky direction.
 *
 * Returns the sky direction of the region of interest centre.
 ***************************************************************************/
inline
const GCTAInstDir& GCTARoi::centre(void) const
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
const double& GCTARoi::radius(void) const
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
void GCTARoi::centre(const GCTAInstDir& centre)
{
    m_centre = centre;
    return;
}


/***********************************************************************//**
 * @brief Checks if RoI is valid
 *
 * @return Bool flag if ROI information is valid.
 ***************************************************************************/
inline
bool GCTARoi::is_valid(void) const
{
    return (m_radius > 0.0);
}


#endif /* GCTAROI_HPP */
