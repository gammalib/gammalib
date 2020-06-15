/***************************************************************************
 *              GSkyRegionRect.hpp - circular sky region class           *
 * ----------------------------------------------------------------------- *
 * copyright (C) 2013-2015 by Michael Mayer                                *
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
 * @file GSkyRegionRect.hpp
 * @brief Circular sky region class interface definition
 * @author Michael Mayer
 */

#ifndef GSKYREGIONRECT_HPP
#define GSKYREGIONRECT_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GMath.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegion.hpp"


/***********************************************************************//**
 * @class GSkyRegionRect
 *
 * @brief Interface for the circular sky region class
 *
 * This class provides an implementation for a circular sky region. The sky
 * region is defined by an array of parameters, the meaning of which is
 * specific to the derived class where the region type or shape is defined.
 * The parameters are GModelPar objects for convenience.
 *
 * The class holds several properties such as solid angle subtended by the
 * region and computed through internal method compute_solid().
 *
 * To be clarified:
 * - Do we want a member relating the region to an observation run ?
 * - Constructor and read/write using XML may not be needed if we use DS9
 *   region file format ?
 * - Replace GModelPar by double for the parameters (GModelPar is overkill) ?
 *
 ***************************************************************************/
class GSkyRegionRect : public GSkyRegion {

public:
    // Constructors and destructors
    GSkyRegionRect(void);
    GSkyRegionRect(const GSkyDir& centre, const double& w, const double& h);
    GSkyRegionRect(const double& ra, const double& dec, const double& w, const double& h);
    explicit GSkyRegionRect(const std::string& line);
    GSkyRegionRect(const GSkyRegionRect& region);
    virtual ~GSkyRegionRect(void);

    // Operators
    GSkyRegionRect& operator=(const GSkyRegionRect& region);

    // Implemented methods
    void              clear(void);
    GSkyRegionRect*   clone(void) const;
    std::string       classname(void) const;
    double            width(void) const;
    void              width(const double& width);
    double            height(void) const;
    void              height(const double& height);
    const GSkyDir&    centre(void) const;
    void              centre(const GSkyDir& centre);
    void              centre(const double& ra,const double& dec);
    double            ra(void) const;
    double            dec(void) const;
    void              read(const std::string& line);
    std::string       write(void) const;
    bool              contains(const GSkyDir& dir) const;
    bool              contains(const GSkyRegion& reg) const;
    bool              overlaps(const GSkyRegion& reg) const;
    std::string       print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSkyRegionRect& region);
    void free_members(void);
    void compute_solid_angle(void);

    // Protected members
    GSkyDir	m_centre;   //!< Centre or reference point of the region
    double 	m_halfwidth;  //!< Half width of the region [deg]
    double  m_halfheight; //!< Half height of the region [deg]
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyRegionRect").
 ***************************************************************************/
inline
std::string GSkyRegionRect::classname(void) const
{
    return ("GSkyRegionRect");
}


/***********************************************************************//**
 * @brief Return region width extension (in degrees)
 *
 * @return Region width [deg].
 *
 * Returns the region width extension in degrees.
 ***************************************************************************/
inline
double GSkyRegionRect::width(void) const
{
    double width = 2*m_halfwidth;
    return (width);
}


/***********************************************************************//**
 * @brief Return region height extension (in degrees)
 *
 * @return Region height [deg].
 *
 * Returns the region height extension in degrees.
 ***************************************************************************/
inline
double GSkyRegionRect::height(void) const
{
    double height = 2*m_halfheight;
    return (height);
}


/***********************************************************************//**
 * @brief Return circular region centre
 *
 * @return Region centre.
 *
 * Returns the region centre.
 ***************************************************************************/
inline
const GSkyDir& GSkyRegionRect::centre(void) const
{
    return (m_centre);
}


/***********************************************************************//**
 * @brief Set circular region centre
 *
 * @param[in] dir Region centre.
 *
 * Sets the centre of the circular region to the specified sky direction.
 ***************************************************************************/
inline
void GSkyRegionRect::centre(const GSkyDir& dir)
{
    // Set centre
    m_centre = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set circular region centre Right Ascension and Declincation
 *
 * @param[in] ra Right Ascension [deg].
 * @param[in] dec Declination [deg].
 *
 * Sets the centre of the circular region to the specified Right Ascension
 * and Declination.
 ***************************************************************************/
inline
void GSkyRegionRect::centre(const double& ra, const double& dec)
{
    // Set centre values
    m_centre.radec_deg(ra,dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return circular region centre Right Ascension
 *
 * @return Region centre Right Ascension [deg].
 *
 * Returns the region centre Right Ascension in degrees.
 ***************************************************************************/
inline
double GSkyRegionRect::ra(void) const
{
    return (m_centre.ra_deg());
}


/***********************************************************************//**
 * @brief Return circular region centre Declination
 *
 * @return Region centre Declination [deg].
 *
 * Returns the region centre Declination in degrees.
 ***************************************************************************/
inline
double GSkyRegionRect::dec(void) const
{
    return (m_centre.dec_deg());
}

#endif /* GSKYREGIONRECT_HPP */
