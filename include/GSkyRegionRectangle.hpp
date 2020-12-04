/***************************************************************************
 *         GSkyRegionRectangle.hpp - Rectangular sky region class          *
 * ----------------------------------------------------------------------- *
 * copyright (C) 2019-2020 by Andreas Specovius                            *
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
 * @file GSkyRegionRectangle.hpp
 * @brief Rectangular sky region class interface definition
 * @author Andreas Specovius
 */

#ifndef GSKYREGIONRECTANGLE_HPP
#define GSKYREGIONRECTANGLE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GSkyDir.hpp"
#include "GSkyRegion.hpp"

/* __ Forward declarations _______________________________________________ */
class GSkyPixel;


/***********************************************************************//**
 * @class GSkyRegionRectangle
 *
 * @brief Interface for the rectangular sky region class
 *
 * This class provides an implementation for a rectangular sky region. The sky
 * region is defined by an array of parameters, the meaning of which is
 * specific to the derived class where the region type or shape is defined.
 *
 * The class holds several properties such as solid angle subtended by the
 * region and computed through internal method compute_solid().
 *
 * The position angle counts counter-clockwise from celestial North and is
 * aligned with the height axis of the rectangle. Note, that height and width
 * describe the diameter of the region.
 *
 * Sky directions can be transformed to the local/global coordinate system via
 * the methods transform_to_local() / transform_to_global().
 * The position (in global coordinates) of the 4 corners of the rectangle can
 * be received via the method get_corner(index).
 *
 ***************************************************************************/
class GSkyRegionRectangle : public GSkyRegion {

public:
    // Constructors and destructors
    GSkyRegionRectangle(void);
    GSkyRegionRectangle(const GSkyDir& centre,
                        const double&  width,
                        const double&  height,
                        const double&  posang);
    GSkyRegionRectangle(const double& ra,
                        const double& dec,
                        const double& width,
                        const double& height,
                        const double& posang);
    explicit GSkyRegionRectangle(const std::string& line);
    GSkyRegionRectangle(const GSkyRegionRectangle& region);
    virtual ~GSkyRegionRectangle(void);

    // Operators
    GSkyRegionRectangle& operator=(const GSkyRegionRectangle& region);

    // Implemented methods
    void                 clear(void);
    GSkyRegionRectangle* clone(void) const;
    std::string          classname(void) const;
    const GSkyDir&       centre(void) const;
    void                 centre(const GSkyDir& centre);
    void                 centre(const double& ra,const double& dec);
    double               ra(void) const;
    double               dec(void) const;
    const double&        width(void) const;
    void                 width(const double& width);
    const double&        height(void) const;
    void                 height(const double& height);
    const double&        posang(void) const;
    void                 posang(const double& posang);
    void                 read(const std::string& line);
    std::string          write(void) const;
    bool                 contains(const GSkyDir& dir) const;
    bool                 contains(const GSkyRegion& reg) const;
    bool                 overlaps(const GSkyRegion& reg) const;
    GSkyDir              corner(const int& index) const;
    std::string          print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void      init_members(void);
    void      copy_members(const GSkyRegionRectangle& region);
    void      free_members(void);
    void      compute_solid_angle(void);
    bool      contains(const GSkyPixel& local) const;
    GSkyPixel dir_to_local(const GSkyDir& dir) const;
    GSkyDir   local_to_dir(const GSkyPixel& local) const;

    // Protected members
    GSkyDir m_centre; //!< Centre or reference point of the region
    double  m_width;  //!< Width of the region (degrees)
    double  m_height; //!< Height of the region (degrees)
    double  m_posang; //!< Position angle, counterclockwise from North (degrees)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyRegionRectangle").
 ***************************************************************************/
inline
std::string GSkyRegionRectangle::classname(void) const
{
    return ("GSkyRegionRectangle");
}


/***********************************************************************//**
 * @brief Return rectangular region centre
 *
 * @return Region centre.
 *
 * Returns the region centre.
 ***************************************************************************/
inline
const GSkyDir& GSkyRegionRectangle::centre(void) const
{
    return (m_centre);
}


/***********************************************************************//**
 * @brief Set rectangular region centre
 *
 * @param[in] dir Region centre.
 *
 * Sets the centre of the rectangular region to the specified sky direction.
 ***************************************************************************/
inline
void GSkyRegionRectangle::centre(const GSkyDir& dir)
{
    m_centre = dir;
    return;
}


/***********************************************************************//**
 * @brief Set rectangular region centre Right Ascension and Declincation
 *
 * @param[in] ra Right Ascension (degrees).
 * @param[in] dec Declination (degrees).
 *
 * Sets the centre of the rectangular region to the specified Right Ascension
 * and Declination.
 ***************************************************************************/
inline
void GSkyRegionRectangle::centre(const double& ra, const double& dec)
{
    m_centre.radec_deg(ra, dec);
    return;
}


/***********************************************************************//**
 * @brief Return rectangular region centre Right Ascension
 *
 * @return Region centre Right Ascension (degrees).
 *
 * Returns the region centre Right Ascension in degrees.
 ***************************************************************************/
inline
double GSkyRegionRectangle::ra(void) const
{
    return (m_centre.ra_deg());
}


/***********************************************************************//**
 * @brief Return rectangular region centre Declination
 *
 * @return Region centre Declination (degrees).
 *
 * Returns the region centre Declination in degrees.
 ***************************************************************************/
inline
double GSkyRegionRectangle::dec(void) const
{
    return (m_centre.dec_deg());
}


/***********************************************************************//**
 * @brief Return region width extension (in degrees)
 *
 * @return Region width (degrees).
 *
 * Returns the region width extension in degrees.
 ***************************************************************************/
inline
const double& GSkyRegionRectangle::width(void) const
{
    return (m_width);
}


/***********************************************************************//**
 * @brief Return region height extension (in degrees)
 *
 * @return Region height (degrees).
 *
 * Returns the region height extension in degrees.
 ***************************************************************************/
inline
const double& GSkyRegionRectangle::height(void) const
{
    return (m_height);
}


/***********************************************************************//**
 * @brief Return region position angle (in degrees)
 *
 * @return Region position angle (degrees).
 *
 * Returns the region position angle in degrees. The position angle is
 * counted counterclockwise from North.
 ***************************************************************************/
inline
const double& GSkyRegionRectangle::posang(void) const
{
    return (m_posang);
}


/***********************************************************************//**
 * @brief Set position angle of rectangular region
 *
 * @param[in] posang Position angle (degrees).
 *
 * Sets the position angle of the rectangular sky region. The position angle
 * is counted counterclockwise from celestial North and is aligned to the
 * rectangle height.
 ***************************************************************************/
inline
void GSkyRegionRectangle::posang(const double& posang)
{
    m_posang = posang;
    return;
}

#endif /* GSKYREGIONRECTANGLE_HPP */
