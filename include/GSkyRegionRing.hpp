/***************************************************************************
 *              GSkyRegionRing.hpp - Ring sky region class                 *
 * ----------------------------------------------------------------------- *
 * copyright (C) 2013 by Maria Krause, Anneli Schulz                       *
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
 * @file GSkyRegionRing.hpp
 * @brief Ring sky region class interface definition
 * @author Maria Krause, Anneli Schulz
 */

#ifndef GSKYREGIONRING_HPP
#define GSKYREGIONRING_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GMath.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegion.hpp"


/***********************************************************************//**
 * @class GSkyRegionRing
 *
 * @brief Interface for the Ring sky region class
 *
 * This class provides an implementation for a circular sky region. The sky
 * region is defined by an array of parameters, the meaning of which is
 * specific to the derived class where the region type or shape is defined.
 *
 * The class holds several properties such as solid angle subtended by the
 * region and computed through internal method compute_solid().
 *
 * To be clarified:
 * - Do we want a member relating the region to an observation run ?
 * - Constructor and read/write using XML may not be needed if we use DS9
 *   region file format ?
 *
 ***************************************************************************/
class GSkyRegionRing : public GSkyRegion {

public:
    // Constructors and destructors
    GSkyRegionRing(void);
    GSkyRegionRing(GSkyDir& centre, const double& radius1, const double& radius2);
    GSkyRegionRing(const double& ra, const double& dec, const double& radius1, const double& radius2);
    explicit GSkyRegionRing(const std::string& line);
    GSkyRegionRing(const GSkyRegionRing& region);
    virtual ~GSkyRegionRing(void);

    // Operators
    GSkyRegionRing& operator=(const GSkyRegionRing& region);

    // Implemented methods
    void              clear(void);
    GSkyRegionRing* clone(void) const;
    const double&     radius1(void) const;
    void              radius1(const double& radius1);
    const double&     radius2(void) const;
    void              radius2(const double& radius2);
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
    void copy_members(const GSkyRegionRing& region);
    void free_members(void);
    void compute_solid_angle(void);

    // Protected members
    GSkyDir	m_centre;   //!< Centre or reference point of the region
    double 	m_radius1; 	//!< Radius of inner ring the region [deg]
    double  m_radius2;  //!< Radius of outer ring the region [deg]
};


/***********************************************************************//**
 * @brief Return ring region inner (radius1) and outer (radius2) region (in degrees)
 *
 * @return Region radius [deg].
 *
 * Returns the region radius in degrees.
 ***************************************************************************/
inline
const double& GSkyRegionRing::radius1(void) const
{
    return (m_radius1);
}

inline
const double& GSkyRegionRing::radius2(void) const
{
    return (m_radius2);
}


/***********************************************************************//**
 * @brief Return circular region centre
 *
 * @return Region centre.
 *
 * Returns the region centre.
 ***************************************************************************/
inline
const GSkyDir& GSkyRegionRing::centre(void) const
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
void GSkyRegionRing::centre(const GSkyDir& dir)
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
void GSkyRegionRing::centre(const double& ra, const double& dec)
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
double GSkyRegionRing::ra(void) const
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
double GSkyRegionRing::dec(void) const
{
    return (m_centre.dec_deg());
}

#endif /* GSKYREGIONRing_HPP */
