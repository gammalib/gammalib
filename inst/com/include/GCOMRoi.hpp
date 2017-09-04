/***************************************************************************
 *              GCOMRoi.hpp - COMPTEL region of interest class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMRoi.hpp
 * @brief COMPTEL region of interest class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMROI_HPP
#define GCOMROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"
#include "GCOMInstDir.hpp"


/***********************************************************************//**
 * @class GCOMRoi
 *
 * @brief COMPTEL region of interest class
 *
 * The COMPTEL region of interest class defines the event direction
 * region that is used for unbinned data analysis.
 ***************************************************************************/
class GCOMRoi : public GRoi {

public:
    // Constructors and destructors
    GCOMRoi(void);
    GCOMRoi(const GCOMRoi& roi);
    GCOMRoi(const GCOMInstDir& centre, const double& radius,
            const double& phibar_min, const double& phibar_max);
    virtual ~GCOMRoi(void);

    // Operators
    GCOMRoi& operator=(const GCOMRoi& roi);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMRoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GCOMInstDir& centre(void) const;
    void               centre(const GCOMInstDir& centre);
    const double&      radius(void) const;
    void               radius(const double& radius);
    const double&      phibar_min(void) const;
    void               phibar_min(const double& phibar_min);
    const double&      phibar_max(void) const;
    void               phibar_max(const double& phibar_max);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMRoi& roi);
    void free_members(void);
    
    // Protected members
    GCOMInstDir m_centre;     //!< Centre of RoI in instrument coordinates
    double      m_radius;     //!< Radius of region of interest
    double      m_phibar_min; //!< Minimum Phibar of region of interest
    double      m_phibar_max; //!< Minimum Phibar of region of interest
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMRoi").
 ***************************************************************************/
inline
std::string GCOMRoi::classname(void) const
{
    return ("GCOMRoi");
}


/***********************************************************************//**
 * @brief Return centre of region of interest
 *
 * @return Centre of region of interest centre.
 *
 * Returns the instrument direction of the centre of the region of interest.
 ***************************************************************************/
inline
const GCOMInstDir& GCOMRoi::centre(void) const
{
    return (m_centre);
}


/***********************************************************************//**
 * @brief Set centre of region of interest
 *
 * @param[in] centre Instrument direction.
 *
 * Set the instrument direction of the centre of the region of interest.
 ***************************************************************************/
inline
void GCOMRoi::centre(const GCOMInstDir& centre)
{
    m_centre = centre;
    return;
}


/***********************************************************************//**
 * @brief Return radius of region of interest
 *
 * @return Radius of region of interest (deg).
 *
 * Returns the radius of the region of interest.
 ***************************************************************************/
inline
const double& GCOMRoi::radius(void) const
{
    return (m_radius);
}


/***********************************************************************//**
 * @brief Set radius of region of interest
 *
 * @param[in] radius Radius of region of interest (deg).
 *
 * Set the radius of the region of interest.
 ***************************************************************************/
inline
void GCOMRoi::radius(const double& radius)
{
    m_radius = radius;
    return;
}


/***********************************************************************//**
 * @brief Return minimum Phibar of region of interest
 *
 * @return Minimum Phibar of region of interest (deg).
 *
 * Returns the minimum Phibar of region of interest.
 ***************************************************************************/
inline
const double& GCOMRoi::phibar_min(void) const
{
    return (m_phibar_min);
}


/***********************************************************************//**
 * @brief Set minimum Phibar of region of interest
 *
 * @param[in] phibar_min Minimum Phibar of region of interest (deg).
 *
 * Set the minimum Phibar of region of interest.
 ***************************************************************************/
inline
void GCOMRoi::phibar_min(const double& phibar_min)
{
    m_phibar_min = phibar_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum Phibar of region of interest
 *
 * @return Maximum Phibar of region of interest (deg).
 *
 * Returns the maximum Phibar of region of interest.
 ***************************************************************************/
inline
const double& GCOMRoi::phibar_max(void) const
{
    return (m_phibar_max);
}


/***********************************************************************//**
 * @brief Set maximum Phibar of region of interest
 *
 * @param[in] phibar_min Maximum Phibar of region of interest (deg).
 *
 * Set the maximum Phibar of region of interest.
 ***************************************************************************/
inline
void GCOMRoi::phibar_max(const double& phibar_max)
{
    m_phibar_max = phibar_max;
    return;
}

#endif /* GCOMROI_HPP */
