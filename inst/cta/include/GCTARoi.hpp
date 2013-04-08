/***************************************************************************
 *                GCTARoi.hpp - CTA region of interest class               *
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
 * @file GCTARoi.hpp
 * @brief GCTARoi class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GCTAROI_HPP
#define GCTAROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"
#include "GCTAInstDir.hpp"


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
    GCTARoi(const GCTARoi& roi);
    virtual ~GCTARoi(void);

    // Operators
    GCTARoi& operator= (const GCTARoi& roi);

    // Implemented pure virtual base class methods
    void         clear(void);
    GCTARoi*     clone(void) const;
    std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GCTAInstDir centre(void) const { return m_centre; }
    double      radius(void) const { return m_radius; }
    void        centre(const GCTAInstDir& centre) { m_centre=centre; return; }
    void        radius(const double& radius) { m_radius=radius; return; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTARoi& roi);
    void free_members(void);
    
    // Protected members
    GCTAInstDir m_centre;   //!< Centre of ROI in instrument coordinates
    double      m_radius;   //!< Radius of ROI in degrees
};

#endif /* GCTAROI_HPP */
