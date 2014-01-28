/***************************************************************************
 *                      GSkyRegionRing.i - Sky region class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Maria Krause, Anneli Schulz                      *
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
 * @file GSkyRegionRing.i
 * @brief Ring sky region class interface definition
 * @author Maria Krause, Anneli Schulz
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegionRing.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSkyRegionRing
 *
 * @brief Interface for the ring sky region class
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

    // Implemented methods
    void              clear(void);
    GSkyRegionCircle* clone(void) const;
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
};


/***********************************************************************//**
 * @brief GSkyRegionRing class extension
 ***************************************************************************/
%extend GSkyRegionRing {
    GSkyRegionRing copy() {
        return (*self);
    }
};
