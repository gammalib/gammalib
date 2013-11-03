/***************************************************************************
 *                      GSkyRegionCircle.i - Sky region class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Michael Mayer                                    *
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
 * @file GSkyRegionCircle.i
 * @brief Circular sky region class interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegionCircle.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSkyRegionCircle
 *
 * @brief GSkyRegionCircle dervied sky region class
 ***************************************************************************/
class GSkyRegionCircle : public GSkyRegion {
public:
    // Constructors and destructors
    GSkyRegionCircle(void);
    GSkyRegionCircle(GSkyDir& centre, const double& radius);
    GSkyRegionCircle(const double& ra, const double& dec, const double& radius);
    explicit GSkyRegionCircle(const std::string& line);
    GSkyRegionCircle(const GSkyRegionCircle& region);
    virtual ~GSkyRegionCircle(void);

     // Implemented methods
    void              clear(void);
    GSkyRegionCircle* clone(void) const;
    double            radius(void) const;
    void              radius(const double& radius);
    GSkyDir           centre(void) const;
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
 * @brief GSkyRegionCircle class extension
 ***************************************************************************/
%extend GSkyRegionCircle {
    GSkyRegionCircle copy() {
        return (*self);
    }
};
