/***************************************************************************
 *                      GSkyRegionCircle.i - Sky region class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Michael Mayer                               *
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
 * @brief Sky region class SWIG file.
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
	GSkyRegionCircle(const GSkyRegionCircle& circle);
    explicit GSkyRegionCircle(const std::string &line);
    GSkyRegionCircle(GSkyDir &centre, const double &radius);
    // GSkyRegionCircle(const double &ra, const double &dec, const double &radius);
    virtual ~GSkyRegionCircle(void);

     // Implemented methods
    void         				clear(void);
    GSkyRegionCircle*  			clone(void) const;
    GSkyDir 					centre(void) const;
    double 						radius(void) const;
    void 						radius(const double& radius);
    void 						centre(const GSkyDir& centre);
    void 						centre(const double& ra,const double& dec);
    void         				read(const std::string line) const;
    std::string  				write() const;
    bool         				contains(const GSkyDir& dir) const;
	bool         				contains(const GSkyRegion& reg) const;
	bool         				overlaps(const GSkyRegion& reg) const;

};
