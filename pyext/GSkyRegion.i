/***************************************************************************
 *                      GSkyRegion.i - Sky region class                    *
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
 * @file GSkyRegion.i
 * @brief Sky region class SWIG file.
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegion.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSkyRegion
 *
 * @brief GSkyRegion base abstract class
 ***************************************************************************/
class GSkyRegion : public GBase {
public:
    // Constructors and destructors
    GSkyRegion(void);
    GSkyRegion(const GSkyRegion& region); 
    virtual ~GSkyRegion(void);
	
    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GSkyRegion* clone(void) const = 0;
    virtual void        read(const std::string& regstring) = 0;
    virtual std::string write(void) const = 0;
	virtual bool        contains(const GSkyDir& dir) const = 0;
    virtual bool        overlaps(const GSkyRegion& reg) const = 0;
	virtual bool        contains(const GSkyRegion& reg) const = 0;
	
    // Implemented methods
    const std::string& type(void) const;
    const std::string& name(void) const;
	void               name(const std::string& name);
    double             solidangle(void) const;
};


/***********************************************************************//**
 * @brief GSkyRegion class extension
 ***************************************************************************/
