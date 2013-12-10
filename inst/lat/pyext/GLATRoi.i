/***************************************************************************
 *              GLATRoi.i - Fermi/LAT region of interest class             *
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
 * @file GLATRoi.i
 * @brief Fermi/LAT region of interest class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATRoi.hpp"
%}


/***********************************************************************//**
 * @class GLATRoi
 *
 * @brief Fermi/LAT region of interest class
 ***************************************************************************/
class GLATRoi : public GRoi {
public:
    // Constructors and destructors
    GLATRoi(void);
    GLATRoi(const GLATInstDir& centre, const double& radius);
    GLATRoi(const GLATRoi& roi);
    virtual ~GLATRoi(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GLATRoi*    clone(void) const;
    virtual bool        contains(const GEvent& event) const;

    // Other methods
    const GLATInstDir&  centre(void) const;
    const double&       radius(void) const;
    void                centre(const GLATInstDir& centre);
    void                radius(const double& radius);
};


/***********************************************************************//**
 * @brief GLATRoi class extension
 ***************************************************************************/
%extend GLATRoi {
    GLATRoi copy() {
        return (*self);
    }
};
