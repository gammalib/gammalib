/***************************************************************************
 *    GModelSpatialDiffuse.i - Abstract diffuse spatial model base class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialDiffuse.i
 * @brief Abstract diffuse spatial model base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialDiffuse.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialDiffuse
 *
 * @brief Abstract diffuse spatial model base class
 ***************************************************************************/
class GModelSpatialDiffuse : public GModelSpatial {
public:
    // Constructors and destructors
    GModelSpatialDiffuse(void);
    GModelSpatialDiffuse(const GModelSpatialDiffuse& model);
    virtual ~GModelSpatialDiffuse(void);

    // Pure virtual methods
    virtual void                  clear(void) = 0;
    virtual GModelSpatialDiffuse* clone(void) const = 0;
    virtual std::string           type(void) const = 0;
    virtual double                eval(const GSkyDir& srcDir) const = 0;
    virtual double                eval_gradients(const GSkyDir& srcDir) const = 0;
    virtual GSkyDir               mc(GRan& ran) const = 0;
    virtual void                  read(const GXmlElement& xml) = 0;
    virtual void                  write(GXmlElement& xml) const = 0;
};


/***********************************************************************//**
 * @brief GModelSpatialDiffuse class extension
 ***************************************************************************/
%extend GModelSpatialDiffuse {
};
