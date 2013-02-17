/***************************************************************************
 *          GModelSpatialCube.i  -  Spatial map cube model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GModelSpatialCube.i
 * @brief Spatial map cube model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialCube.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialCube
 *
 * @brief Spatial map cube model
 ***************************************************************************/
class GModelSpatialCube  : public GModelSpatial {
public:
    // Constructors and destructors
    GModelSpatialCube(void);
    explicit GModelSpatialCube(const GXmlElement& xml);
    GModelSpatialCube(const GModelSpatialCube& model);
    virtual ~GModelSpatialCube(void);

    // Implemented pure virtual methods
    virtual void               clear(void);
    virtual GModelSpatialCube* clone(void) const;
    virtual std::string        type(void) const;
    virtual double             eval(const GSkyDir& srcDir) const;
    virtual double             eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir            mc(GRan& ran) const;
    virtual void               read(const GXmlElement& xml);
    virtual void               write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelSpatialCube class extension
 ***************************************************************************/
%extend GModelSpatialCube {
    GModelSpatialCube copy() {
        return (*self);
    }
};
