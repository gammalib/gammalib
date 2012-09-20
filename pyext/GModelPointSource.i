/***************************************************************************
 *       GModelPointSource.i  -  Point source model class python I/F       *
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
 * @file GModelPointSource.i
 * @brief GModelPointSource class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelPointSource.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelPointSource
 *
 * @brief Point source model class python interface defintion
 ***************************************************************************/
class GModelPointSource : public GModelSky {

public:
    // Constructors and destructors
    GModelPointSource(void);
    explicit GModelPointSource(const GXmlElement& xml);
    explicit GModelPointSource(const GModelSpatialPtsrc& ptsrc, const GModelSpectral& spectral);
    explicit GModelPointSource(const GXmlElement& ptsrc, const GXmlElement& spectral);
    GModelPointSource(const GModelPointSource& model);
    virtual ~GModelPointSource(void);

    // Implemented pure virtual methods
    void               clear(void);
    GModelPointSource* clone(void) const;
    std::string        type(void) const;
};


/***********************************************************************//**
 * @brief GModelPointSource class extension
 ***************************************************************************/
%extend GModelPointSource {
    GModelPointSource copy() {
        return (*self);
    }
};
