/***************************************************************************
 *         GModelExtendedSource.i  -  Extended source model class          *
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
 * @file GModelExtendedSource.i
 * @brief Extended source model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelExtendedSource.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelExtendedSource
 *
 * @brief Extended source model class
 ***************************************************************************/
class GModelExtendedSource : public GModelSky {
public:
    // Constructors and destructors
    GModelExtendedSource(void);
    explicit GModelExtendedSource(const GXmlElement& xml);
    explicit GModelExtendedSource(const GModelRadial& radial, const GModelSpectral& spectral);
    explicit GModelExtendedSource(const GXmlElement& radial, const GXmlElement& spectral);
    GModelExtendedSource(const GModelExtendedSource& model);
    virtual ~GModelExtendedSource(void);

    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GModelExtendedSource* clone(void) const;
    virtual std::string           type(void) const;
    virtual void                  read(const GXmlElement& xml);
    virtual void                  write(GXmlElement& xml) const;

    // Other methods
    GModelRadial* radial(void) const;
    GSkyDir       dir(void) const;
};


/***********************************************************************//**
 * @brief GModelExtendedSource class extension
 ***************************************************************************/
%extend GModelExtendedSource {
    GModelExtendedSource copy() {
        return (*self);
    }
};
