/***************************************************************************
 *          GModelDiffuseSource.hpp  -  Diffuse source model class         *
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
 * @file GModelDiffuseSource.hpp
 * @brief GModelDiffuseSource class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELDIFFUSESOURCE_HPP
#define GMODELDIFFUSESOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSky.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelDiffuseSource
 *
 * @brief Diffuse source model class interface defintion
 *
 * This class implements a factorised diffuse source model.
 ***************************************************************************/
class GModelDiffuseSource : public GModelSky {

public:
    // Constructors and destructors
    GModelDiffuseSource(void);
    explicit GModelDiffuseSource(const GXmlElement& xml);
    explicit GModelDiffuseSource(const GModelSpatial& spatial, const GModelSpectral& spectral);
    explicit GModelDiffuseSource(const GXmlElement& spatial, const GXmlElement& spectral);
    GModelDiffuseSource(const GModelDiffuseSource& model);
    virtual ~GModelDiffuseSource(void);

    // Operators
    GModelDiffuseSource& operator= (const GModelDiffuseSource& model);

    // Implemented pure virtual methods
    void                 clear(void);
    GModelDiffuseSource* clone(void) const;
    std::string          type(void) const { return "DiffuseSource"; }
    std::string          print(void) const;
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelDiffuseSource& model);
    void free_members(void);
};

#endif /* GMODELDIFFUSESOURCE_HPP */
