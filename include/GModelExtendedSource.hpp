/***************************************************************************
 *        GModelExtendedSource.hpp  -  Extended source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GModelExtendedSource.hpp
 * @brief Extended source model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELEXTENDEDSOURCE_HPP
#define GMODELEXTENDEDSOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSky.hpp"
#include "GModelRadial.hpp"
#include "GModelSpectral.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelExtendedSource
 *
 * @brief Point source model class interface defintion
 *
 * This class implements a factorised point source model.
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

    // Operators
    virtual GModelExtendedSource& operator= (const GModelExtendedSource& model);

    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GModelExtendedSource* clone(void) const;
    virtual std::string           type(void) const { return "ExtendedSource"; }
    virtual std::string           print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelExtendedSource& model);
    void free_members(void);
};

#endif /* GMODELEXTENDEDSOURCE_HPP */
