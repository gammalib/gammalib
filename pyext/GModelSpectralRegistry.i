/***************************************************************************
 *        GModelSpectralRegistry.i - Spectral model registry class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralRegistry.i
 * @brief Spectral model registry class definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralRegistry.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralRegistry
 *
 * @brief Interface definition for the spectral model registry class
 ***************************************************************************/
class GModelSpectralRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelSpectralRegistry(void);
    GModelSpectralRegistry(const GModelSpectral* model);
    GModelSpectralRegistry(const GModelSpectralRegistry& registry);
    virtual ~GModelSpectralRegistry(void);

    // Methods
    std::string     classname(void) const;
    int             size(void) const;
    GModelSpectral* alloc(const GXmlElement& xml) const;
    std::string     name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelSpectralRegistry class extension
 ***************************************************************************/
%extend GModelSpectralRegistry {
};
