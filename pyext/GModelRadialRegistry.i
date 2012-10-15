/***************************************************************************
 *           GModelRadialRegistry.i - Radial model registry class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @file GModelRadialRegistry.i
 * @brief Radial model registry class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelRadialRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelRadialRegistry
 *
 * @brief Radial model registry class
 ***************************************************************************/
class GModelRadialRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelRadialRegistry(void);
    explicit GModelRadialRegistry(const GModelRadial* model);
    GModelRadialRegistry(const GModelRadialRegistry& registry);
    virtual ~GModelRadialRegistry(void);

    // Methods
    int           size(void) const;
    GModelRadial* alloc(const std::string& name) const;
    std::string   name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelRadialRegistry class extension
 ***************************************************************************/
%extend GModelRadialRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
