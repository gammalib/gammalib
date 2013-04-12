/***************************************************************************
 *         GModelSpatialRegistry.i - Spatial model registry class          *
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
 * @file GModelSpatialRegistry.i
 * @brief Spatial model registry class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialRegistry
 *
 * @brief Interface definition for the spatial model registry class
 ***************************************************************************/
class GModelSpatialRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelSpatialRegistry(void);
    GModelSpatialRegistry(const GModelSpatial* model);
    GModelSpatialRegistry(const GModelSpatialRegistry& registry);
    virtual ~GModelSpatialRegistry(void);

    // Methods
    int            size(void) const { return m_number; }
    GModelSpatial* alloc(const std::string& name) const;
    std::string    name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelSpatialRegistry class extension
 ***************************************************************************/
%extend GModelSpatialRegistry {
    char *__str__() {
        return gammalib::tochar(self->print());
    }
};
