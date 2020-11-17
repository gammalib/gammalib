/***************************************************************************
 *         GCTAModelSpatialRegistry.i - Spatial model registry class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAModelSpatialRegistry.i
 * @brief Spatial model registry class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelSpatialRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelSpatialRegistry
 *
 * @brief Interface definition for the spatial model registry class
 ***************************************************************************/
class GCTAModelSpatialRegistry : public GRegistry {

public:
    // Constructors and destructors
    GCTAModelSpatialRegistry(void);
    GCTAModelSpatialRegistry(const GCTAModelSpatial* model);
    GCTAModelSpatialRegistry(const GCTAModelSpatialRegistry& registry);
    virtual ~GCTAModelSpatialRegistry(void);

    // Methods
    std::string       classname(void) const;
    int               size(void) const;
    GCTAModelSpatial* alloc(const GXmlElement& xml) const;
    std::string       name(const int& index) const;
};


/***********************************************************************//**
 * @brief GCTAModelSpatialRegistry class extension
 ***************************************************************************/
%extend GCTAModelSpatialRegistry {
};
