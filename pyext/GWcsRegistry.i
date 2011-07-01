/***************************************************************************
 *  GWcsRegistry.i  -  Temporal model registry class python I/F  *
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
 * @file GWcsRegistry.i
 * @brief World Coordinate Projection registry class interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GWcsRegistry
 *
 * @brief Interface definition for the WCS registry class
 ***************************************************************************/
class GWcsRegistry {
public:
    // Constructors and destructors
    GWcsRegistry(void);
    GWcsRegistry(const GWcs* prj);
    GWcsRegistry(const GWcsRegistry& registry);
    virtual ~GWcsRegistry(void);

    // Methods
    int         size(void) const { return m_number; }
    GWcs*       alloc(const std::string& code) const;
    std::string code(const int& index) const;
    std::string name(const int& index) const;
    std::string list(void) const;
};


/***********************************************************************//**
 * @brief GWcsRegistry class extension
 ***************************************************************************/
%extend GWcsRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
