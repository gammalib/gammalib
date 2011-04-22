/***************************************************************************
 *    GObservationRegistry.i  -  Observation registry class python I/F     *
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
 * @file GObservationRegistry.i
 * @brief GObservationRegistry class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GObservationRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GObservationRegistry
 *
 * @brief Interface definition for the observation registry class
 ***************************************************************************/
class GObservationRegistry {
public:
    // Constructors and destructors
    GObservationRegistry(void);
    GObservationRegistry(const GObservation* obs);
    GObservationRegistry(const GObservationRegistry& registry);
    virtual ~GObservationRegistry(void);

    // Methods
    int           size(void) const { return m_number; }
    GObservation* alloc(const std::string& type) const;
    std::string   type(const int& index) const;
};


/***********************************************************************//**
 * @brief GObservationRegistry class extension
 ***************************************************************************/
%extend GObservationRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
