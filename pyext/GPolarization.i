/***************************************************************************
 *                    GPolarization.i - Polarization class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GPolarization.i
 * @brief GPolarization class python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPolarization.hpp"
%}


/***********************************************************************//**
 * @class GPolarization
 *
 * @brief Class that handles polarizations.
 ***************************************************************************/
class GPolarization : public GBase {
public:
    // Constructors and destructors
    GPolarization(void);
    GPolarization(const GPolarization& polarization);
    virtual ~GPolarization(void);
 
    // Operators

    // Methods
    void           clear(void);
    GPolarization* clone(void) const;
    std::string    classname(void) const;
};


/***********************************************************************//**
 * @brief GPolarization class extension
 ***************************************************************************/
%extend GPolarization {
    GPolarization copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = {}
        return state
    def __setstate__(self, state):
        self.__init__()
}
};
