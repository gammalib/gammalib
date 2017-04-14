/***************************************************************************
 *                   GPhases.i - Phase intervals class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GPhases.i
 * @brief Phase intervals class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPhases.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPhases
 *
 * @brief Phase intervals class
 ***************************************************************************/
class GPhases : public GContainer {

public:
    // Constructors and destructors
    GPhases(void);
    GPhases(const GPhases& phases);
    GPhases(const double& pmin, const double& pmax);
    virtual ~GPhases(void);

    // Methods
    void        clear(void);
    GPhases*    clone(void) const;
    std::string classname(void) const;
    int         size(void) const;
    bool        is_empty(void) const;
    void        append(const double& pmin, const double& pmax);
    void        remove(const int& index);
    void        reserve(const int& num);
    void        extend(const GPhases& phases);
    double      pmin(const int& index) const;
    double      pmax(const int& index) const;
};


/***********************************************************************//**
 * @brief GPhases class extension
 ***************************************************************************/
%extend GPhases {
    int __len__() {
        return (self->size());
    }
    GPhases copy() {
        return (*self);
    }
};
