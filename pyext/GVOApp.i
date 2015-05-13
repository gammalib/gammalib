/***************************************************************************
 *                       GVOApp.i - VO SAMP Hub class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Thierry Louge                                    *
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
 * @file GVOApp.i
 * @brief SAMP hub class interface definition
 * @author Thierry Louge
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GVOApp.hpp"
%}


/***********************************************************************//**
 * @class GVOApp
 *
 * @brief VO Hub class
 ***************************************************************************/
class GVOApp : public GBase {

public:
    // Constructors and destructors
    GVOApp(void);
    GVOApp(const GVOApp& hub);
    virtual ~GVOApp(void);

    // Operators
    GVOApp& operator=(const GVOApp& hub);

    // Methods
    void    clear(void);
    GVOApp* clone(void) const;
    void    start(void);
};


/***********************************************************************//**
 * @brief GVOApp class extension
 ***************************************************************************/
%extend GVOApp {
    GVOApp copy() {
        return (*self);
    }
};


