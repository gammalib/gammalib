/***************************************************************************
 *                       GVOHub.i - VO SAMP Hub class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Thierry Louge                                    *
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
 * @file GVOHub.i
 * @brief SAMP hub class interface definition
 * @author Thierry Louge
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GVOHub.hpp"
%}


/***********************************************************************//**
 * @class GVOHub
 *
 * @brief VO Hub class
 ***************************************************************************/
class GVOHub {
public:
    // Constructors and destructors
    // Constructors and destructors
    GVOHub(void);
    GVOHub(const GVOHub& hub);
    virtual ~GVOHub(void);

    // Methods
    void        clear(void);
    GVOHub*  clone(void) const;
    //void        connect(void);
    //void        disconnect(void);
    //bool        has_hub(void) const;
    //bool        is_connected(void) const;
    //GXml        response(void) const;
};


/***********************************************************************//**
 * @brief GVOHub class extension
 ***************************************************************************/
/*
%extend GVOHub {
    GVOHub copy() {
        return (*self);
    }
};
*/

