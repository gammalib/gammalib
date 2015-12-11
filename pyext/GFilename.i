/***************************************************************************
 *                       GFilename.i - Filename class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Juergen Knoedlseder                              *
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
 * @file GFilename.i
 * @brief Filename class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFilename.hpp"
%}


/***********************************************************************//**
 * @class GFilename
 *
 * @brief Filename class
 ***************************************************************************/
class GFilename : public GBase {

public:
    // Constructors and destructors
    GFilename(void);
    GFilename(const std::string& filename);
    GFilename(const GFilename& filename);
    GFilename(const char* filename);
    virtual ~GFilename(void);

    // Methods
    void               clear(void);
    GFilename*         clone(void) const;
    std::string        classname(void) const;
    bool               empty(void) const;
    int                size(void) const;
    int                length(void) const;
    const std::string& filename(void) const;
    std::string        extname(const std::string& defaultname = "") const;
    const int&         extno(void) const;
    const int&         extver(void) const;
    bool               has_extname(void) const;
    bool               has_extno(void) const;
    bool               has_extver(void) const;
};


/***********************************************************************//**
 * @brief GFilename class extension
 ***************************************************************************/
%extend GFilename {
    GFilename copy() {
        return (*self);
    }
};
