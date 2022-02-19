/***************************************************************************
 *                         GDaemon.i - Daemon class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GDaemon.i
 * @brief Daemon class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GDaemon.hpp"
%}


/***********************************************************************//**
 * @class GDaemon
 *
 * @brief Daemon class
 ***************************************************************************/
class GDaemon : public GBase {

public:
    // Constructors and destructors
    GDaemon(void);
    //GDaemon(const GDaemon& daemon);
    virtual ~GDaemon(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GDaemon*    clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    void start(void);
    bool alive(void) const;
};


/***********************************************************************//**
 * @brief GDaemon class extension
 ***************************************************************************/
%extend GDaemon {
    GDaemon copy() {
        return (*self);
    }
};
