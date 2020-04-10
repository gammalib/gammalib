/***************************************************************************
 *           GInstDir.i - Abstract instrument direction base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GInstDir.i
 * @brief Abstract instrument direction base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GInstDir.hpp"
%}


/***********************************************************************//**
 * @class GInstDir
 *
 * @brief Abstract interface for the instrument direction of an event.
 *
 * The instrument direction of an event is the equivalent of the sky
 * direction (implemented by GSkyDir) but in the instrument data space.
 * The instrument direction may be any kind of position or direction
 * information encoded in the data space, such as incident event
 * reconstructions for imaging devices or detector numbers etc. for
 * non-imaging devices.
 ***************************************************************************/
class GInstDir : public GBase {

public:
    // Constructors and destructors
    GInstDir(void);
    GInstDir(const GInstDir& dir);
    virtual ~GInstDir(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GInstDir*   clone(void) const = 0;
    virtual u_int64_t   hash(void) const = 0;
    virtual std::string classname(void) const = 0;
};


/***********************************************************************//**
 * @brief GInstDir class extension
 ***************************************************************************/
%extend GInstDir {
};
