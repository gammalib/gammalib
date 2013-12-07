/***************************************************************************
 *            GEventList.i - Abstract event atom container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GEventList.i
 * @brief Abstract event atom container class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEventList.hpp"
%}


/***********************************************************************//**
 * @class GEventList
 *
 * @brief Abstract event atom container class interface
 ***************************************************************************/
class GEventList : public GEvents {
public:
    // Constructors and destructors
    GEventList(void);
    GEventList(const GEventList& list);
    virtual ~GEventList(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventList* clone(void) const = 0;
    virtual int         size(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual void        save(const std::string& filename,
                             const bool& clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;
    virtual void        roi(const GRoi& roi) = 0;
    virtual const GRoi& roi(void) const = 0;
};


/***********************************************************************//**
 * @brief GEventList class extension
 ***************************************************************************/
%extend GEventList {
};
