/***************************************************************************
 *           GSPIEventCube.hpp  -  SPI event bin container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GSPIEventCube.hpp
 * @brief SPI event bin container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIEVENTCUBE_HPP
#define GSPIEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GSPIEventBin.hpp"
#include "GSPIInstDir.hpp"


/***********************************************************************//**
 * @class GSPIEventCube
 *
 * @brief SPI event bin container class
 *
 * This class is a container class for SPI event bins.
 ***************************************************************************/
class GSPIEventCube : public GEventCube {

public:
    // Constructors and destructors
    GSPIEventCube(void);
    GSPIEventCube(const GSPIEventCube& cube);
    explicit GSPIEventCube(const std::string& filename);
    virtual ~GSPIEventCube(void);

    // Operators
    virtual GSPIEventCube&      operator=(const GSPIEventCube& cube);
    virtual GSPIEventBin*       operator[](const int& index);
    virtual const GSPIEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GSPIEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(int axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename, bool clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual std::string    print(void) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GSPIEventCube& cube);
    void         free_members(void);
    virtual void set_energies(void);
    virtual void set_times(void);

    // Protected members
};

#endif /* GSPIEVENTCUBE_HPP */
