/***************************************************************************
 *          GWcsGLS.hpp - Global Sinusoidal (GLS) projection class         *
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
 * @file GWcsGLS.hpp
 * @brief Global Sinusoidal (GLS) projection class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCSGLS_HPP
#define GWCSGLS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GWcsSFL.hpp"


/***********************************************************************//**
 * @class GWcsGLS
 *
 * @brief Global Sinusoidal (GLS) projection class definition
 *
 * This class implements the Global Sinusoidal projection for the World
 * Coordinate System. The Global Sinusoidal projection is an alias of the
 * Sanson-Flamsteed projection, and the GWcsGLS therefore derives from the
 * GWcsSFL class.
 ***************************************************************************/
class GWcsGLS : public GWcsSFL {

public:
    // Constructors and destructors
    GWcsGLS(void);
    GWcsGLS(const std::string& coords,
            const double& crval1, const double& crval2,
            const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2);
    GWcsGLS(const GWcsGLS& wcs);
    virtual ~GWcsGLS(void);

    // Operators
    GWcsGLS& operator=(const GWcsGLS& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GWcsGLS*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Private methods
    void init_members(void);
    void copy_members(const GWcsGLS& wcs);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GWcsGLS").
 ***************************************************************************/
inline
std::string GWcsGLS::classname(void) const
{
    return ("GWcsGLS");
}


/***********************************************************************//**
 * @brief Return projection code
 *
 * @return Projection code.
 *
 * Returns the projection code "GLS".
 ***************************************************************************/
inline
std::string GWcsGLS::code(void) const
{
    return "GLS";
}


/***********************************************************************//**
 * @brief Return projection name
 *
 * @return Projection name.
 *
 * Returns the projection name.
 ***************************************************************************/
inline
std::string GWcsGLS::name(void) const
{
    return "Global Sinusoidal";
}

#endif /* GWCSGLS_HPP */
