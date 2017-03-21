/***************************************************************************
 *          GWcsSFL.hpp - Sanson-Flamsteed (SFL) projection class          *
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
 * @file GWcsSFL.hpp
 * @brief Sanson-Flamsteed (SFL) projection class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCSSFL_HPP
#define GWCSSFL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GWcs.hpp"


/***********************************************************************//**
 * @class GWcsSFL
 *
 * @brief Sanson-Flamsteed (SFL) projection class definition
 *
 * This class implements the Sanson-Flamsteed projection for the World
 * Coordinate System.
 ***************************************************************************/
class GWcsSFL : public GWcs {

public:
    // Constructors and destructors
    GWcsSFL(void);
    GWcsSFL(const std::string& coords,
            const double& crval1, const double& crval2,
            const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2);
    GWcsSFL(const GWcsSFL& wcs);
    virtual ~GWcsSFL(void);

    // Operators
    GWcsSFL& operator=(const GWcsSFL& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GWcsSFL*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GWcsSFL& wcs);
    void free_members(void);
    void prj_set(void) const;
    void prj_x2s(int nx, int ny, int sxy, int spt, 
                 const double* x, const double* y,
                 double* phi, double* theta, int* stat) const;
    void prj_s2x(int nphi, int ntheta, int spt, int sxy,
                 const double* phi, const double* theta,
                 double* x, double* y, int* stat) const;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GWcsSFL").
 ***************************************************************************/
inline
std::string GWcsSFL::classname(void) const
{
    return ("GWcsSFL");
}


/***********************************************************************//**
 * @brief Return projection code
 *
 * @return Projection code.
 *
 * Returns the projection code "SFL".
 ***************************************************************************/
inline
std::string GWcsSFL::code(void) const
{
    return "SFL";
}


/***********************************************************************//**
 * @brief Return projection name
 *
 * @return Projection name.
 *
 * Returns the projection name.
 ***************************************************************************/
inline
std::string GWcsSFL::name(void) const
{
    return "Sanson-Flamsteed";
}

#endif /* GWCSSFL_HPP */
