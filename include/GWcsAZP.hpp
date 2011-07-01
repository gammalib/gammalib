/***************************************************************************
 *  GWcsAZP.hpp  -  Zenithal/azimuthal perspective (AZP) projection class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GWcsAZP.hpp
 * @brief Zenithal/azimuthal perspective (AZP) projection class definition
 * @author J. Knodlseder
 */

#ifndef GWCSAZP_HPP
#define GWCSAZP_HPP

/* __ Includes ___________________________________________________________ */
#include "GWcslib.hpp"


/***********************************************************************//**
 * @class GWcsAZP
 *
 * @brief Zenithal/azimuthal perspective (AZP) projection class definition
 *
 * This class implements the "zenithal/azimuthal perspective" projection for
 * the World Coordinate System.
 ***************************************************************************/
class GWcsAZP : public GWcslib {

public:
    // Constructors and destructors
    GWcsAZP(void);
    explicit GWcsAZP(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2);
    GWcsAZP(const GWcsAZP& wcs);
    virtual ~GWcsAZP(void);

    // Operators
    GWcsAZP& operator= (const GWcsAZP& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GWcsAZP*    clone(void) const;
    virtual std::string code(void) const { return "AZP"; }
    virtual std::string name(void) const { return "zenithal/azimuthal perspective"; }
    virtual std::string print(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GWcsAZP& wcs);
    void free_members(void);
    void prj_set(void) const;
    void prj_x2s(int nx, int ny, int sxy, int spt, 
                 const double* x, const double* y,
                 double* phi, double* theta, int* stat) const;
    void prj_s2x(int nphi, int ntheta, int spt, int sxy,
                 const double* phi, const double* theta,
                 double* x, double* y, int* stat) const;
};

#endif /* GWCSAZP_HPP */
