/***************************************************************************
 *               GWcsTAN.hpp - Gnomonic (TAN) projection class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GWcsTAN.hpp
 * @brief Gnomonic (TAN) projection class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCSTAN_HPP
#define GWCSTAN_HPP

/* __ Includes ___________________________________________________________ */
#include "GWcs.hpp"


/***********************************************************************//**
 * @class GWcsTAN
 *
 * @brief Gnomonic (TAN) projection class definition
 *
 * This class implements the "gnomonic" projection for the World Coordinate
 * System.
 ***************************************************************************/
class GWcsTAN : public GWcs {

public:
    // Constructors and destructors
    GWcsTAN(void);
    explicit GWcsTAN(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2);
    GWcsTAN(const GWcsTAN& wcs);
    virtual ~GWcsTAN(void);

    // Operators
    GWcsTAN& operator= (const GWcsTAN& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GWcsTAN*    clone(void) const;
    virtual std::string code(void) const { return "TAN"; }
    virtual std::string name(void) const { return "gnomonic"; }
    virtual std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GWcsTAN& wcs);
    void free_members(void);
    void prj_set(void) const;
    void prj_x2s(int nx, int ny, int sxy, int spt, 
                 const double* x, const double* y,
                 double* phi, double* theta, int* stat) const;
    void prj_s2x(int nphi, int ntheta, int spt, int sxy,
                 const double* phi, const double* theta,
                 double* x, double* y, int* stat) const;
};

#endif /* GWCSTAN_HPP */
