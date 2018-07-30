/***************************************************************************
 *   GWcsARC.hpp - Zenithal/azimuthal equidistant (ARC) projection class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GWcsARC.hpp
 * @brief Zenithal/azimuthal equidistant (ARC) projection class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCSARC_HPP
#define GWCSARC_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GWcs.hpp"


/***********************************************************************//**
 * @class GWcsARC
 *
 * @brief Zenithal/azimuthal equidistant (ARC) projection class definition
 *
 * This class implements the Zenithal/azimuthal equidistant projection for
 * the World Coordinate System.
 ***************************************************************************/
class GWcsARC : public GWcs {

public:
    // Constructors and destructors
    GWcsARC(void);
    GWcsARC(const std::string& coords,
            const double& crval1, const double& crval2,
            const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2);
    GWcsARC(const GWcsARC& wcs);
    virtual ~GWcsARC(void);

    // Operators
    GWcsARC& operator=(const GWcsARC& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GWcsARC*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GWcsARC& wcs);
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
 * @return String containing the class name ("GWcsARC").
 ***************************************************************************/
inline
std::string GWcsARC::classname(void) const
{
    return ("GWcsARC");
}


/***********************************************************************//**
 * @brief Return projection code
 *
 * @return Projection code.
 *
 * Returns the projection code "ARC".
 ***************************************************************************/
inline
std::string GWcsARC::code(void) const
{
    return "ARC";
}


/***********************************************************************//**
 * @brief Return projection name
 *
 * @return Projection name.
 *
 * Returns the projection name.
 ***************************************************************************/
inline
std::string GWcsARC::name(void) const
{
    return "Zenithal/azimuthal equidistant";
}

#endif /* GWCSARC_HPP */
