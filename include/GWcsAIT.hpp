/***************************************************************************
 *                GWcsAIT.hpp - Aitoff (AIT) projection class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Juergen Knoedlseder                         *
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
 * @file GWcsAIT.hpp
 * @brief Aitoff (AIT) projection class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCSAIT_HPP
#define GWCSAIT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GWcs.hpp"


/***********************************************************************//**
 * @class GWcsAIT
 *
 * @brief Aitoff (AIT) projection class definition
 *
 * This class implements the Aitoff projection for the World Coordinate
 * System.
 ***************************************************************************/
class GWcsAIT : public GWcs {

public:
    // Constructors and destructors
    GWcsAIT(void);
    explicit GWcsAIT(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2);
    GWcsAIT(const GWcsAIT& wcs);
    virtual ~GWcsAIT(void);

    // Operators
    GWcsAIT& operator=(const GWcsAIT& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GWcsAIT*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GWcsAIT& wcs);
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
 * @return String containing the class name ("GWcsAIT").
 ***************************************************************************/
inline
std::string GWcsAIT::classname(void) const
{
    return ("GWcsAIT");
}


/***********************************************************************//**
 * @brief Return projection code
 *
 * @return Projection code.
 *
 * Returns the projection code "AIT".
 ***************************************************************************/
inline
std::string GWcsAIT::code(void) const
{
    return "AIT";
}


/***********************************************************************//**
 * @brief Return projection name
 *
 * @return Projection name.
 *
 * Returns the projection name.
 ***************************************************************************/
inline
std::string GWcsAIT::name(void) const
{
    return "Aitoff";
}

#endif /* GWCSAIT_HPP */
