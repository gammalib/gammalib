/***************************************************************************
 *                 GCTAAeff.hpp - CTA effective area base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAAeff.hpp
 * @brief CTA effective area base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAAEFF_HPP
#define GCTAAEFF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"


/***********************************************************************//**
 * @class GCTAAeff
 *
 * @brief Abstract base class for the CTA effective area
 *
 * This class implements the abstract base class for the CTA effective area.
 * The effective area is accessed using the following arguments
 *    logE - log10 of true (or measured) energy in TeV
 *    theta - Offset angle of true photon direction in camera system (rad)
 *    phi - Azimuth angle of true photon direction in camera system (rad)
 *    zenith - Zenith angle of pointing in Earth system (rad)
 *    azimuth - Azimuth angle of pointing in Earth system (rad)
 *    etrue - Use true or measured energy
 ***************************************************************************/
class GCTAAeff : public GBase {

public:
    // Constructors and destructors
    GCTAAeff(void);
    GCTAAeff(const GCTAAeff& aeff);
    virtual ~GCTAAeff(void);

    // Pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& theta = 0.0, 
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0,
                              const bool&   etrue = true) const = 0;

    // Operators
    GCTAAeff& operator=(const GCTAAeff& aeff);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GCTAAeff*   clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual std::string filename(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAAeff& aeff);
    void free_members(void);
};

#endif /* GCTAAEFF_HPP */
