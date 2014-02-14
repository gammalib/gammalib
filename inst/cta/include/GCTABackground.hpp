/***************************************************************************
 *            GCTABackground.hpp - CTA background model base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackground.hpp
 * @brief CTA background model base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTABACKGROUND_HPP
#define GCTABACKGROUND_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GModelSpectralNodes.hpp"
#include "GCTAInstDir.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GCTABackground
 *
 * @brief Abstract base class for the CTA background model
 ***************************************************************************/
class GCTABackground : public GBase {

public:
    // Constructors and destructors
    GCTABackground(void);
    GCTABackground(const GCTABackground& bgd);
    virtual ~GCTABackground(void);

    // Pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety,
                              const bool&   etrue = false) const = 0;

    // Operators
    GCTABackground& operator=(const GCTABackground& bgd);

    // Pure virtual methods
    virtual void                       clear(void) = 0;
    virtual GCTABackground*            clone(void) const = 0;
    virtual void                       load(const std::string& filename) = 0;
    virtual std::string                filename(void) const = 0;
    virtual GCTAInstDir                mc(const GEnergy& energy,
                                          const GTime& time,
                                          GRan& ran) const = 0;
    virtual const GModelSpectralNodes& spectrum(void) const = 0;
    virtual std::string                print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTABackground& bgd);
    void free_members(void);
};

#endif /* GCTABACKGROUND_HPP */
