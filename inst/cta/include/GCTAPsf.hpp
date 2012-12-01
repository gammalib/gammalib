/***************************************************************************
 *            GCTAPsf.hpp - CTA point spread function base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCTAPsf.hpp
 * @brief CTA point spread function base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAPSF_HPP
#define GCTAPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GCTAPsf
 *
 * @brief Abstract base class for the CTA point spread function
 *
 * This class implements the abstract base class for the CTA point spread
 * function.
 ***************************************************************************/
class GCTAPsf : public GBase {

public:
    // Constructors and destructors
    GCTAPsf(void);
    GCTAPsf(const GCTAPsf& psf);
    virtual ~GCTAPsf(void);

    // Pure virtual operators
    virtual double operator()(const double& delta,
                              const double& logE, 
                              const double& theta = 0.0, 
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0,
                              const bool&   etrue = true) const = 0;

    // Operators
    GCTAPsf& operator=(const GCTAPsf& psf);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GCTAPsf*    clone(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual std::string filename(void) const = 0;
    virtual double      mc(GRan&         ran,
                           const double& logE, 
                           const double& theta = 0.0, 
                           const double& phi = 0.0,
                           const double& zenith = 0.0,
                           const double& azimuth = 0.0,
                           const bool&   etrue = true) const = 0;
    virtual double      delta_max(const double& logE, 
                                  const double& theta = 0.0, 
                                  const double& phi = 0.0,
                                  const double& zenith = 0.0,
                                  const double& azimuth = 0.0,
                                  const bool&   etrue = true) const = 0;
    virtual std::string print(void) const = 0;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAPsf& psf);
    void free_members(void);
};

#endif /* GCTAPSF_HPP */
