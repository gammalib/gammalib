/***************************************************************************
 *         GCTAEdisp.hpp - Abstract CTA energy dispersion base class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAEdisp.hpp
 * @brief Abstract CTA energy dispersion base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAEDISP_HPP
#define GCTAEDISP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GEbounds.hpp"
#include "GFilename.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GEnergy;


/***********************************************************************//**
 * @class GCTAEdisp
 *
 * @brief Abstract base class for the CTA energy dispersion
 *
 * This class implements the abstract base class for the CTA energy
 * dispersion.
 ***************************************************************************/
class GCTAEdisp : public GBase {

public:
    // Constructors and destructors
    GCTAEdisp(void);
    GCTAEdisp(const GCTAEdisp& edisp);
    virtual ~GCTAEdisp(void);

    // Pure virtual operators
    virtual double operator()(const double& logEobs, 
                              const double& logEsrc, 
                              const double& theta = 0.0, 
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0) const = 0;

    // Operators
    GCTAEdisp& operator=(const GCTAEdisp& edisp);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GCTAEdisp*  clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual void        load(const GFilename& filename) = 0;
    virtual GFilename   filename(void) const = 0;
    virtual GEnergy     mc(GRan&         ran,
                           const double& logE,
                           const double& theta = 0.0,
                           const double& phi = 0.0,
                           const double& zenith = 0.0,
                           const double& azimuth = 0.0) const = 0;
    virtual GEbounds    ebounds_obs(const double& logEsrc,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0) const = 0;
    virtual GEbounds    ebounds_src(const double& logEobs,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0) const = 0;
    virtual double      prob_erecobin(const GEnergy& ereco_min,
                                      const GEnergy& ereco_max,
                                      const GEnergy& etrue,
                                      const double&  theta) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdisp& edisp);
    void free_members(void);
};

#endif /* GCTAEDISP_HPP */
