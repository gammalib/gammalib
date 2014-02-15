/***************************************************************************
 *             GCTAEdispRmf.hpp - CTA RMF energy dispersion class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Christoph Deil & Ellis Owen                      *
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
 * @file GCTAEdispRmf.hpp
 * @brief CTA RMF energy dispersion class definition
 * @author Christoph Deil & Ellis Owen
 */

#ifndef GCTAEDISPRMF_HPP
#define GCTAEDISPRMF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRan.hpp"
#include "GRmf.hpp"
#include "GVector.hpp"
#include "GCTAEdisp.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GEnergy;
class GEbounds;


/***********************************************************************//**
 * @class GCTAEdispRmf
 *
 * @brief CTA Redistribution Matrix File (RMF) energy dispersion class
 ***************************************************************************/
class GCTAEdispRmf : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispRmf(void);
    explicit GCTAEdispRmf(const std::string& filename);
    GCTAEdispRmf(const GCTAEdispRmf& edisp);
    virtual ~GCTAEdispRmf(void);

    // Operators
    GCTAEdispRmf& operator=(const GCTAEdispRmf& edisp);
    double operator()(const double& logEobs,
                      const double& logEsrc,
                      const double& theta = 0.0,
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;

    // Implemented pure virtual methods
    void          clear(void);
    GCTAEdispRmf* clone(void) const;
    void          load(const std::string& filename);
    std::string   filename(void) const;
    GEnergy       mc(GRan& ran,
                     const double& logEsrc,
                     const double& theta = 0.0,
                     const double& phi = 0.0,
                     const double& zenith = 0.0,
                     const double& azimuth = 0.0) const;
    GEbounds      ebounds_obs(const double& logEsrc,
                              const double& theta = 0.0,
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0) const;
    GEbounds      ebounds_src(const double& logEobs,
                              const double& theta = 0.0,
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int   size(void) const;
    const GRmf& rmf(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdispRmf& psf);
    void free_members(void);
    void set_mc_cache(void);

    // Members
    std::string m_filename;  //!< Name of response file
    GRmf        m_rmf;       //!< Redistribution matrix file

    // Monte Carlo cache
    mutable std::vector<int>     m_mc_measured_start;
    mutable std::vector<GVector> m_mc_measured_cdf;
};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the Redistribution Matrix was loaded.
 ***************************************************************************/
inline
std::string GCTAEdispRmf::filename(void) const
{
    return m_filename;
}


/***********************************************************************//**
 * @brief Return Redistribution Matrix File
 *
 * @return Returns Redistribution Matrix File.
 ***************************************************************************/
inline
const GRmf& GCTAEdispRmf::rmf(void) const
{
    return m_rmf;
}

#endif /* GCTAEDISPRMF_HPP */
