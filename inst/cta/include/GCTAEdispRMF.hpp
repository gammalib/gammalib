/***************************************************************************
 *  GCTAEdispRMF.hpp - CTA RMF energy dispersion class *
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
 * @file GCTAEdispRMF.hpp
 * @brief CTA RMF energy dispersion class definition
 * @author Christoph Deil & Ellis Owen
 */

#ifndef GCTAEDISPRMF_HPP
#define GCTAEDISPRMF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GRan.hpp"
#include "GMatrixSparse.hpp"
#include "GEbounds.hpp"
#include "GCTAEdisp.hpp"

/* __ Forward declarations _______________________________________________ */
class GCTAResponse;


/***********************************************************************//**
 * @class GCTAEdispRMF
 *
 * @brief CTA performance table energy dispersion class
 *
 * This class implements the CTA energy dispersion response as function
 * of energy as determined from a performance table. The performance table is
 * an ASCII file that specifies the CTA performance parameters in a simple
 * way.
 ***************************************************************************/
class GCTAEdispRMF : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispRMF(void);
    explicit GCTAEdispRMF(const std::string& filename);
    GCTAEdispRMF(const GCTAEdispRMF& edisp);
    virtual ~GCTAEdispRMF(void);

    // Operators
    GCTAEdispRMF& operator=(const GCTAEdispRMF& edisp);
    double operator()(const double& logEobs,
                      const double& logEsrc,
                      const double& theta = 0.0,
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;

    // Implemented pure virtual methods
    void                clear(void);
    GCTAEdispRMF* 		clone(void) const;
    void                load(const std::string& filename);
    std::string         filename(void) const;
    GEnergy             mc(GRan&         ran,
                               const double& logE,
                               const double& theta = 0.0,
                               const double& phi = 0.0,
                               const double& zenith = 0.0,
                               const double& azimuth = 0.0) const;
    GEbounds            ebounds_obs(const double& logEsrc,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0) const;
    GEbounds            ebounds_src(const double& logEobs,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0) const;
    std::string         print(const GChatter& chatter = NORMAL) const;

    // Methods
    int           size(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdispRMF& psf);
    void free_members(void);
    void convert_cdf(void);

    // Members
    std::string         m_filename;  //!< Name of response file

    // For now, let's just store the plain RMF info and implement it
    // as a step function model.
    GEbounds            m_ebounds_src;     //!< Source energy boundaries
    GEbounds            m_ebounds_obs;     //!< Observed energy boundaries
    GMatrixSparse       m_matrix;          //!< Sparse redistribution matrix
    GVector 			m_mc_cache;        //!< Monte Carlo cache
    std::vector<GVector> m_cdf_cache;     //!< Vector of GVectors from RMF file
    // TODO: optionally implement interpolated (linear or spline) model
    // and store interpolation object here.
    //GNodeArray          m_logE;      //!< log(E) nodes for interpolation
    //std::vector<double> m_aeff;      //!< Effective area in cm2

};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the energy resolution was loaded
 ***************************************************************************/
inline
std::string GCTAEdispRMF::filename(void) const
{
    return m_filename;
}

#endif /* GCTAEDISPRMF_HPP */
