/***************************************************************************
 *  GCTAEdispPerfTable.hpp - CTA performance table energy dispersion class *
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
 * @file GCTAEdispPerfTable.hpp
 * @brief CTA performance table energy dispersion class definition
 * @author Christoph Deil & Ellis Owen
 */

#ifndef GCTAEDISPPERFTABLE_HPP
#define GCTAEDISPPERFTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GRan.hpp"
#include "GNodeArray.hpp"
#include "GCTAEdisp.hpp"


/***********************************************************************//**
 * @class GCTAEdispPerfTable
 *
 * @brief CTA performance table energy dispersion class
 *
 * This class implements the CTA energy dispersion response as function
 * of energy as determined from a performance table. The performance table is
 * an ASCII file that specifies the CTA performance parameters in a simple
 * way.
 ***************************************************************************/
class GCTAEdispPerfTable : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispPerfTable(void);
    explicit GCTAEdispPerfTable(const std::string& filename);
    GCTAEdispPerfTable(const GCTAEdispPerfTable& psf);
    virtual ~GCTAEdispPerfTable(void);

    // Operators
    GCTAEdispPerfTable& operator=(const GCTAEdispPerfTable& psf);
    double operator()(const double& logEobs,
                      const double& logEsrc,
                      const double& theta = 0.0,
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;

    // Implemented pure virtual methods
    void                clear(void);
    GCTAEdispPerfTable* clone(void) const;
    std::string         classname(void) const;
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

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdispPerfTable& psf);
    void free_members(void);
    void update(const double& logE) const;

    // Members
    std::string         m_filename;  //!< Name of response file
    GNodeArray          m_logE;      //!< log(E) nodes for interpolation
    std::vector<double> m_sigma;     //!< Sigma value (rms) of energy resolution

    // Precomputation cache
    mutable double      m_par_logE;  //!< Energy for which precomputation is done
    mutable double      m_par_scale; //!< Gaussian normalization
    mutable double      m_par_sigma; //!< Gaussian sigma
    mutable double      m_par_width; //!< Gaussian width parameter
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEdispPerfTable").
 ***************************************************************************/
inline
std::string GCTAEdispPerfTable::classname(void) const
{
    return ("GCTAEdispPerfTable");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the energy resolution was loaded
 ***************************************************************************/
inline
std::string GCTAEdispPerfTable::filename(void) const
{
    return m_filename;
}

#endif /* GCTAEDISPPERFTABLE_HPP */
