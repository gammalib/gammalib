/***************************************************************************
 *            GCTAEdisp2D.hpp - CTA 2D energy dispersion class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Florent Forest                                   *
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
 * @file GCTAEdisp2D.hpp
 * @brief CTA 2D energy dispersion class definition
 * @author Florent Forest
 */

#ifndef GCTAEDISP2D_HPP
#define GCTAEDISP2D_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
//#include "GFits.hpp"
#include "GCTAEdisp.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GFits;
class GFitsBinTable;


/***********************************************************************//**
 * @class GCTAEdisp2D
 *
 * @brief CTA 2D energy dispersion class
 *
 * This class implements the energy dispersion for the CTA 2D response.
 ***************************************************************************/
class GCTAEdisp2D : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdisp2D(void);
    explicit GCTAEdisp2D(const std::string& filename);
    GCTAEdisp2D(const GCTAEdisp2D& edisp);
    virtual ~GCTAEdisp2D(void);

    // Operators
    GCTAEdisp2D& operator=(const GCTAEdisp2D& edisp);
    double operator()(const double& logEobs, 
                      const double& logEsrc, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;

    // Implemented methods
    void         clear(void);
    GCTAEdisp2D* clone(void) const;
    std::string  classname(void) const;
    void         load(const std::string& filename);
    std::string  filename(void) const;
    GEnergy      mc(GRan&         ran,
                    const double& logE,
                    const double& theta = 0.0,
                    const double& phi = 0.0,
                    const double& zenith = 0.0,
                    const double& azimuth = 0.0) const;
    GEbounds     ebounds_obs(const double& logEsrc,
                             const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    GEbounds     ebounds_src(const double& logEobs,
                             const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    std::string  print(const GChatter& chatter = NORMAL) const;

    // Methods
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFits& file);
    void                     write(GFitsBinTable& hdu) const;
    void                     save(const std::string& filename,
                                  const bool& clobber = false) const;
private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdisp2D& edisp);
    void free_members(void);
    void update(const double& logEobs,
                const double& logEsrc,
                const double& theta) const;
    void compute_cumul(const double& theta = 0.0,
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;
    void compute_ebounds_obs(const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    void compute_ebounds_src(const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    void set_fmax(void) const;

    // Members
    std::string       m_filename;  //!< Name of Edisp response file
    GCTAResponseTable m_edisp;     //!< Edisp response table

    mutable bool      m_cdf_computed;
    mutable double    m_theta;
    mutable bool      m_ebounds_obs_computed;
    mutable bool      m_ebounds_src_computed;
    mutable double    m_logEsrc;
    mutable double    m_logEobs;
    mutable int       m_index_obs;
    mutable int       m_index_src;
    mutable double    m_fmax;

    mutable std::vector<GEbounds> m_ebounds_obs;
    mutable std::vector<GEbounds> m_ebounds_src;
    mutable std::vector<std::vector<std::pair<double, double> > > m_cumul;
    mutable std::vector<double> m_eobs_axis;
    mutable std::vector<std::vector<double> > m_temp;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEdisp2D").
 ***************************************************************************/
inline
std::string GCTAEdisp2D::classname(void) const
{
    return ("GCTAEdisp2D");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Name of FITS file from which energy dispersion was loaded.
 ***************************************************************************/
inline
std::string GCTAEdisp2D::filename(void) const
{
    // Return filename
    return (m_filename);
}


/***********************************************************************//**
 * @brief Return response table
 *
 * @return Reference to response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAEdisp2D::table(void) const
{
    return (m_edisp);
}


/***********************************************************************//**
 * @brief Set response table
 *
 * @param[in] Response table.
 ***************************************************************************/
inline
void GCTAEdisp2D::table(const GCTAResponseTable& table)
{
    m_edisp = table;
    return;
}

#endif /* GCTAEDISP2D_HPP */
