/***************************************************************************
 *            GCTAEdisp2D.hpp - CTA 2D energy dispersion class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2017 by Florent Forest                              *
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
#include "GFilename.hpp"
#include "GFunction.hpp"
#include "GCTAEdisp.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GFits;
class GFitsBinTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_edisp2d = "ENERGY DISPERSION";
}


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
    explicit GCTAEdisp2D(const GFilename& filename);
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
    void         load(const GFilename& filename);
    GFilename    filename(void) const;
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
    void                     fetch(void) const;
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdisp2D& edisp);
    void free_members(void);
    void update(const double& logEobs,
                const double& logEsrc,
                const double& theta) const;
    void compute_ebounds_obs(const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    void compute_ebounds_src(const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    void set_indices(void);
    void set_boundaries(void);
    void set_max_edisp(void) const;
    void normalize_table(void);

    // Protected classes
    class edisp_kern : public GFunction {
    public:
        edisp_kern(const GCTAEdisp2D* parent,
                   const double&      logEsrc,
                   const double&      theta) :
                   m_parent(parent),
                   m_logEsrc(logEsrc),
                   m_theta(theta) { }
        double eval(const double& x);
    protected:
        const GCTAEdisp2D* m_parent;  //!< Pointer to parent class
        double             m_logEsrc; //!< True photon energy
        double             m_theta;   //!< Offset angle
    };

    // Members
    mutable GFilename m_filename;    //!< Name of Edisp response file
    GCTAResponseTable m_edisp;       //!< Edisp response table
    mutable bool      m_fetched;     //!< Signals that Edisp has been fetched
    int               m_inx_etrue;   //!< True energy index
    int               m_inx_migra;   //!< Migration index
    int               m_inx_theta;   //!< Theta index
    int               m_inx_matrix;  //!< Matrix
    double            m_logEsrc_min; //!< Minimum logE (log10(E/TeV))
    double            m_logEsrc_max; //!< Maximum logE (log10(E/TeV))
    double            m_migra_min;   //!< Minimum migration
    double            m_migra_max;   //!< Maximum migration
    double            m_theta_min;   //!< Minimum theta (radians)
    double            m_theta_max;   //!< Maximum theta (radians)

    // Computation cache
    mutable bool                  m_ebounds_obs_computed;
    mutable bool                  m_ebounds_src_computed;
    mutable double                m_last_theta_obs;
    mutable double                m_last_theta_src;
    mutable double                m_last_logEsrc;
    mutable double                m_last_logEobs;
    mutable int                   m_index_obs;
    mutable int                   m_index_src;
    mutable double                m_max_edisp;
    mutable std::vector<GEbounds> m_ebounds_obs;
    mutable std::vector<GEbounds> m_ebounds_src;
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
GFilename GCTAEdisp2D::filename(void) const
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
    fetch();
    return (m_edisp);
}

#endif /* GCTAEDISP2D_HPP */
