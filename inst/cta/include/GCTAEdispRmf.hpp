/***************************************************************************
 *             GCTAEdispRmf.hpp - CTA RMF energy dispersion class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Christoph Deil & Ellis Owen                 *
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
#include "GRmf.hpp"
#include "GVector.hpp"
#include "GFunction.hpp"
#include "GCTAEdisp.hpp"
#include "GNodeArray.hpp"

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
    std::string   classname(void) const;
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
    int         size(void) const;
    const GRmf& rmf(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdispRmf& psf);
    void free_members(void);
    void set_matrix(void);
    void set_cache(void) const;
    void set_max_edisp(void) const;
    void update(const double& arg1, const double& arg2) const;
    void compute_ebounds_obs(const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    void compute_ebounds_src(const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;

    // Protected classes
    class edisp_kern : public GFunction {
    public:
        edisp_kern(const GCTAEdispRmf* parent,
                   const double&       logEsrc,
                   const double&       theta) :
                   m_parent(parent),
                   m_logEsrc(logEsrc),
                   m_theta(theta) { }
        double eval(const double& x);
    protected:
        const GCTAEdispRmf* m_parent;  //!< Pointer to parent class
        double              m_logEsrc; //!< True photon energy
        double              m_theta;   //!< Offset angle
    };

    // Members
    std::string   m_filename;  //!< Name of response file
    GRmf          m_rmf;       //!< Redistribution matrix file
    GMatrixSparse m_matrix;    //!< Normalised redistribution matrix

    // Interpolation cache
    mutable GNodeArray m_etrue;          //!< Array of log10(Etrue)
    mutable GNodeArray m_emeasured;      //!< Array of log10(Emeasured)
    mutable double     m_last_etrue;     //!< Last log10(Etrue)
    mutable double     m_last_emeasured; //!< Last log10(Emeasured)
    mutable int        m_itrue1;         //!< Index of left Etrue
    mutable int        m_itrue2;         //!< Index of right Etrue
    mutable int        m_imeas1;         //!< Index of left Emeasured
    mutable int        m_imeas2;         //!< Index of right Emeasured
    mutable double     m_wgt1;           //!< Weight of lower left node
    mutable double     m_wgt2;           //!< Weight of upper left node
    mutable double     m_wgt3;           //!< Weight of lower right node
    mutable double     m_wgt4;           //!< Weight of upper right node

    // Monte Carlo cache
    mutable double                m_max_edisp;
    mutable double                m_last_theta_obs;
    mutable double                m_last_theta_src;
    mutable double                m_last_logEsrc;
    mutable double                m_last_logEobs;
    mutable int                   m_index_obs;
    mutable int                   m_index_src;
    mutable bool                  m_ebounds_obs_computed;
    mutable bool                  m_ebounds_src_computed;
    mutable std::vector<GEbounds> m_ebounds_obs;
    mutable std::vector<GEbounds> m_ebounds_src;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEdispRmf").
 ***************************************************************************/
inline
std::string GCTAEdispRmf::classname(void) const
{
    return ("GCTAEdispRmf");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the Redistribution Matrix was loaded.
 ***************************************************************************/
inline
std::string GCTAEdispRmf::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Return Redistribution Matrix File
 *
 * @return Reference to Redistribution Matrix File.
 ***************************************************************************/
inline
const GRmf& GCTAEdispRmf::rmf(void) const
{
    return (m_rmf);
}

#endif /* GCTAEDISPRMF_HPP */
