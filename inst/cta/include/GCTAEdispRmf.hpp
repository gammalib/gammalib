/***************************************************************************
 *             GCTAEdispRmf.hpp - CTA RMF energy dispersion class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Christoph Deil & Ellis Owen                 *
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

/* __ Includes (members required) ________________________________________ */
#include <string>
#include "GFilename.hpp"
#include "GRmf.hpp"
#include "GMatrixSparse.hpp"
#include "GNodeArray.hpp"
#include "GCTAEdisp.hpp"

/* __ Forward declarations _______________________________________________ */
class GFunction;
class GRan;
class GEnergy;
class GEbounds;


/***********************************************************************//**
 * @class GCTAEdispRmf
 *
 * @brief CTA Redistribution Matrix File (RMF) energy dispersion class
 *
 * The energy dispersion is defined as
 *
 * \f[
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true})
 * \f]
 *
 * in units of MeV\f$^{-1}\f$ where
 * \f$E_{\rm reco}\f$ is the reconstructed energy, and
 * \f$E_{\rm true}\f$ is the true energy.
 ***************************************************************************/
class GCTAEdispRmf : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispRmf(void);
    explicit GCTAEdispRmf(const GFilename& filename);
    GCTAEdispRmf(const GCTAEdispRmf& edisp);
    virtual ~GCTAEdispRmf(void);

    // Operators
    GCTAEdispRmf& operator=(const GCTAEdispRmf& edisp);
    double operator()(const GEnergy& ereco,
                      const GEnergy& etrue,
                      const double&  theta = 0.0,
                      const double&  phi = 0.0,
                      const double&  zenith = 0.0,
                      const double&  azimuth = 0.0) const;

    // Implemented pure virtual methods
    void          clear(void);
    GCTAEdispRmf* clone(void) const;
    std::string   classname(void) const;
    void          load(const GFilename& filename);
    GFilename     filename(void) const;
    GEnergy       mc(GRan& ran,
                     const GEnergy& etrue,
                     const double&  theta = 0.0,
                     const double&  phi = 0.0,
                     const double&  zenith = 0.0,
                     const double&  azimuth = 0.0) const;
    GEbounds      ereco_bounds(const GEnergy& etrue,
                               const double&  theta = 0.0,
                               const double&  phi = 0.0,
                               const double&  zenith = 0.0,
                               const double&  azimuth = 0.0) const;
    GEbounds      etrue_bounds(const GEnergy& ereco,
                               const double&  theta = 0.0,
                               const double&  phi = 0.0,
                               const double&  zenith = 0.0,
                               const double&  azimuth = 0.0) const;
    double        prob_erecobin(const GEnergy& ereco_min,
                                const GEnergy& ereco_max,
                                const GEnergy& etrue,
                                const double&  theta) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GRmf& rmf(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdispRmf& psf);
    void free_members(void);
    void set_matrix(void);
    void set_max_edisp(void);
    void set_cache(void) const;
    void update(const GEnergy& ereco, const GEnergy& etrue) const;
    void compute_ereco_bounds(void) const;
    void compute_etrue_bounds(void) const;

    // Members
    GFilename     m_filename;  //!< Name of response file
    GRmf          m_rmf;       //!< Redistribution matrix file
    GMatrixSparse m_matrix;    //!< Normalised redistribution matrix
    double        m_max_edisp; //!< Maximum energy dispersion value for MC

    // Interpolation cache
    mutable GNodeArray m_etrue;      //!< Array of log10(Etrue)
    mutable GNodeArray m_ereco;      //!< Array of log10(Ereco)
    mutable GEnergy    m_last_etrue; //!< Last true energy
    mutable GEnergy    m_last_ereco; //!< Last reconstructed energy
    mutable int        m_itrue1;     //!< Index of left Etrue
    mutable int        m_itrue2;     //!< Index of right Etrue
    mutable int        m_ireco1;     //!< Index of left Ereco
    mutable int        m_ireco2;     //!< Index of right Ereco
    mutable double     m_wgt1;       //!< Weight of lower left node
    mutable double     m_wgt2;       //!< Weight of upper left node
    mutable double     m_wgt3;       //!< Weight of lower right node
    mutable double     m_wgt4;       //!< Weight of upper right node

    // Computation cache
    mutable bool                  m_ereco_bounds_computed;
    mutable bool                  m_etrue_bounds_computed;
    mutable GEnergy               m_last_etrue_bounds;
    mutable GEnergy               m_last_ereco_bounds;
    mutable int                   m_index_ereco;
    mutable int                   m_index_etrue;
    mutable std::vector<GEbounds> m_ereco_bounds;
    mutable std::vector<GEbounds> m_etrue_bounds;
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
GFilename GCTAEdispRmf::filename(void) const
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
