/***************************************************************************
 *           GCTAPsf2D.hpp - CTA 2D point spread function class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAPsf2D.hpp
 * @brief CTA 2D point spread function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAPSF2D_HPP
#define GCTAPSF2D_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GRan.hpp"
#include "GCTAPsf.hpp"
#include "GCTAResponseTable.hpp"


/***********************************************************************//**
 * @class GCTAPsf2D
 *
 * @brief CTA 2D point spread function class
 *
 * This class implements the CTA point spread function response as function
 * of energy and offset angle.
 ***************************************************************************/
class GCTAPsf2D : public GCTAPsf {

public:
    // Constructors and destructors
    GCTAPsf2D(void);
    GCTAPsf2D(const std::string& filename);
    GCTAPsf2D(const GCTAPsf2D& psf);
    virtual ~GCTAPsf2D(void);

    // Operators
    GCTAPsf2D& operator=(const GCTAPsf2D& psf);
    double operator()(const double& delta,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void        clear(void);
    GCTAPsf2D*  clone(void) const;
    void        load(const std::string& filename);
    std::string filename(void) const;
    double      mc(GRan&         ran,
                   const double& logE, 
                   const double& theta = 0.0, 
                   const double& phi = 0.0,
                   const double& zenith = 0.0,
                   const double& azimuth = 0.0,
                   const bool&   etrue = true) const;
    double      delta_max(const double& logE, 
                          const double& theta = 0.0, 
                          const double& phi = 0.0,
                          const double& zenith = 0.0,
                          const double& azimuth = 0.0,
                          const bool&   etrue = true) const;
    std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAPsf2D& psf);
    void free_members(void);
    void update(const double& logE, const double& theta) const;

    // Members
    std::string       m_filename;   //!< Name of Aeff response file
    GCTAResponseTable m_psf;        //!< PSF response table

    // Precomputation cache
    mutable double    m_par_logE;   //!< Cache energy
    mutable double    m_par_theta;  //!< Cache offset angle
    mutable double    m_norm;       //!< Global normalization
    mutable double    m_norm2;      //!< Gaussian 2 normalization
    mutable double    m_norm3;      //!< Gaussian 3 normalization
    mutable double    m_sigma1;     //!< Gaussian 1 sigma
    mutable double    m_sigma2;     //!< Gaussian 2 sigma
    mutable double    m_sigma3;     //!< Gaussian 3 sigma
    mutable double    m_width1;     //!< Gaussian 1 width
    mutable double    m_width2;     //!< Gaussian 2 width
    mutable double    m_width3;     //!< Gaussian 3 width
};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which point spread function was loaded
 ***************************************************************************/
inline
std::string GCTAPsf2D::filename(void) const
{
    return m_filename;
}

#endif /* GCTAPSF2D_HPP */
