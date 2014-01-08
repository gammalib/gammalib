/***************************************************************************
 *                GCTAAeffArf.hpp - CTA ARF effective area class           *
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
 * @file GCTAAeffArf.hpp
 * @brief CTA ARF effective area class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAAEFFARF_HPP
#define GCTAAEFFARF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GFits.hpp"
#include "GNodeArray.hpp"
#include "GFitsTable.hpp"
#include "GCTAAeff.hpp"

/* __ Forward declarations _______________________________________________ */
class GCTAResponse;


/***********************************************************************//**
 * @class GCTAAeffArf
 *
 * @brief CTA ARF effective area class
 *
 * This class implements the CTA effective area response that is defined by
 * an auxiliary response function (ARF) file.
 ***************************************************************************/
class GCTAAeffArf : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeffArf(void);
    explicit GCTAAeffArf(const std::string& filename);
    GCTAAeffArf(const GCTAAeffArf& cta);
    virtual ~GCTAAeffArf(void);

    // Operators
    GCTAAeffArf& operator=(const GCTAAeffArf& aeff);
    double operator()(const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void         clear(void);
    GCTAAeffArf* clone(void) const;
    void         load(const std::string& filename);
    std::string  filename(void) const;
    std::string  print(const GChatter& chatter = NORMAL) const;

    // Methods
    int           size(void) const;
    void          sigma(const double& sigma);
    const double& sigma(void) const;
    void          thetacut(const double& thetacut);
    const double& thetacut(void) const;
    void          scale(const double& scale);
    const double& scale(void) const;
    void          read(const GFitsTable& hdu);
    void          remove_thetacut(const GCTAResponse& rsp);
    
private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAAeffArf& aeff);
    void free_members(void);

    // Members
    std::string         m_filename;  //!< Name of Aeff response file
    GNodeArray          m_logE;      //!< log(E) nodes for Aeff interpolation
    std::vector<double> m_aeff;      //!< Effective area in cm2
    double              m_sigma;     //!< Sigma for offset angle computation (0=none)
    double              m_thetacut;  //!< Theta cut for ARF
    double              m_scale;     //!< Scale for ARF
};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which effective area was loaded.
 ***************************************************************************/
inline
std::string GCTAAeffArf::filename(void) const
{
    return m_filename;
}


/***********************************************************************//**
 * @brief Return number of node energies in response
 *
 * @return Number of node energies.
 ***************************************************************************/
inline
int GCTAAeffArf::size(void) const
{
    return (m_logE.size());
}


/***********************************************************************//**
 * @brief Set sigma for offset angle dependence
 *
 * @param[in] sigma Sigma for offset angle dependence.
 *
 * Sets the sigma parameter for the offset angle dependence of the effective
 * area. If @p sigma is 0 (which is the default value), no offset angle
 * dependency will be assumed.
 ***************************************************************************/
inline
void GCTAAeffArf::sigma(const double& sigma)
{
    m_sigma = sigma;
    return;
}


/***********************************************************************//**
 * @brief Return sigma for offset angle dependence
 *
 * @return Sigma for offset angle dependence.
 ***************************************************************************/
inline
const double& GCTAAeffArf::sigma(void) const
{
    return (m_sigma);
}


/***********************************************************************//**
 * @brief Set theta cut angle
 *
 * @param[in] thetacut Set theta cut angle.
 *
 * Sets the theta cut angle which defines the energy independent cut that
 * has been assumed in deriving the ARF values. If @p thetacut os 0 (which
 * is the default value), not thetacut will be applied.
 ***************************************************************************/
inline
void GCTAAeffArf::thetacut(const double& thetacut)
{
    m_thetacut = thetacut;
    return;
}


/***********************************************************************//**
 * @brief Return theta cut angle
 *
 * @return Theta cut angle.
 ***************************************************************************/
inline
const double& GCTAAeffArf::thetacut(void) const
{
    return (m_thetacut);
}


/***********************************************************************//**
 * @brief Set effective area scaling factor
 *
 * @param[in] scale Set effective area scaling factor.
 *
 * Sets the scaling factor for the effective area (by default is 1).
 ***************************************************************************/
inline
void GCTAAeffArf::scale(const double& scale)
{
    m_scale = scale;
    return;
}


/***********************************************************************//**
 * @brief Return effective area scaling factor
 *
 * @return Effective area scaling factor.
 ***************************************************************************/
inline
const double& GCTAAeffArf::scale(void) const
{
    return (m_scale);
}

#endif /* GCTAAEFFARF_HPP */
