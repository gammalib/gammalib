/***************************************************************************
 *                GCTAAeffArf.hpp - CTA ARF effective area class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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


/***********************************************************************//**
 * @class GCTAAeffArf
 *
 * @brief CTA ARF effective area class
 *
 * This class implements the CTA effective area response that is defined by
 * an auxilliary response function (ARF) file.
 ***************************************************************************/
class GCTAAeffArf : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeffArf(void);
    GCTAAeffArf(const std::string& filename);
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
    std::string  print(void) const;

    // Methods
    int           size(void) const { return m_logE.size(); }
    void          sigma(const double& sigma) { m_sigma=sigma; }
    const double& sigma(void) const { return m_sigma; }
    void          thetacut(const double& thetacut) { m_thetacut=thetacut; }
    const double& thetacut(void) const { return m_thetacut; }
    void          scale(const double& scale) { m_scale=scale; }
    const double& scale(void) const { return m_scale; }
    void          read_arf(const GFitsTable* hdu);
    
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

#endif /* GCTAAEFFARF_HPP */
