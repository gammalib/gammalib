/***************************************************************************
 *    GCTAAeffPerfTable.hpp - CTA performance table effective area class   *
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
 * @file GCTAAeffPerfTable.hpp
 * @brief CTA performance table effective area class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAAEFFPERFTABLE_HPP
#define GCTAAEFFPERFTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GCTAAeff.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GCTAAeffPerfTable
 *
 * @brief CTA performance table effective area class
 *
 * This class implements the CTA effective area response as function of
 * energy as determined from a performance table. The performance table is
 * an ASCII file that specifies the CTA performance parameters in a simple
 * way.
 ***************************************************************************/
class GCTAAeffPerfTable : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeffPerfTable(void);
    GCTAAeffPerfTable(const std::string& filename);
    GCTAAeffPerfTable(const GCTAAeffPerfTable& cta);
    virtual ~GCTAAeffPerfTable(void);

    // Operators
    GCTAAeffPerfTable& operator=(const GCTAAeffPerfTable& aeff);
    double operator()(const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void               clear(void);
    GCTAAeffPerfTable* clone(void) const;
    void               load(const std::string& filename);
    std::string        filename(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

    // Methods
    int           size(void) const { return m_logE.size(); }
    void          sigma(const double& sigma) { m_sigma=sigma; }
    const double& sigma(void) const { return m_sigma; }

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAAeffPerfTable& aeff);
    void free_members(void);

    // Members
    std::string         m_filename;  //!< Name of Aeff response file
    GNodeArray          m_logE;      //!< log(E) nodes for Aeff interpolation
    std::vector<double> m_aeff;      //!< Effective area in cm2
    double              m_sigma;     //!< Sigma for offset angle computation (0=none)
};

#endif /* GCTAAEFFPERFTABLE_HPP */
