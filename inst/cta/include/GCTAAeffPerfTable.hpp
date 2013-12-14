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
    explicit GCTAAeffPerfTable(const std::string& filename);
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
    int           size(void) const;
    void          sigma(const double& sigma);
    const double& sigma(void) const;

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


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which effective area was loaded.
 ***************************************************************************/
inline
std::string GCTAAeffPerfTable::filename(void) const
{
    return m_filename;
}


/***********************************************************************//**
 * @brief Return number of node energies in response
 *
 * @return Number of node energies.
 ***************************************************************************/
inline
int GCTAAeffPerfTable::size(void) const
{
    return (m_logE.size());
}


/***********************************************************************//**
 * @brief Set sigma for offset angle dependence
 *
 * @param[in] sigma Sigma for offset angle dependence.
 *
 * Sets the sigma parameter for the offset angle dependence of the effective
 * area. If @p sigma is 0, no offset angle dependency will be assumed. By
 * default, @p sigma = 3.
 ***************************************************************************/
inline
void GCTAAeffPerfTable::sigma(const double& sigma)
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
const double& GCTAAeffPerfTable::sigma(void) const
{
    return (m_sigma);
}

#endif /* GCTAAEFFPERFTABLE_HPP */
