/***************************************************************************
 *   GCTABackgroundPerfTable.hpp - CTA performance table background class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackgroundPerfTable.hpp
 * @brief CTA performance table background class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTABACKGROUNDPERFTABLE_HPP
#define GCTABACKGROUNDPERFTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GCTABackground.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;


/***********************************************************************//**
 * @class GCTABackgroundPerfTable
 *
 * @brief CTA 3D background class
 ***************************************************************************/
class GCTABackgroundPerfTable : public GCTABackground {

public:
    // Constructors and destructors
    GCTABackgroundPerfTable(void);
    explicit GCTABackgroundPerfTable(const std::string& filename);
    GCTABackgroundPerfTable(const GCTABackgroundPerfTable& bgd);
    virtual ~GCTABackgroundPerfTable(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety,
                              const bool&   etrue = false) const;

    // Operators
    GCTABackgroundPerfTable& operator=(const GCTABackgroundPerfTable& bgd);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackgroundPerfTable*   clone(void) const;
    void                       load(const std::string& filename);
    std::string                filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    std::string                print(const GChatter& chatter = NORMAL) const;

    // Methods
    int           size(void) const;
    void          sigma(const double& sigma);
    const double& sigma(void) const;
    
private:
    // Methods
    void   init_members(void);
    void   copy_members(const GCTABackgroundPerfTable& bgd);
    void   free_members(void);
    double solidangle(void) const;
    void   init_mc_cache(void) const;

    // Radial integration class (used by solidangle() method). Note that
    // the sigma parameter is given in rad^2
    class integrand : public GFunction {
    public:
        integrand(const double& sigma) : m_sigma(sigma) { }
        double eval(const double& x) {
            double arg  = x * x / m_sigma;
            double arg2 = arg * arg;
            double f    = std::exp(-0.5 * arg2);
            return (f*std::sin(x));
        }
    private:
        const double& m_sigma;
    };

    // Members
    std::string         m_filename;   //!< Name of background response file
    GNodeArray          m_logE;       //!< log(E) nodes for background interpolation
    std::vector<double> m_background; //!< Background rate
    double              m_sigma;      //!< Sigma for offset angle computation (0=none)

    // Monte Carlo cache
    mutable GModelSpectralNodes m_mc_spectrum; //!< Background spectrum
};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the background was loaded.
 ***************************************************************************/
inline
std::string GCTABackgroundPerfTable::filename(void) const
{
    // Return filename
    return m_filename;
}


/***********************************************************************//**
 * @brief Return number of node energies in response
 *
 * @return Number of node energies.
 ***************************************************************************/
inline
int GCTABackgroundPerfTable::size(void) const
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
void GCTABackgroundPerfTable::sigma(const double& sigma)
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
const double& GCTABackgroundPerfTable::sigma(void) const
{
    return (m_sigma);
}


/***********************************************************************//**
 * @brief Get response cube spectrum
 *
 * @return Response cube spectrum.
 *
 * Returns the response cube spectrum.
 ***************************************************************************/
inline
const GModelSpectralNodes& GCTABackgroundPerfTable::spectrum(void) const
{
    if (m_mc_spectrum.nodes() == 0) {
        init_mc_cache();
    }
    return (m_mc_spectrum);
}

#endif /* GCTABACKGROUNDPERFTABLE_HPP */
