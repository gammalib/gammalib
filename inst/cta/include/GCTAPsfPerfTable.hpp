/***************************************************************************
 *          GCTAPsfPerfTable.hpp - CTA performance table PSF class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAPsfPerfTable.hpp
 * @brief CTA performance table point spread function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAPSFPERFTABLE_HPP
#define GCTAPSFPERFTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GRan.hpp"
#include "GNodeArray.hpp"
#include "GCTAPsf.hpp"


/***********************************************************************//**
 * @class GCTAPsfPerfTable
 *
 * @brief CTA performance table point spread function class
 *
 * This class implements the CTA point spread function response as function
 * of energy as determined from a performance table. The performance table is
 * an ASCII file that specifies the CTA performance parameters in a simple
 * way.
 ***************************************************************************/
class GCTAPsfPerfTable : public GCTAPsf {

public:
    // Constructors and destructors
    GCTAPsfPerfTable(void);
    GCTAPsfPerfTable(const std::string& filename);
    GCTAPsfPerfTable(const GCTAPsfPerfTable& psf);
    virtual ~GCTAPsfPerfTable(void);

    // Operators
    GCTAPsfPerfTable& operator=(const GCTAPsfPerfTable& psf);
    double operator()(const double& delta,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void              clear(void);
    GCTAPsfPerfTable* clone(void) const;
    std::string       classname(void) const;
    void              load(const std::string& filename);
    std::string       filename(void) const;
    double            mc(GRan&         ran,
                         const double& logE, 
                         const double& theta = 0.0, 
                         const double& phi = 0.0,
                         const double& zenith = 0.0,
                         const double& azimuth = 0.0,
                         const bool&   etrue = true) const;
    double            delta_max(const double& logE, 
                                const double& theta = 0.0, 
                                const double& phi = 0.0,
                                const double& zenith = 0.0,
                                const double& azimuth = 0.0,
                                const bool&   etrue = true) const;
    std::string       print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAPsfPerfTable& psf);
    void free_members(void);
    void update(const double& logE) const;

    // Members
    std::string         m_filename;  //!< Name of Aeff response file
    GNodeArray          m_logE;      //!< log(E) nodes for Aeff interpolation
    std::vector<double> m_r68;       //!< 68% containment radius of PSF in degrees
    std::vector<double> m_r80;       //!< 80% containment radius of PSF in degrees
    std::vector<double> m_sigma;     //!< Sigma value of PSF in radians

    // Precomputation cache
    mutable double      m_par_logE;  //!< Energy for which precomputation is done
    mutable double      m_par_scale; //!< Gaussian normalization
    mutable double      m_par_sigma; //!< Gaussian sigma (radians)
    mutable double      m_par_width; //!< Gaussian width parameter
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAPsfPerfTable").
 ***************************************************************************/
inline
std::string GCTAPsfPerfTable::classname(void) const
{
    return ("GCTAPsfPerfTable");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which point spread function was loaded
 ***************************************************************************/
inline
std::string GCTAPsfPerfTable::filename(void) const
{
    return m_filename;
}

#endif /* GCTAPSFPERFTABLE_HPP */
