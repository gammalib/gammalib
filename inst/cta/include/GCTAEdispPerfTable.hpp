/***************************************************************************
 *  GCTAEdispPerfTable.hpp - CTA performance table energy dispersion class *
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
 * @file GCTAEdispPerfTable.hpp
 * @brief CTA performance table energy dispersion class definition
 * @author Christoph Deil & Ellis Owen
 */

#ifndef GCTAEDISPPERFTABLE_HPP
#define GCTAEDISPPERFTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GFilename.hpp"
#include "GNodeArray.hpp"
#include "GCTAEdisp.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;


/***********************************************************************//**
 * @class GCTAEdispPerfTable
 *
 * @brief CTA performance table energy dispersion class
 *
 * This class implements the CTA energy dispersion response as function
 * of true energy as determined from a performance table. A performance
 * table is an ASCII file that represent a table, specifying the CTA
 * performance parameters as function of true photon energy.
 *
 * The energy dispersion is defined as
 *
 * \f[
 *    E_{\rm disp}(E_{\rm true}, E_{\rm reco}) =
 *    \frac{1}{\sqrt{2\pi}\sigma(E_{\rm true})}
 *    \exp \left(\frac{-(\log_{10} E_{\rm reco} - \log_{10} E_{\rm true})^2}
 *                    {2 \sigma(E_{\rm true})^2} \right) \times
 *    \frac{1}{\log_{10} E_{\rm reco}}
 * \f]
 *
 * and given in units of MeV\f$^{-1}\f$, where
 * \f$E_{\rm reco}\f$ is the reconstructed energy in units of MeV,
 * \f$E_{\rm true}\f$ is the true energy in units of MeV, and
 * \f$\sigma(E_{\rm true})\f$ is the standard deviation of the energy
 * dispersion that depends on the true photon energy.
 ***************************************************************************/
class GCTAEdispPerfTable : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispPerfTable(void);
    explicit GCTAEdispPerfTable(const GFilename& filename);
    GCTAEdispPerfTable(const GCTAEdispPerfTable& psf);
    virtual ~GCTAEdispPerfTable(void);

    // Operators
    GCTAEdispPerfTable& operator=(const GCTAEdispPerfTable& psf);
    double operator()(const GEnergy& ereco,
                      const GEnergy& etrue,
                      const double&  theta = 0.0,
                      const double&  phi = 0.0,
                      const double&  zenith = 0.0,
                      const double&  azimuth = 0.0) const;

    // Implemented pure virtual methods
    void                clear(void);
    GCTAEdispPerfTable* clone(void) const;
    std::string         classname(void) const;
    void                load(const GFilename& filename);
    GFilename           filename(void) const;
    GEnergy             mc(GRan&          ran,
                           const GEnergy& etrue,
                           const double&  theta = 0.0,
                           const double&  phi = 0.0,
                           const double&  zenith = 0.0,
                           const double&  azimuth = 0.0) const;
    GEbounds            ereco_bounds(const GEnergy& etrue,
                                     const double&  theta = 0.0,
                                     const double&  phi = 0.0,
                                     const double&  zenith = 0.0,
                                     const double&  azimuth = 0.0) const;
    GEbounds            etrue_bounds(const GEnergy& ereco,
                                     const double&  theta = 0.0,
                                     const double&  phi = 0.0,
                                     const double&  zenith = 0.0,
                                     const double&  azimuth = 0.0) const;
    double              prob_erecobin(const GEnergy& ereco_min,
                                      const GEnergy& ereco_max,
                                      const GEnergy& etrue,
                                      const double&  theta) const;
    std::string         print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAEdispPerfTable& psf);
    void free_members(void);
    void update(const double& logE) const;

    // Members
    mutable GFilename   m_filename;  //!< Name of response file
    GNodeArray          m_logE;      //!< log(E) nodes for interpolation
    std::vector<double> m_sigma;     //!< Sigma value (rms) of energy resolution

    // Precomputation cache
    mutable double      m_par_logE;  //!< Energy for which precomputation is done
    mutable double      m_par_scale; //!< Gaussian normalization
    mutable double      m_par_sigma; //!< Gaussian sigma
    mutable double      m_par_width; //!< Gaussian width parameter
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEdispPerfTable").
 ***************************************************************************/
inline
std::string GCTAEdispPerfTable::classname(void) const
{
    return ("GCTAEdispPerfTable");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the energy resolution was loaded
 ***************************************************************************/
inline
GFilename GCTAEdispPerfTable::filename(void) const
{
    return m_filename;
}

#endif /* GCTAEDISPPERFTABLE_HPP */
