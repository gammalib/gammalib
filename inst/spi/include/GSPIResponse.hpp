/***************************************************************************
 *              GSPIResponse.hpp - INTEGRAL/SPI response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIResponse.hpp
 * @brief INTEGRAL/SPI instrument response function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIRESPONSE_HPP
#define GSPIRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GResponse.hpp"
#include "GSkyMap.hpp"
#include "GNodeArray.hpp"
#include "GEbounds.hpp"
#include "GFilename.hpp"

/* __ Forward declaration ________________________________________________ */
class GEvent;
class GPhoton;
class GEnergy;
class GTime;
class GSource;
class GObservation;
class GModelSky;
class GEbounds;
class GSPIObservation;
class GSPIEventCube;
class GSPIEventBin;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIResponse
 *
 * @brief INTEGRAL/SPI instrument response function class
 *
 * The INTEGRAL/SPI instrument response function class defines the function
 * that translates from physical quantities to measured quantities.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GSPIResponse : public GResponse {

public:
    // Constructors and destructors
    GSPIResponse(void);
    GSPIResponse(const GSPIResponse& rsp);
    explicit GSPIResponse(const GFilename& rspname);
    virtual ~GSPIResponse(void);

    // Operators
    virtual GSPIResponse& operator=(const GSPIResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GSPIResponse* clone(void) const;
    virtual std::string   classname(void) const;
    virtual bool          use_edisp(void) const;
    virtual bool          use_tdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        nroi(const GModelSky&    model,
                               const GEnergy&      obsEng,
                               const GTime&        obsTime,
                               const GObservation& obs) const;
    virtual GEbounds      ebounds(const GEnergy& obsEnergy) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other Methods
    void             rspname(const GFilename& rspname);
    const GFilename& rspname(void) const;
    bool             is_precomputed(void) const;
    const double&    energy_keV(void) const;
    const double&    dlogE(void) const;
    const double&    gamma(void) const;
    void             set(const GSPIObservation& obs,
                         const GEnergy&         energy = GEnergy());
    double           irf_value(const GSkyDir&      srcDir,
                               const GSPIEventBin& bin,
                               const int&          ireg) const;
    double           zenith(const int& ipt, const GSkyDir& dir) const;
    double           azimuth(const int& ipt, const GSkyDir& dir) const;
    void             read(const GFits& fits);
    void             write(GFits& fits) const;
    void             load(const GFilename& filename);
    void             save(const GFilename& filename,
                          const bool&      clobber = false) const;

private:
    // Private methods
    void    init_members(void);
    void    copy_members(const GSPIResponse& rsp);
    void    free_members(void);
    void    read_detids(const GFits& fits);
    void    read_energies(const GFits& fits);
    void    write_detids(GFits& fits) const;
    void    write_energies(GFits& fits) const;
    void    load_irfs(const int& region);
    GSkyMap load_irf(const GFilename& rspname) const;
    GSkyMap compute_irf(const double& emin, const double& emax) const;
    void    set_wcs(const GFitsImage* image);
    void    set_detids(const GSPIEventCube* cube);
    void    set_cache(const GSPIEventCube* cube);
    int     irf_detid(const int& detid) const;
    double  irf_weight(const double& beta,
                       const double& emin,
                       const double& emax) const;

    // Private data members
    GFilename            m_rspname;    //!< File name of response group
    std::vector<int>     m_detids;     //!< Vector of detector IDs
    GNodeArray           m_energies;   //!< Node array of IRF energies
    GEbounds             m_ebounds;    //!< Energy bounaries of IRF
    GSkyMap              m_irfs;       //!< IRFs stored as sky map
    double               m_energy_keV; //!< IRF line energy (optional)
    double               m_dlogE;      //!< Logarithmic energy step for IRF band
    double               m_gamma;      //!< Power-law spectral index for IRF band

    // Private cache
    std::vector<GSkyDir> m_spix;         //!< SPI pointing direction
    std::vector<double>  m_posang;       //!< Position angle of Y axis (CEL, radians)
    mutable bool         m_has_wcs;      //!< Has WCS information
    mutable double       m_wcs_xmin;     //!< Minimum X value (radians)
    mutable double       m_wcs_ymin;     //!< Minimum Y value (radians)
    mutable double       m_wcs_xmax;     //!< Maximum X value (radians)
    mutable double       m_wcs_ymax;     //!< Maximum Y value (radians)
    mutable double       m_wcs_xbin;     //!< X value bin size (radians)
    mutable double       m_wcs_ybin;     //!< Y value bin size (radians)
    mutable double       m_wcs_xpix_max; //!< Maximum X pixel index
    mutable double       m_wcs_ypix_max; //!< Maximum Y pixel index
    mutable double       m_max_zenith;   //!< Maximum zenith angle (radians)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIResponse").
 ***************************************************************************/
inline
std::string GSPIResponse::classname(void) const
{
    return ("GSPIResponse");
}


/***********************************************************************//**
 * @brief Signal if energy dispersion will be used
 *
 * @return False.
 *
 * @todo Implement method as needed.
 ***************************************************************************/
inline
bool GSPIResponse::use_edisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if time dispersion will be used
 *
 * @return False.
 *
 * @todo Implement method as needed.
 ***************************************************************************/
inline
bool GSPIResponse::use_tdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Set response name
 *
 * @param[in] rspname Response group file name.
 *
 * Sets the response group file name.
 ***************************************************************************/
inline
void GSPIResponse::rspname(const GFilename& rspname)
{
    m_rspname = rspname;
    return;
}


/***********************************************************************//**
 * @brief Get response group file name
 *
 * @return Response group file name.
 *
 * Returns the response group file name.
 ***************************************************************************/
inline
const GFilename& GSPIResponse::rspname(void) const
{
    return m_rspname;
}


/***********************************************************************//**
 * @brief Signals if response is precomputed
 *
 * @return True if response is precomputed.
 *
 * Signals if the response was precomputed.
 ***************************************************************************/
inline
bool GSPIResponse::is_precomputed(void) const
{
    return (!m_ebounds.is_empty());
}


/***********************************************************************//**
 * @brief Return line IRF energy in keV
 *
 * @return Line IRF energy (keV).
 *
 * Returns the energy in keV for a line IRF. If the IRF is a continuum IRF
 * the method returns 0.
 ***************************************************************************/
inline
const double& GSPIResponse::energy_keV(void) const
{
    return (m_energy_keV);
}


/***********************************************************************//**
 * @brief Return logarithmic step size for continuum IRFs
 *
 * @return Logarithmic step size for continuum IRFs.
 *
 * Returns the logarithmic step size for the computation of continuum IRFs.
 ***************************************************************************/
inline
const double& GSPIResponse::dlogE(void) const
{
    return (m_dlogE);
}


/***********************************************************************//**
 * @brief Return power-law index for continuum IRFs
 *
 * @return Power-law index for continuum IRFs.
 *
 * Returns the power-law index for the computation of continuum IRFs.
 ***************************************************************************/
inline
const double& GSPIResponse::gamma(void) const
{
    return (m_gamma);
}


/***********************************************************************//**
 * @brief Return zenith angle of sky direction for pointing in radians
 *
 * @param[in] ipt Pointing index.
 * @param[in] dir Sky direction.
 * @return Zenith angle (radians).
 *
 * Returns zenith angle of sky direction for pointing in radians.
 ***************************************************************************/
inline
double GSPIResponse::zenith(const int& ipt, const GSkyDir& dir) const
{
    return (m_spix[ipt].dist(dir));
}


/***********************************************************************//**
 * @brief Return azimuth angle of sky direction for pointing in radians
 *
 * @param[in] ipt Pointing index.
 * @param[in] dir Sky direction.
 * @return Azimuth angle (radians).
 *
 * Returns azimuth angle of sky direction for pointing in radians.
 ***************************************************************************/
inline
double GSPIResponse::azimuth(const int& ipt, const GSkyDir& dir) const
{
    double azimuth = m_posang[ipt] - m_spix[ipt].posang(dir); // Celestial system
    if (azimuth < 0.0) {
        azimuth += gammalib::twopi;
    }
    return (azimuth);
}

#endif /* GSPIRESPONSE_HPP */
