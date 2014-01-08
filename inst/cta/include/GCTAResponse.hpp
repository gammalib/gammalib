/***************************************************************************
 *                  GCTAResponse.hpp - CTA Response class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAResponse.hpp
 * @brief CTA instrument response function class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSE_HPP
#define GCTARESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GResponse.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GSkyDir;
class GPhoton;
class GEbounds;
class GEvent;
class GObservation;
class GCTAObservation;
class GCTAPointing;
class GCTAEventAtom;
class GCTARoi;
class GCTAInstDir;
class GCTAAeff;
class GCTAPsf;
class GCTAEdisp;


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief CTA instrument response function class
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    GCTAResponse(const std::string& rspname, const std::string& caldb = "");
    virtual ~GCTAResponse(void);

    // Operators
    virtual GCTAResponse& operator=(const GCTAResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GCTAResponse* clone(void) const;
    virtual bool          has_edisp(void) const;
    virtual bool          has_tdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Overload virtual base class methods
    virtual double irf_radial(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const;
    virtual double irf_elliptical(const GEvent&       event,
                                  const GSource&      source,
                                  const GObservation& obs) const;
    virtual double irf_diffuse(const GEvent&       event,
                               const GSource&      source,
                               const GObservation& obs) const;
    virtual double npred_radial(const GSource&      source,
                                const GObservation& obs) const;
    virtual double npred_elliptical(const GSource&      source,
                                    const GObservation& obs) const;
    virtual double npred_diffuse(const GSource&      source,
                                 const GObservation& obs) const;

    // Other Methods
    GCTAEventAtom*     mc(const double& area, const GPhoton& photon,
                          const GObservation& obs, GRan& ran) const;
    void               caldb(const std::string& caldb);
    const std::string& caldb(void) const;
    void               load(const std::string& rspname);
    void               eps(const double& eps);
    const double&      eps(void) const;
    const std::string& rmffile(void) const;
    void               load_aeff(const std::string& filename);
    void               load_psf(const std::string& filename);
    void               load_edisp(const std::string& filename);
    void               offset_sigma(const double& sigma);
    double             offset_sigma(void) const;
    const GCTAAeff*    aeff(void) const;
    void               aeff(GCTAAeff* aeff);
    const GCTAPsf*     psf(void) const;
    void               psf(GCTAPsf* psf);
    const GCTAEdisp*   edisp(void) const;
    void               edisp(GCTAEdisp* edisp);

    // Low-level response methods
    double aeff(const double& theta,
                const double& phi,
                const double& zenith,
                const double& azimuth,
                const double& srcLogEng) const;
    double psf(const double& delta,
               const double& theta,
               const double& phi,
               const double& zenith,
               const double& azimuth,
               const double& srcLogEng) const;
    double psf_delta_max(const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const double& srcLogEng) const;
    double edisp(const double& obsLogEng,
                 const double& theta,
                 const double& phi,
                 const double& zenith,
                 const double& azimuth,
                 const double& srcLogEng) const;
    double npsf(const GSkyDir&      srcDir,
                const double&       srcLogEng,
                const GTime&        srcTime,
                const GCTAPointing& pnt,
                const GCTARoi&      roi) const;
    double nedisp(const GSkyDir&      srcDir,
                  const GEnergy&      srcEng,
                  const GTime&        srcTime,
                  const GCTAPointing& pnt,
                  const GEbounds&     ebds) const;

private:
    // Private methods
    void                   init_members(void);
    void                   copy_members(const GCTAResponse& rsp);
    void                   free_members(void);
    const GCTAObservation& retrieve_obs(const std::string& origin,
                                        const GObservation& obs) const;
    const GCTAPointing&    retrieve_pnt(const std::string& origin,
                                        const GObservation& obs) const;
    const GCTARoi&         retrieve_roi(const std::string& origin,
                                        const GObservation& obs) const;
    const GCTAInstDir&     retrieve_dir(const std::string& origin,
                                        const GEvent&      event) const;

    // Private data members
    std::string         m_caldb;    //!< Name of or path to the calibration database
    std::string         m_rspname;  //!< Name of the instrument response
    std::string         m_rmffile;  //!< Name of RMF file
    double              m_eps;      //!< Integration precision
    GCTAAeff*           m_aeff;     //!< Effective area
    GCTAPsf*            m_psf;      //!< Point spread function
    GCTAEdisp*          m_edisp;    //!< Energy dispersion

    // Npred cache
    mutable std::vector<std::string> m_npred_names;    //!< Model names
    mutable std::vector<GEnergy>     m_npred_energies; //!< Model energy
    mutable std::vector<GTime>       m_npred_times;    //!< Model time
    mutable std::vector<double>      m_npred_values;   //!< Model values
};


/***********************************************************************//**
 * @brief Signal if energy dispersion has been implemented
 *
 * @return False.
 ***************************************************************************/
inline
bool GCTAResponse::has_edisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if time dispersion has been implemented
 *
 * @return False.
 ***************************************************************************/
inline
bool GCTAResponse::has_tdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Return calibration database
 *
 * @return Calibration database.
 ***************************************************************************/
inline
const std::string& GCTAResponse::caldb(void) const
{
    return m_caldb;
}


/***********************************************************************//**
 * @brief Set computation precision
 *
 * @param[in] eps Computation precision.
 ***************************************************************************/
inline
void GCTAResponse::eps(const double& eps)
{
    m_eps = eps;
    return;
}


/***********************************************************************//**
 * @brief Return computation precision
 *
 * @return Computation precision.
 ***************************************************************************/
inline
const double& GCTAResponse::eps(void) const
{
    return m_eps;
}


/***********************************************************************//**
 * @brief Return Redestribution Matrix File name
 *
 * @return Redestribution Matrix File name.
 ***************************************************************************/
inline
const std::string& GCTAResponse::rmffile(void) const
{
    return m_rmffile;
}


/***********************************************************************//**
 * @brief Load energy dispersion information
 *
 * @param[in] filename Energy dispersion file name.
 ***************************************************************************/
inline
void GCTAResponse::load_edisp(const std::string& filename)
{
    return;
}


/***********************************************************************//**
 * @brief Return pointer to effective area response
 *
 * @return Pointer to effective area response.
 ***************************************************************************/
inline
const GCTAAeff* GCTAResponse::aeff(void) const
{
    return m_aeff;
}


/***********************************************************************//**
 * @brief Set pointer to effective area response
 *
 * @param[in] aeff Pointer to effective area response.
 ***************************************************************************/
inline
void GCTAResponse::aeff(GCTAAeff* aeff)
{
    m_aeff = aeff;
    return;
}


/***********************************************************************//**
 * @brief Return pointer to point spread function
 *
 * @return Pointer to point spread function.
 ***************************************************************************/
inline
const GCTAPsf* GCTAResponse::psf(void) const
{
    return m_psf;
}


/***********************************************************************//**
 * @brief Set pointer to point spread function
 *
 * @param[in] psf Pointer to point spread function.
 ***************************************************************************/
inline
void GCTAResponse::psf(GCTAPsf* psf)
{
    m_psf = psf;
    return;
}


/***********************************************************************//**
 * @brief Return pointer to energy dispersion
 *
 * @return Pointer to energy dispersion.
 ***************************************************************************/
inline
const GCTAEdisp* GCTAResponse::edisp(void) const
{
    return m_edisp;
}


/***********************************************************************//**
 * @brief Set pointer to energy dispersion
 *
 * @param[in] edisp Pointer to energy dispersion.
 ***************************************************************************/
inline
void GCTAResponse::edisp(GCTAEdisp* edisp)
{
    m_edisp = edisp;
    return;
}

#endif /* GCTARESPONSE_HPP */
