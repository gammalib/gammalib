/***************************************************************************
 *                  GCTAResponse.hpp  -  CTA Response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
 */

#ifndef GCTARESPONSE_HPP
#define GCTARESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include <vector>
#include "GMatrix.hpp"
#include "GEvent.hpp"
#include "GModelSky.hpp"
#include "GModelPointSource.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelDiffuseSource.hpp"
#include "GModelRadial.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"
#include "GInstDir.hpp"
#include "GRoi.hpp"
#include "GGti.hpp"
#include "GEbounds.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPhoton.hpp"
#include "GNodeArray.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAPointing.hpp"
#include "GCTARoi.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTADir.hpp"
#include "GCTAResponseTable.hpp"
#include "GCTAAeff.hpp"

/* __ Type definitions ___________________________________________________ */
typedef std::vector<double> GCTAPsfPars;

/* __ Forward declaration ________________________________________________ */
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief Interface for the CTA instrument response function
 *
 * The CTA response is still in the prototyping stage. To for allow
 * evolutions in the response formats, versions numbers have been defined.
 *
 * For the PSF, versions are stored in the m_psf_version member. The
 * following versions numbers are actually defined:
 * -99 : Unknown PSF version
 * -10 : PSF is given by ASCII file
 * -9  : PSF is given by vector
 * -8  : PSF is given by response table
 * Note that negative values are used so far as we're still in the
 * prototyping stage for the IRF definition. Once the IRFs are defined we'll
 * start with positive version numbers.
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    explicit GCTAResponse(const std::string& rspname, const std::string& caldb = "");
    virtual ~GCTAResponse(void);

    // Operators
    virtual GCTAResponse& operator= (const GCTAResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GCTAResponse* clone(void) const;
    virtual bool          hasedisp(void) const { return false; }
    virtual bool          hastdisp(void) const { return false; }
    virtual double        irf(const GInstDir&     obsDir,
                              const GEnergy&      obsEng,
                              const GTime&        obsTime,
                              const GSkyDir&      srcDir,
                              const GEnergy&      srcEng,
                              const GTime&        srcTime,
                              const GObservation& obs) const;
    virtual double        npred(const GSkyDir&      srcDir,
                                const GEnergy&      srcEng,
                                const GTime&        srcTime,
                                const GObservation& obs) const;
    virtual std::string   print(void) const;

    // Overload virtual base class methods
    virtual double irf_extended(const GInstDir&             obsDir,
                                const GEnergy&              obsEng,
                                const GTime&                obsTime,
                                const GModelExtendedSource& model,
                                const GEnergy&              srcEng,
                                const GTime&                srcTime,
                                const GObservation&         obs) const;
    virtual double irf_diffuse(const GInstDir&            obsDir,
                               const GEnergy&             obsEng,
                               const GTime&               obsTime,
                               const GModelDiffuseSource& model,
                               const GEnergy&             srcEng,
                               const GTime&               srcTime,
                               const GObservation&        obs) const;
    virtual double npred_extended(const GModelExtendedSource& model,
                                  const GEnergy&              srcEng,
                                  const GTime&                srcTime,
                                  const GObservation&         obs) const;
    virtual double npred_diffuse(const GModelDiffuseSource& model,
                                 const GEnergy&             srcEng,
                                 const GTime&               srcTime,
                                 const GObservation&        obs) const;

    // Other Methods
    GCTAEventAtom*  mc(const double& area, const GPhoton& photon,
                       const GObservation& obs, GRan& ran) const;
    void            caldb(const std::string& caldb);
    std::string     caldb(void) const { return m_caldb; }
    void            load(const std::string& rspname);
    void            eps(const double& eps) { m_eps=eps; }
    const double&   eps(void) const { return m_eps; }
    void            offset_sigma(const double& sigma) { m_offset_sigma=sigma; }
    const double&   offset_sigma(void) const { return m_offset_sigma; }
    std::string     rmffile(void) const { return m_rmffile; }
    std::string     psffile(void) const { return m_psffile; }
    void            read_psf(const GFitsTable* hdu);
    
    // New factorisation
    void            load_aeff(const std::string& filename);
    void            load_psf(const std::string& filename);
    void            load_edisp(const std::string& filename) {}
    const GCTAAeff* aeff(void) const { return m_aeff; }
    void            aeff(GCTAAeff* aeff) { m_aeff=aeff; }

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

    // Analytical PSF implementation
    double      psf_dummy(const double& delta,
                          const GCTAPsfPars& pars) const;
    GCTAPsfPars psf_dummy_sigma(const double& srcLogEng,
                                const double& theta) const;
    double      psf_dummy_max(const GCTAPsfPars& pars) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCTAResponse& rsp);
    void free_members(void);
    void read_performance_table(const std::string& filename);

    // Private data members
    std::string         m_caldb;        //!< Name of or path to the calibration database
    std::string         m_rspname;      //!< Name of the instrument response
    std::string         m_rmffile;      //!< Name of RMF file
    std::string         m_psffile;      //!< Name of PSF file
    GNodeArray          m_psf_logE;     //!< log(E) nodes for PSF interpolation
    std::vector<double> m_logE;         //!< log(E) = log10(E/TeV) - bin centre
    std::vector<double> m_r68;          //!< 68% containment radius of PSF post cuts in degrees
    std::vector<double> m_r80;          //!< 80% containment radius of PSF post cuts in degrees
    double              m_eps;          //!< Integration precision
    double              m_offset_sigma; //!< Sigma for offset angle computation (0=none)

    // New PSF handling
    int                 m_psf_version;  //!< PSF version
    GCTAResponseTable   m_psf_table;    //!< PSF response table

    // New PSF factorisation
    GCTAAeff*           m_aeff;         //!< Effective area component
};

#endif /* GCTARESPONSE_HPP */
