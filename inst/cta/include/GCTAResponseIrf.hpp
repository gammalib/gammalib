/***************************************************************************
 *        GCTAResponseIrf.hpp - CTA instrument response function class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseIrf.hpp
 * @brief CTA instrument response function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSEIRF_HPP
#define GCTARESPONSEIRF_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAResponse.hpp"
#include "GCaldb.hpp"

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
class GCTABackground;


/***********************************************************************//**
 * @class GCTAResponseIrf
 *
 * @brief CTA instrument response function class
 ***************************************************************************/
class GCTAResponseIrf : public GCTAResponse {

public:
    // Constructors and destructors
    GCTAResponseIrf(void);
    GCTAResponseIrf(const GCTAResponseIrf& rsp);
    GCTAResponseIrf(const std::string& rspname, const GCaldb& caldb);
    virtual ~GCTAResponseIrf(void);

    // Operators
    virtual GCTAResponseIrf& operator=(const GCTAResponseIrf & rsp);

    // Implement pure virtual base class methods
    virtual void             clear(void);
    virtual GCTAResponseIrf* clone(void) const;
    virtual bool             use_edisp(void) const;
    virtual bool             use_tdisp(void) const;
    virtual double           irf(const GEvent&       event,
                                 const GPhoton&      photon,
                                 const GObservation& obs) const;
    virtual double           npred(const GPhoton&      photon,
                                   const GObservation& obs) const;
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;
    virtual std::string      print(const GChatter& chatter = NORMAL) const;

    // Overload virtual base class methods
    virtual double   irf_radial(const GEvent&       event,
                                const GSource&      source,
                                const GObservation& obs) const;
    virtual double   irf_elliptical(const GEvent&       event,
                                    const GSource&      source,
                                    const GObservation& obs) const;
    virtual double   irf_diffuse(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const;
    virtual double   npred_radial(const GSource&      source,
                                  const GObservation& obs) const;
    virtual double   npred_elliptical(const GSource&      source,
                                      const GObservation& obs) const;
    virtual double   npred_diffuse(const GSource&      source,
                                   const GObservation& obs) const;
    virtual GEbounds ebounds_src(const GEnergy& obsEnergy) const;

    // Other Methods
    bool                  apply_edisp(void) const;
    void                  apply_edisp(const bool& apply_edisp) const;
    GCTAEventAtom*        mc(const double& area, const GPhoton& photon,
                             const GObservation& obs, GRan& ran) const;
    void                  caldb(const GCaldb& caldb);
    const GCaldb&         caldb(void) const;
    void                  load(const std::string& rspname);
    void                  eps(const double& eps);
    const double&         eps(void) const;
    void                  load_aeff(const std::string& filename);
    void                  load_psf(const std::string& filename);
    void                  load_edisp(const std::string& filename);
    void                  load_background(const std::string& filename);
    void                  offset_sigma(const double& sigma);
    double                offset_sigma(void) const;
    const GCTAAeff*       aeff(void) const;
    void                  aeff(GCTAAeff* aeff);
    const GCTAPsf*        psf(void) const;
    void                  psf(GCTAPsf* psf);
    const GCTAEdisp*      edisp(void) const;
    void                  edisp(GCTAEdisp* edisp);
    const GCTABackground* background(void) const;
    void                  background(GCTABackground* background);

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
    double edisp(const GEnergy& obsEng,
                 const double&  theta,
                 const double&  phi,
                 const double&  zenith,
                 const double&  azimuth,
                 const double&  srcLogEng) const;
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
    void        init_members(void);
    void        copy_members(const GCTAResponseIrf& rsp);
    void        free_members(void);
    std::string irf_filename(const std::string& filename) const;

    // Private data members
    GCaldb          m_caldb;          //!< Calibration database
    std::string     m_rspname;        //!< Name of the instrument response
    double          m_eps;            //!< Integration precision
    GCTAAeff*       m_aeff;           //!< Effective area
    GCTAPsf*        m_psf;            //!< Point spread function
    GCTAEdisp*      m_edisp;          //!< Energy dispersion
    GCTABackground* m_background;     //!< Energy dispersion
    mutable bool    m_apply_edisp;    //!< Apply energy dispersion

    // XML response filename
    std::string     m_xml_caldb;      //!< Calibration database string in XML file
    std::string     m_xml_rspname;    //!< Response name in XML file
    std::string     m_xml_aeff;       //!< Aeff file name in XML file
    std::string     m_xml_psf;        //!< PSF file name in XML file
    std::string     m_xml_edisp;      //!< Edisp file name in XML file
    std::string     m_xml_background; //!< Background file name in XML file

    // Npred cache
    mutable std::vector<std::string> m_npred_names;    //!< Model names
    mutable std::vector<GEnergy>     m_npred_energies; //!< Model energy
    mutable std::vector<GTime>       m_npred_times;    //!< Model time
    mutable std::vector<double>      m_npred_values;   //!< Model values
};


/***********************************************************************//**
 * @brief Signal if response uses energy dispersion
 *
 * @return True if response uses energy dispersion
 *
 * Signals if the response uses energy dispersion. This implies that the
 * apply_edisp flag has been set to true and that energy dispersion response
 * information is available.
 ***************************************************************************/
inline
bool GCTAResponseIrf::use_edisp(void) const
{
    return (m_apply_edisp && (m_edisp != NULL));
}


/***********************************************************************//**
 * @brief Signal if energy dispersion should be applied
 *
 * @return True if energy dispersion should be applied
 ***************************************************************************/
inline
bool GCTAResponseIrf::apply_edisp(void) const
{
    return m_apply_edisp;
}

/***********************************************************************//**
 * @brief Signal if energy dispersion should be applied
 *
 * @param[in] apply_edisp Set true if energy dispersion should be applied
 ***************************************************************************/
inline
void GCTAResponseIrf::apply_edisp(const bool& apply_edisp) const
{
    m_apply_edisp = apply_edisp;
    return;
}


/***********************************************************************//**
 * @brief Signal if time dispersion will be used
 *
 * @return False.
 ***************************************************************************/
inline
bool GCTAResponseIrf::use_tdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Return calibration database
 *
 * @return Calibration database.
 ***************************************************************************/
inline
const GCaldb& GCTAResponseIrf::caldb(void) const
{
    return m_caldb;
}


/***********************************************************************//**
 * @brief Set path to the calibration database
 *
 * @param[in] caldb Calibration database.
 *
 * Sets the calibration database for the CTA response.
 ***************************************************************************/
inline
void GCTAResponseIrf::caldb(const GCaldb& caldb)
{
    m_caldb = caldb;
    return;
}


/***********************************************************************//**
 * @brief Set computation precision
 *
 * @param[in] eps Computation precision.
 ***************************************************************************/
inline
void GCTAResponseIrf::eps(const double& eps)
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
const double& GCTAResponseIrf::eps(void) const
{
    return m_eps;
}


/***********************************************************************//**
 * @brief Return pointer to effective area response
 *
 * @return Pointer to effective area response.
 ***************************************************************************/
inline
const GCTAAeff* GCTAResponseIrf::aeff(void) const
{
    return m_aeff;
}


/***********************************************************************//**
 * @brief Set pointer to effective area response
 *
 * @param[in] aeff Pointer to effective area response.
 ***************************************************************************/
inline
void GCTAResponseIrf::aeff(GCTAAeff* aeff)
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
const GCTAPsf* GCTAResponseIrf::psf(void) const
{
    return m_psf;
}


/***********************************************************************//**
 * @brief Set pointer to point spread function
 *
 * @param[in] psf Pointer to point spread function.
 ***************************************************************************/
inline
void GCTAResponseIrf::psf(GCTAPsf* psf)
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
const GCTAEdisp* GCTAResponseIrf::edisp(void) const
{
    return m_edisp;
}


/***********************************************************************//**
 * @brief Set pointer to energy dispersion
 *
 * @param[in] edisp Pointer to energy dispersion.
 ***************************************************************************/
inline
void GCTAResponseIrf::edisp(GCTAEdisp* edisp)
{
    m_edisp = edisp;
    return;
}




/***********************************************************************//**
 * @brief Return pointer to background model
 *
 * @return Pointer to background model.
 ***************************************************************************/
inline
const GCTABackground* GCTAResponseIrf::background(void) const
{
    return m_background;
}


/***********************************************************************//**
 * @brief Set pointer to background model
 *
 * @param[in] background Pointer to background model.
 ***************************************************************************/
inline
void GCTAResponseIrf::background(GCTABackground* background)
{
    m_background = background;
    return;
}

#endif /* GCTARESPONSEIRF_HPP */
