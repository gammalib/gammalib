/***************************************************************************
 *     GCTAResponseCube.hpp - CTA cube analysis response function class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseCube.hpp
 * @brief CTA cube-style response function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSECUBE_HPP
#define GCTARESPONSECUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCTAResponse.hpp"
#include "GCTACubeExposure.hpp"
#include "GCTACubePsf.hpp"
#include "GCTACubeBackground.hpp"
#include "GCTACubeSource.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GPhoton;
class GEvent;
class GObservation;
class GCTAObservation;
class GCTAInstDir;
class GModelSpatialRadial;
class GModelSpatialElliptical;


/***********************************************************************//**
 * @class GCTAResponseCube
 *
 * @brief CTA cube-style response function class
 ***************************************************************************/
class GCTAResponseCube : public GCTAResponse {

public:
    // Constructors and destructors
    GCTAResponseCube(void);
    GCTAResponseCube(const GCTAResponseCube& rsp);
    explicit GCTAResponseCube(const GXmlElement& xml);
    GCTAResponseCube(const GCTACubeExposure&   exposure,
                     const GCTACubePsf&        psf,
                     const GCTACubeBackground& background);
    virtual ~GCTAResponseCube(void);

    // Operators
    virtual GCTAResponseCube& operator=(const GCTAResponseCube & rsp);

    // Implement pure virtual base class methods
    virtual void              clear(void);
    virtual GCTAResponseCube* clone(void) const;
    virtual std::string       classname(void) const;
    virtual bool              is_valid(void) const;
    virtual bool              use_edisp(void) const;
    virtual bool              use_tdisp(void) const;
    virtual bool              apply_edisp(void) const;
    virtual void              apply_edisp(const bool& apply_edisp) const;
    virtual double            irf(const GEvent&       event,
                                  const GPhoton&      photon,
                                  const GObservation& obs) const;
    virtual double            npred(const GPhoton&      photon,
                                    const GObservation& obs) const;
    virtual void              read(const GXmlElement& xml);
    virtual void              write(GXmlElement& xml) const;
    virtual std::string       print(const GChatter& chatter = NORMAL) const;

    // New pure virtual methods
    virtual double   convolve(const GModelSky&    model,
                              const GEvent&       event,
                              const GObservation& obs) const;
    virtual double   nroi(const GModelSky&    model,
                          const GEnergy&      obsEng,
                          const GTime&        obsTime,
                          const GObservation& obs) const;

    // Overload base class methods
    virtual double irf(const GEvent&       event,
                       const GSource&      source,
                       const GObservation& obs) const;
    virtual double irf_ptsrc(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const;
    virtual double irf_radial(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const;
    virtual double irf_elliptical(const GEvent&       event,
                                  const GSource&      source,
                                  const GObservation& obs) const;
    virtual double irf_diffuse(const GEvent&       event,
                               const GSource&      source,
                               const GObservation& obs) const;

    // Other Methods
    const GCTACubeExposure&   exposure(void) const;
    void                      exposure(const GCTACubeExposure& exposure);
    const GCTACubePsf&        psf(void) const;
    void                      psf(const GCTACubePsf& psf);
    const GCTACubeBackground& background(void) const;
    void                      background(const GCTACubeBackground& background);

private:
    // Private methods
    void   init_members(void);
    void   copy_members(const GCTAResponseCube& rsp);
    void   free_members(void);
    int    cache_index(const std::string& name) const;
    double psf_radial(const GModelSpatialRadial* model,
                      const double&              rho_obs,
                      const GSkyDir&             obsDir,
                      const GEnergy&             srcEng,
                      const GTime&               srcTime) const;
    double psf_elliptical(const GModelSpatialElliptical* model,
                          const double&                  rho_obs,
                          const double&                  posangle_obs,
                          const GSkyDir&                 obsDir,
                          const GEnergy&                 srcEng,
                          const GTime&                   srcTime) const;

    // Private data members
    GCTACubeExposure   m_exposure;    //!< Exposure cube
    GCTACubePsf        m_psf;         //!< Mean point spread function
    GCTACubeBackground m_background;  //!< Background cube
    mutable bool       m_apply_edisp; //!< Apply energy dispersion

    // Response cache
    mutable std::vector<GCTACubeSource*> m_cache; //!< Response cache
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAResponseCube").
 ***************************************************************************/
inline
std::string GCTAResponseCube::classname(void) const
{
    return ("GCTAResponseCube");
}


/***********************************************************************//**
 * @brief Signal if response is valid
 *
 * @return True if response is valid
 *
 * @todo: To be implemented (check if GCTACubeExposure and GCTACubePsf is loaded)
 ***************************************************************************/
inline
bool GCTAResponseCube::is_valid(void) const
{
    return (true);
}


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
bool GCTAResponseCube::use_edisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if time dispersion will be used
 *
 * @return False.
 ***************************************************************************/
inline
bool GCTAResponseCube::use_tdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if energy dispersion should be applied
 *
 * @return True if energy dispersion should be applied
 ***************************************************************************/
inline
bool GCTAResponseCube::apply_edisp(void) const
{
    return m_apply_edisp;
}

/***********************************************************************//**
 * @brief Signal if energy dispersion should be applied
 *
 * @param[in] apply_edisp Set true if energy dispersion should be applied
 ***************************************************************************/
inline
void GCTAResponseCube::apply_edisp(const bool& apply_edisp) const
{
    m_apply_edisp = apply_edisp;
    return;
}


/***********************************************************************//**
 * @brief Return exposure cube
 *
 * @return Reference to exposure cube.
 ***************************************************************************/
inline
const GCTACubeExposure& GCTAResponseCube::exposure(void) const
{
    return (m_exposure);
}


/***********************************************************************//**
 * @brief Set exposure cube
 *
 * @param[in] exposure Exposure cube.
 ***************************************************************************/
inline
void GCTAResponseCube::exposure(const GCTACubeExposure& exposure)
{
    m_exposure = exposure;
    return;
}


/***********************************************************************//**
 * @brief Return cube analysis point spread function
 *
 * @return Reference to cube analysis point spread function.
 ***************************************************************************/
inline
const GCTACubePsf& GCTAResponseCube::psf(void) const
{
    return (m_psf);
}


/***********************************************************************//**
 * @brief Set cube analysis point spread function cube
 *
 * @param[in] psf Cube analysis point spread function.
 ***************************************************************************/
inline
void GCTAResponseCube::psf(const GCTACubePsf& psf)
{
    m_psf = psf;
    return;
}


/***********************************************************************//**
 * @brief Set cube analysis background cube
 *
 * @param[in] background Cube analysis background cube.
 ***************************************************************************/
inline
void GCTAResponseCube::background(const GCTACubeBackground& background)
{
    m_background = background;
    return;
}


/***********************************************************************//**
 * @brief Return cube analysis background cube
 *
 * @return Reference to cube analysis background cube.
 ***************************************************************************/
inline
const GCTACubeBackground& GCTAResponseCube::background(void) const
{
    return (m_background);
}

#endif /* GCTARESPONSECUBE_HPP */
