/***************************************************************************
 *      GCTACubeResponse.hpp - CTA cube-style response function class      *
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
 * @file GCTACubeResponse.hpp
 * @brief CTA cube-style response function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTACUBERESPONSE_HPP
#define GCTACUBERESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GResponse.hpp"
#include "GCTAExposure.hpp"
#include "GCTAMeanPsf.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GPhoton;
class GEvent;
class GObservation;
class GCTAObservation;
class GCTAInstDir;


/***********************************************************************//**
 * @class GCTACubeResponse
 *
 * @brief CTA cube-style response function class
 ***************************************************************************/
class GCTACubeResponse : public GResponse {

public:
    // Constructors and destructors
    GCTACubeResponse(void);
    GCTACubeResponse(const GCTACubeResponse& rsp);
    virtual ~GCTACubeResponse(void);

    // Operators
    virtual GCTACubeResponse& operator=(const GCTACubeResponse & rsp);

    // Implement pure virtual base class methods
    virtual void              clear(void);
    virtual GCTACubeResponse* clone(void) const;
    virtual bool              use_edisp(void) const;
    virtual bool              use_tdisp(void) const;
    virtual bool              apply_edisp(void) const;
    virtual void              apply_edisp(const bool& apply_edisp) const;
    virtual double            irf(const GEvent&       event,
                                  const GPhoton&      photon,
                                  const GObservation& obs) const;
    virtual double            npred(const GPhoton&      photon,
                                    const GObservation& obs) const;
    virtual std::string       print(const GChatter& chatter = NORMAL) const;

    // Other Methods
    const GCTAExposure&       exposure(void) const;
    void                      exposure(const GCTAExposure& exposure);
    const GCTAMeanPsf&        psf(void) const;
    void                      psf(const GCTAMeanPsf& psf);

private:
    // Private methods
    void                   init_members(void);
    void                   copy_members(const GCTACubeResponse& rsp);
    void                   free_members(void);
    const GCTAObservation& retrieve_obs(const std::string& origin,
                                        const GObservation& obs) const;
    const GCTAInstDir&     retrieve_dir(const std::string& origin,
                                        const GEvent&      event) const;

    // Private data members
    GCTAExposure m_exposure;    //!< Exposure cube
    GCTAMeanPsf  m_psf;         //!< Mean point spread function
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
bool GCTACubeResponse::use_edisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if energy dispersion should be applied
 *
 * @return True if energy dispersion should be applied
 ***************************************************************************/
inline
bool GCTACubeResponse::apply_edisp(void) const
{
    return false;
}

/***********************************************************************//**
 * @brief Signal if energy dispersion should be applied
 *
 * @param[in] apply_edisp Set true if energy dispersion should be applied
 ***************************************************************************/
inline
void GCTACubeResponse::apply_edisp(const bool& apply_edisp) const
{
    return;
}


/***********************************************************************//**
 * @brief Signal if time dispersion will be used
 *
 * @return False.
 ***************************************************************************/
inline
bool GCTACubeResponse::use_tdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Return exposure cube
 *
 * @return Reference to exposure cube.
 ***************************************************************************/
inline
const GCTAExposure& GCTACubeResponse::exposure(void) const
{
    return m_exposure;
}


/***********************************************************************//**
 * @brief Set exposure cube
 *
 * @param[in] exposure Exposure cube.
 ***************************************************************************/
inline
void GCTACubeResponse::exposure(const GCTAExposure& exposure)
{
    m_exposure = exposure;
    return;
}


/***********************************************************************//**
 * @brief Return point spread function cube
 *
 * @return Reference to point spread function cube.
 ***************************************************************************/
inline
const GCTAMeanPsf& GCTACubeResponse::psf(void) const
{
    return m_psf;
}


/***********************************************************************//**
 * @brief Set pointer to point spread function cube
 *
 * @param[in] psf Pointer to point spread function cube.
 ***************************************************************************/
inline
void GCTACubeResponse::psf(const GCTAMeanPsf& psf)
{
    m_psf = psf;
    return;
}

#endif /* GCTACUBERESPONSE_HPP */
