/***************************************************************************
 *                GLATResponse.hpp - Fermi LAT Response class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
 * @file GLATResponse.hpp
 * @brief Fermi LAT Response class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATRESPONSE_HPP
#define GLATRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GResponse.hpp"
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"
#include "GLATEdisp.hpp"
#include "GLATMeanPsf.hpp"

/* __ Forward declarations _______________________________________________ */
class GSource;
class GEvent;
class GPhoton;
class GObservation;
class GEnergy;
class GTime;
class GModelSky;
class GLATEventAtom;
class GLATEventBin;


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Fermi/LAT Response class
 ***************************************************************************/
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Operators
    GLATResponse& operator=(const GLATResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GLATResponse* clone(void) const;
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

    // Overloaded virtual base class methods
    virtual double        irf_spatial(const GEvent&       event,
                                      const GSource&      source,
                                      const GObservation& obs) const;

    // Other Methods
    int                size(void) const;
    const std::string& rspname(void) const;
    void               load(const std::string& rspname);
    void               save(const std::string& rspname) const;
    const bool&        force_mean(void) const;
    void               force_mean(const bool& value);
    GLATAeff*          aeff(const int& index) const;
    GLATPsf*           psf(const int& index) const;
    GLATEdisp*         edisp(const int& index) const;

    // Reponse methods
    double irf_spatial_atom(const GLATEventAtom& event,
                            const GSource&       source,
                            const GObservation&  obs) const;
    double irf_spatial_bin(const GLATEventBin& event,
                           const GSource&      source,
                           const GObservation& obs) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GLATResponse& rsp);
    void free_members(void);

    // Private members
    std::string               m_rspname;    //!< Name of the instrument response
    bool                      m_force_mean; //!< Use mean PSF in any case
    std::vector<GLATAeff*>    m_aeff;       //!< Effective areas
    std::vector<GLATPsf*>     m_psf;        //!< Point spread functions
    std::vector<GLATEdisp*>   m_edisp;      //!< Energy dispersions
    std::vector<GLATMeanPsf*> m_ptsrc;      //!< Mean PSFs for point sources
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATResponse").
 ***************************************************************************/
inline
std::string GLATResponse::classname(void) const
{
    return ("GLATResponse");
}


/***********************************************************************//**
 * @brief Signal if response uses energy dispersion
 *
 * @return True if response uses energy dispersion.
 ***************************************************************************/
inline
bool GLATResponse::use_edisp(void) const
{
    return (false);
}


/***********************************************************************//**
 * @brief Signal if response uses time dispersion
 *
 * @return True if response uses time dispersion.
 ***************************************************************************/
inline
bool GLATResponse::use_tdisp(void) const
{
    return (false);
}


/***********************************************************************//**
 * @brief Return number of event types
 *
 * @return Number of event types.
 *
 * Returns the number of events types in this reponse.
 ***************************************************************************/
inline
int GLATResponse::size(void) const
{
    return (int)m_aeff.size();
}


/***********************************************************************//**
 * @brief Return response name
 *
 * @return Response name.
 ***************************************************************************/
inline
const std::string& GLATResponse::rspname(void) const
{
    return (m_rspname);
}


/***********************************************************************//**
 * @brief Signal if mean PSF should be used for response computation
 *
 * @return True if mean PSF should be used for response computation.
 ***************************************************************************/
inline
const bool& GLATResponse::force_mean(void) const
{
    return (m_force_mean);
}


/***********************************************************************//**
 * @brief Set if mean PSF should be used for response computation
 *
 * @param[in] value True if mean PSF should be used for response computation.
 ***************************************************************************/
inline
void GLATResponse::force_mean(const bool& value)
{
    m_force_mean = value;
    return;
}

#endif /* GLATRESPONSE_HPP */
