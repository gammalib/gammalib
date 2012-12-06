/***************************************************************************
 *               GLATResponse.hpp  -  Fermi/LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @brief Fermi/LAT Response class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATRESPONSE_HPP
#define GLATRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"
#include "GLATEdisp.hpp"
#include "GLATMeanPsf.hpp"
#include "GEvent.hpp"
#include "GModel.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Interface for the Fermi LAT instrument response function.
 ***************************************************************************/
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Operators
    GLATResponse& operator= (const GLATResponse & rsp);

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GLATResponse* clone(void) const;
    virtual bool          hasedisp(void) const { return false; }
    virtual bool          hastdisp(void) const { return false; }
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;
    virtual std::string   print(void) const;

    // Implemented virtual methods
    virtual double irf(const GEvent&       event,
                       const GSource&      source,
                       const GObservation& obs) const;

    // Other Methods
    void        caldb(const std::string& caldb);
    std::string caldb(void) const { return m_caldb; }
    std::string rspname(void) const { return m_rspname; }
    void        load(const std::string& rspname);
    int         size(void) const { return m_aeff.size(); }
    GLATAeff*   aeff(const int& index) const;
    GLATPsf*    psf(const int& index) const;
    GLATEdisp*  edisp(const int& index) const;
    void        save(const std::string& rspname) const;
    bool        force_mean(void) { return m_force_mean; }
    void        force_mean(const bool& value) { m_force_mean=value; }

    // Reponse methods
    double irf(const GLATEventAtom& event,
               const GSource&       source,
               const GObservation&  obs) const;
    double irf(const GLATEventBin& event,
               const GSource&      source,
               const GObservation& obs) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GLATResponse& rsp);
    void free_members(void);

    // Private members
    std::string               m_caldb;      //!< Name of or path to the calibration database
    std::string               m_rspname;    //!< Name of the instrument response
    bool                      m_hasfront;   //!< Front IRF loaded?
    bool                      m_hasback;    //!< Back IRF loaded?
    bool                      m_force_mean; //!< Use mean PSF in any case
    std::vector<GLATAeff*>    m_aeff;       //!< Effective areas
    std::vector<GLATPsf*>     m_psf;        //!< Point spread functions
    std::vector<GLATEdisp*>   m_edisp;      //!< Energy dispersions
    std::vector<GLATMeanPsf*> m_ptsrc;      //!< Mean PSFs for point sources
};

#endif /* GLATRESPONSE_HPP */
