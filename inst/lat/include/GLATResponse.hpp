/***************************************************************************
 *               GLATResponse.hpp  -  Fermi LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATResponse.hpp
 * @brief Fermi LAT Response class definition
 * @author J. Knodlseder
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
#include "GEnergy.hpp"
#include "GTime.hpp"


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

    // Implement pure virtual base class methods
    void          clear(void);
    GLATResponse* clone(void) const;
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    double        irf(const GEvent& event, const GModel& model,
                      const GEnergy& srcEng, const GTime& srcTime,
                      const GObservation& obs) const;
    double        npred(const GModel& model, const GEnergy& srcEng,
                        const GTime& srcTime,
                        const GObservation& obs) const;
    std::string   print(void) const;

    // Other Methods
    void        caldb(const std::string& caldb);
    std::string caldb(void) const { return m_caldb; }
    void        load(const std::string& rspname);
    int         size(void) const { return m_aeff.size(); }
    GLATAeff*   aeff(const int& index) const;
    GLATPsf*    psf(const int& index) const;
    GLATEdisp*  edisp(const int& index) const;
    void        save(const std::string& rspname) const;

    // Reponse methods
    double irf(const GLATEventAtom& event, const GModel& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GObservation& obs) const;
    double irf(const GLATEventBin& event, const GModel& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GObservation& obs) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GLATResponse& rsp);
    void free_members(void);

    // Private members
    std::string               m_caldb;    //!< Name of or path to the calibration database
    std::string               m_rspname;  //!< Name of the instrument response
    bool                      m_hasfront; //!< Front IRF loaded?
    bool                      m_hasback;  //!< Back IRF loaded?
    std::vector<GLATAeff*>    m_aeff;     //!< Effective areas
    std::vector<GLATPsf*>     m_psf;      //!< Point spread functions
    std::vector<GLATEdisp*>   m_edisp;    //!< Energy dispersions
    std::vector<GLATMeanPsf*> m_ptsrc;    //!< Mean PSFs for point sources
};

#endif /* GLATRESPONSE_HPP */
