/***************************************************************************
 *                GCTAObservation.hpp - CTA Observation class              *
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
 * @file GCTAObservation.hpp
 * @brief CTA observation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAOBSERVATION_HPP
#define GCTAOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GObservation.hpp"
#include "GCTAResponse.hpp"
#include "GCTAPointing.hpp"

/* __ Forward declarations _______________________________________________ */
class GTime;
class GFits;
class GFitsHDU;
class GResponse;
class GXmlElement;
class GCaldb;


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief CTA observation class
 *
 * This class implements a CTA observation.
 ***************************************************************************/
class GCTAObservation : public GObservation {

    // Friends
    friend class GCTAModelCubeBackground;

public:
    // Constructors and destructors
    GCTAObservation(void);
    explicit GCTAObservation(const std::string& instrument);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Operators
    GCTAObservation& operator=(const GCTAObservation& obs);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GCTAObservation*    clone(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GCTAResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    bool                hasresponse(void) const;
    void                read(const GFits& fits);
    void                write(GFits& fits) const;
    void                load(const std::string& filename);
    void                save(const std::string& filename,
                             const bool& clobber = false) const;
    void                response(const std::string& rspname,
                                 const GCaldb& caldb);
    void                pointing(const GCTAPointing& pointing);
    const GCTAPointing& pointing(void) const;
    void                obs_id(const int& id);
    const int&          obs_id(void) const;
    void                ra_obj(const double& ra);
    const double&       ra_obj(void) const;
    void                dec_obj(const double& dec);
    const double&       dec_obj(void) const;
    void                ontime(const double& ontime);
    void                livetime(const double& livetime);
    void                deadc(const double& deadc);
    void                eventfile(const std::string& filename);
    const std::string&  eventfile(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAObservation& obs);
    void free_members(void);
    void read_attributes(const GFitsHDU& hdu);
    void write_attributes(GFitsHDU& hdu) const;

    // Protected members
    std::string   m_instrument;  //!< Instrument name
    std::string   m_eventfile;   //!< Event filename
    GCTAResponse* m_response;    //!< Pointer to instrument response functions
    GCTAPointing  m_pointing;    //!< Pointing direction
    int           m_obs_id;      //!< Observation ID
    double        m_ontime;      //!< Ontime (seconds)
    double        m_livetime;    //!< Livetime (seconds)
    double        m_deadc;       //!< Deadtime correction (livetime/ontime)
    double        m_ra_obj;      //!< Right Ascension of object (degrees)
    double        m_dec_obj;     //!< Declination of object (degrees)

    // Special protected member for GCTAModelCubeBackground friend
    std::string  m_bgdfile;      //!< Background filename
};


/***********************************************************************//**
 * @brief Return instrument name
 *
 * @return Instrument name.
 ***************************************************************************/
inline
std::string GCTAObservation::instrument(void) const
{
    return m_instrument;
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime in seconds.
 ***************************************************************************/
inline
double GCTAObservation::ontime(void) const
{
    return m_ontime;
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime in seconds.
 ***************************************************************************/
inline
double GCTAObservation::livetime(void) const
{
    return m_livetime;
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @param[in] time Time.
 * @return Deadtime correction factor.
 *
 * Returns the deadtime correction factor as function of time. The deadtime
 * correction factor is defined by the livetime divided by the ontime.
 ***************************************************************************/
inline
double GCTAObservation::deadc(const GTime& time) const
{
    return m_deadc;
}


/***********************************************************************//**
 * @brief Signal if CTA observation contains response information
 *
 * @return True if CTA observation contains response information.
 ***************************************************************************/
inline
bool GCTAObservation::hasresponse(void) const
{
    return (m_response != NULL);
}


/***********************************************************************//**
 * @brief Return CTA pointing
 *
 * @return CTA pointing
 ***************************************************************************/
inline
const GCTAPointing& GCTAObservation::pointing(void) const
{
    return m_pointing;
}


/***********************************************************************//**
 * @brief Set CTA pointing
 *
 * @param[in] pointing CTA pointing.
 ***************************************************************************/
inline
void GCTAObservation::pointing(const GCTAPointing& pointing)
{
    m_pointing = pointing;
    return;
}


/***********************************************************************//**
 * @brief Set CTA observation identifier
 *
 * @param[in] id Observation identifier.
 ***************************************************************************/
inline
void GCTAObservation::obs_id(const int& id)
{
    m_obs_id = id;
    return;
}


/***********************************************************************//**
 * @brief Return CTA observation identifier
 *
 * @return Observation identifier.
 ***************************************************************************/
inline
const int& GCTAObservation::obs_id(void) const
{
    return m_obs_id;
}


/***********************************************************************//**
 * @brief Set CTA object Right Ascension
 *
 * @param[in] ra Object Right Ascension.
 ***************************************************************************/
inline
void GCTAObservation::ra_obj(const double& ra)
{
    m_ra_obj = ra;
    return;
}


/***********************************************************************//**
 * @brief Return CTA object Right Ascension
 *
 * @return Object Right Ascension.
 ***************************************************************************/
inline
const double& GCTAObservation::ra_obj(void) const
{
    return m_ra_obj;
}


/***********************************************************************//**
 * @brief Set CTA object Declination
 *
 * @param[in] dec Object Declination.
 ***************************************************************************/
inline
void GCTAObservation::dec_obj(const double& dec)
{
    m_dec_obj = dec;
    return;
}


/***********************************************************************//**
 * @brief Return CTA object Declination
 *
 * @return Object Declination.
 ***************************************************************************/
inline
const double& GCTAObservation::dec_obj(void) const
{
    return m_dec_obj;
}


/***********************************************************************//**
 * @brief Set ontime
 *
 * @param[in] ontime Ontime.
 ***************************************************************************/
inline
void GCTAObservation::ontime(const double& ontime)
{
    m_ontime = ontime;
    return;
}


/***********************************************************************//**
 * @brief Set livetime
 *
 * @param[in] livetime Livetime.
 ***************************************************************************/
inline
void GCTAObservation::livetime(const double& livetime)
{
    m_livetime = livetime;
    return;
}


/***********************************************************************//**
 * @brief Set deadtime correction
 *
 * @param[in] deadc Deadtime correction.
 ***************************************************************************/
inline
void GCTAObservation::deadc(const double& deadc)
{
    m_deadc = deadc;
    return;
}


/***********************************************************************//**
 * @brief Set event file name
 *
 * @param[in] filename Event file name.
 ***************************************************************************/
inline
void GCTAObservation::eventfile(const std::string& filename)
{
    m_eventfile = filename;
    return;
}


/***********************************************************************//**
 * @brief Return event file name
 *
 * @return Event file name.
 ***************************************************************************/
inline
const std::string& GCTAObservation::eventfile(void) const
{
    return m_eventfile;
}

#endif /* GCTAOBSERVATION_HPP */
