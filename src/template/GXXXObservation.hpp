/***************************************************************************
 *           GXXXObservation.hpp - [INSTRUMENT] observation class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXObservation.hpp
 * @brief [INSTRUMENT] observation class definition
 * @author [AUTHOR]
 */

#ifndef GXXXOBSERVATION_HPP
#define GXXXOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GObservation.hpp"
#include "GXXXResponse.hpp"

/* __ Forward declarations _______________________________________________ */
class GResponse;
class GXmlElement;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GXXXObservation
 *
 * @brief [INSTRUMENT] observation class
 *
 * The [INSTRUMENT] observation class defines an observation.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GXXXObservation : public GObservation {

public:
    // Constructors and destructors
    GXXXObservation(void);
    GXXXObservation(const GXXXObservation& obs);
    virtual ~GXXXObservation(void);

    // Operators
    virtual GXXXObservation& operator=(const GXXXObservation& obs);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GXXXObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GXXXResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void ontime(const double& ontime);
    void livetime(const double& livetime);
    void deadc(const double& deadc);
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXObservation& obs);
    void free_members(void);

    // Protected members
    GXXXResponse m_response; //!< Response functions
    double       m_ontime;   //!< Ontime (sec)
    double       m_livetime; //!< Livetime (sec)
    double       m_deadc;    //!< Deadtime correction
    // TODO: Add any data members that are necessary. Note that the event
    // list or cubes (type GEvents) are stored in the GObservation base
    // class
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXObservation").
 ***************************************************************************/
inline
std::string GXXXObservation::classname(void) const
{
    return ("GXXXObservation");
}


/***********************************************************************//**
 * @brief Return pointer to response function
 *
 * @return Response function pointer.
 ***************************************************************************/
inline
const GXXXResponse* GXXXObservation::response(void) const
{
    // Return response pointer
    return &m_response;
}


/***********************************************************************//**
 * @brief Return instrument
 *
 * @return Instrument name.
 ***************************************************************************/
inline
std::string GXXXObservation::instrument(void) const
{
    // Return instrument
    return ("XXX");
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime (seconds).
 ***************************************************************************/
inline
double GXXXObservation::ontime(void) const
{
    // Return ontime
    return (m_ontime);
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime (seconds).
 ***************************************************************************/
inline
double GXXXObservation::livetime(void) const
{
    // Return livetime
    return (m_livetime);
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @param[in] time Time.
 *
 * @return Deadtime correction factor.
 ***************************************************************************/
inline
double GXXXObservation::deadc(const GTime& time) const
{
    // Return livetime
    return (m_deadc);
}


/***********************************************************************//**
 * @brief Set ontime
 *
 * @param[in] ontime Ontime.
 ***************************************************************************/
inline
void GXXXObservation::ontime(const double& ontime)
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
void GXXXObservation::livetime(const double& livetime)
{
    m_livetime = livetime;
    return;
}


/***********************************************************************//**
 * @brief Set deadtime correction factor
 *
 * @param[in] deadc Deadtime correction factor.
 ***************************************************************************/
inline
void GXXXObservation::deadc(const double& deadc)
{
    m_deadc = deadc;
    return;
}

#endif /* GXXXOBSERVATION_HPP */
