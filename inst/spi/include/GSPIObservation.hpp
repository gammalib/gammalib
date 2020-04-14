/***************************************************************************
 *           GSPIObservation.hpp - INTEGRAL/SPI observation class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIObservation.hpp
 * @brief INTEGRAL/SPI observation class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIOBSERVATION_HPP
#define GSPIOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFilename.hpp"
#include "GObservation.hpp"
#include "GSPIResponse.hpp"

/* __ Forward declarations _______________________________________________ */
class GTime;
class GFits;
class GResponse;
class GXmlElement;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIObservation
 *
 * @brief INTEGRAL/SPI observation class
 *
 * The INTEGRAL/SPI observation class defines an observation.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GSPIObservation : public GObservation {

public:
    // Constructors and destructors
    GSPIObservation(void);
    explicit GSPIObservation(const GXmlElement& xml);
    explicit GSPIObservation(const GFilename& filename);
    GSPIObservation(const GSPIObservation& obs);
    virtual ~GSPIObservation(void);

    // Operators
    virtual GSPIObservation& operator=(const GSPIObservation& obs);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GSPIObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GSPIResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void read(const GFits& fits);
    void load(const GFilename& filename);
    void ontime(const double& ontime);
    void livetime(const double& livetime);
    void deadc(const double& deadc);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSPIObservation& obs);
    void free_members(void);

    // Protected members
    std::string  m_instrument;   //!< Instrument name
    GFilename    m_filename;     //!< OG FITS filename
    GFilename    m_rsp_grpname;  //!< Response group FITS filename (optional)
    GFilename    m_rsp_filename; //!< Response FITS filename (optional)
    GSPIResponse m_response;     //!< Response functions
    double       m_ontime;       //!< Ontime (sec)
    double       m_livetime;     //!< Livetime (sec)
    double       m_deadc;        //!< Deadtime correction
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIObservation").
 ***************************************************************************/
inline
std::string GSPIObservation::classname(void) const
{
    return ("GSPIObservation");
}


/***********************************************************************//**
 * @brief Return pointer to response function
 *
 * @return Response function pointer.
 ***************************************************************************/
inline
const GSPIResponse* GSPIObservation::response(void) const
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
std::string GSPIObservation::instrument(void) const
{
    // Return instrument
    return (m_instrument);
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime (seconds).
 ***************************************************************************/
inline
double GSPIObservation::ontime(void) const
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
double GSPIObservation::livetime(void) const
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
double GSPIObservation::deadc(const GTime& time) const
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
void GSPIObservation::ontime(const double& ontime)
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
void GSPIObservation::livetime(const double& livetime)
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
void GSPIObservation::deadc(const double& deadc)
{
    m_deadc = deadc;
    return;
}

#endif /* GSPIOBSERVATION_HPP */
