/***************************************************************************
 *             GLATObservation.hpp - Fermi/LAT observation class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GLATObservation.hpp
 * @brief Fermi/LAT observation class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATOBSERVATION_HPP
#define GLATOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GObservation.hpp"
#include "GLATResponse.hpp"
#include "GLATLtCube.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Fermi/LAT observation class
 *
 * This class implements a Fermi/LAT observation.
 ***************************************************************************/
class GLATObservation : public GObservation {

public:
    // Constructors and destructors
    GLATObservation(void);
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation(void);

    // Operators
    GLATObservation& operator=(const GLATObservation& obs);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GLATObservation*    clone(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GLATResponse& response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void              load_unbinned(const std::string& ft1name,
                                    const std::string& ft2name,
                                    const std::string& ltcube_name);
    void              load_binned(const std::string& cntmap_name,
                                  const std::string& expmap_name,
                                  const std::string& ltcube_name);
    void              response(const std::string& irfname,
                               const std::string& caldb = "");    
    const GLATLtCube* ltcube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATObservation& obs);
    void free_members(void);

    // Protected members
    std::string  m_ft1file;      //!< FT1 filename
    std::string  m_ft2file;      //!< FT2 filename
    std::string  m_ltfile;       //!< Lifetime cube filename
    std::string  m_cntfile;      //!< Counts map filename
    std::string  m_expfile;      //!< Exposure map filename
    GLATResponse m_response;     //!< Instrument response functions
    GLATLtCube*  m_ltcube;       //!< Pointer to lifetime cube
};


/***********************************************************************//**
 * @brief Return Fermi/LAT response function
 *
 * @return Fermi/LAT response function
 ***************************************************************************/
inline
const GLATResponse& GLATObservation::response(void) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Return Fermi/LAT livetime cube
 *
 * @return Fermi/LAT livetime cube
 ***************************************************************************/
inline
const GLATLtCube* GLATObservation::ltcube(void) const
{
    // Return livetime cube pointer
    return m_ltcube;
}


/***********************************************************************//**
 * @brief Return instrument name
 *
 * @return Instrument name ("LAT")
 ***************************************************************************/
inline
std::string GLATObservation::instrument(void) const
{
    // Return instrument name
    return ("LAT");
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime
 ***************************************************************************/
inline
double GLATObservation::ontime(void) const
{
    // Return ontime
    return 0.0;
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime
 ***************************************************************************/
inline
double GLATObservation::livetime(void) const
{
    // Return livetime
    return 0.0;
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @return Deadtime correction factor
 ***************************************************************************/
inline
double GLATObservation::deadc(const GTime& time) const
{
    // Return deadtime correction factor
    return 0.0;
}

#endif /* GLATOBSERVATION_HPP */
