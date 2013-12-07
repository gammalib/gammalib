/***************************************************************************
 *             GLATObservation.hpp - Fermi LAT Observation class           *
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
 * @brief LAT Observation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATOBSERVATION_HPP
#define GLATOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GObservation.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GLATLtCube.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Interface for the LAT observation class
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
    GLATObservation& operator= (const GLATObservation& obs);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GLATObservation*    clone(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GLATResponse& response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const { return 0.0; }
    virtual double              livetime(void) const { return 0.0; }
    virtual double              deadc(const GTime& time) const { return 0.0; }
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GLATPointing& pointing(void) const;
    void                load_unbinned(const std::string& ft1name,
                                      const std::string& ft2name,
                                      const std::string& ltcube_name);
    void                load_binned(const std::string& cntmap_name,
                                    const std::string& expmap_name,
                                    const std::string& ltcube_name);
    void                response(const std::string& irfname,
                                 std::string caldb = "");    
    GLATLtCube*         ltcube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATObservation& obs);
    void free_members(void);

    // Protected members
    std::string   m_ft1file;      //!< FT1 filename
    std::string   m_ft2file;      //!< FT2 filename
    std::string   m_ltfile;       //!< Lifetime cube filename
    std::string   m_cntfile;      //!< Counts map filename
    std::string   m_expfile;      //!< Exposure map filename
    GLATResponse* m_response;     //!< Pointer to instrument response functions
    GLATPointing* m_pointing;     //!< Pointer to pointing direction
    GLATLtCube*   m_ltcube;       //!< Pointer to lifetime cube
};


/***********************************************************************//**
 * @brief Returns pointer to LAT response function
 ***************************************************************************/
inline
const GLATResponse& GLATObservation::response(void) const
{
    // Return response pointer
    return *m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to LAT pointing direction
 ***************************************************************************/
inline
const GLATPointing& GLATObservation::pointing(void) const
{
    // Return response pointer
    return *m_pointing;
}


/***********************************************************************//**
 * @brief Returns pointer to LAT livetime cube
 ***************************************************************************/
inline
GLATLtCube* GLATObservation::ltcube(void) const
{
    // Return livetime cube pointer
    return m_ltcube;
}


/***********************************************************************//**
 * @brief Returns instrument name
 ***************************************************************************/
inline
std::string GLATObservation::instrument(void) const
{
    // Return instument name
    return ("LAT");
}

#endif /* GLATOBSERVATION_HPP */
