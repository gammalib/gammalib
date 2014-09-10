/***************************************************************************
 *        GMWLObservation.hpp - Multi-wavelength observation class         *
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
 * @file GMWLObservation.hpp
 * @brief Multi-wavelength observation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMWLOBSERVATION_HPP
#define GMWLOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GMWLResponse.hpp"


/***********************************************************************//**
 * @class GMWLObservation
 *
 * @brief Interface class for multi-wavelength observations
 *
 * This class implements a multi-wavelength observation. A multi-wavelength
 * observation contains spectral points obtained with an unspecified
 * instrument. The spectral points are given in physical units.
 ***************************************************************************/
class GMWLObservation : public GObservation {

public:
    // Constructors and destructors
    GMWLObservation(void);
    explicit GMWLObservation(const std::string& filename);
    explicit GMWLObservation(const std::string& filename, const int& extno);
    explicit GMWLObservation(const std::string& filename, const std::string& extname);
    GMWLObservation(const GMWLObservation& obs);
    virtual ~GMWLObservation(void);

    // Operators
    virtual GMWLObservation& operator= (const GMWLObservation& obs);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GMWLObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GMWLResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void               load(const std::string& filename);
    void               load(const std::string& filename, const int& extno);
    void               load(const std::string& filename, const std::string& extname);
    const std::string& filename(void) const;
    const std::string& extno(void) const;
    const std::string& extname(void) const;
    void               filename(const std::string& filename);
    void               extno(const std::string& extno);
    void               extname(const std::string& extname);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLObservation& obs);
    void free_members(void);

    // Protected members
    std::string  m_instrument;   //!< Instrument name
    std::string  m_filename;     //!< Filename
    std::string  m_extno;        //!< Extension number
    std::string  m_extname;      //!< Extension name
    GMWLResponse m_response;     //!< Response function
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GMWLObservation").
 ***************************************************************************/
inline
std::string GMWLObservation::classname(void) const
{
    return ("GMWLObservation");
}


/***********************************************************************//**
 * @brief Return response
 *
 * @return Response.
 ***************************************************************************/
inline
const GMWLResponse* GMWLObservation::response(void) const
{
    return &m_response;
}


/***********************************************************************//**
 * @brief Return instrument name
 *
 * @return Instrument name.
 ***************************************************************************/
inline
std::string GMWLObservation::instrument(void) const
{
    return "MWL";
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime (always 1).
 ***************************************************************************/
inline
double GMWLObservation::ontime(void) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime (always 1).
 ***************************************************************************/
inline
double GMWLObservation::livetime(void) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @return Deadtime correction factor (always 1).
 ***************************************************************************/
inline
double GMWLObservation::deadc(const GTime& time) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Filename.
 ***************************************************************************/
inline
const std::string& GMWLObservation::filename(void) const
{
    return m_filename;
}


/***********************************************************************//**
 * @brief Return extension number
 *
 * @return Extension number.
 ***************************************************************************/
inline
const std::string& GMWLObservation::extno(void) const
{
    return m_extno;
}


/***********************************************************************//**
 * @brief Return extension name
 *
 * @return Extension name.
 ***************************************************************************/
inline
const std::string& GMWLObservation::extname(void) const
{
    return m_extname;
}


/***********************************************************************//**
 * @brief Set filename
 *
 * @param[in] filename Filename.
 ***************************************************************************/
inline
void GMWLObservation::filename(const std::string& filename)
{
    m_filename = filename;
    return;
}


/***********************************************************************//**
 * @brief Set extension number
 *
 * @param[in] extno Extension number.
 ***************************************************************************/
inline
void GMWLObservation::extno(const std::string& extno)
{
    m_extno = extno;
    return;
}


/***********************************************************************//**
 * @brief Set extension name
 *
 * @param[in] extname Extension name..
 ***************************************************************************/
inline
void GMWLObservation::extname(const std::string& extname)
{
    m_extname = extname;
    return;
}

#endif /* GMWLOBSERVATION_HPP */
