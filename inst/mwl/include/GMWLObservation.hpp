/***************************************************************************
 *       GMWLObservation.hpp  -  Multi-wavelength observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
 */

#ifndef GMWLOBSERVATION_HPP
#define GMWLOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GMWLResponse.hpp"
#include "GMWLPointing.hpp"


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
    virtual void             clear(void);
    virtual GMWLObservation* clone(void) const;
    virtual void             response(const GResponse& rsp);
    virtual GMWLResponse*    response(void) const;
    virtual GMWLPointing*    pointing(const GTime& time) const;
    virtual std::string      instrument(void) const { return "MWL"; }
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;
    virtual std::string      print(void) const;

    // Other methods
    void        load(const std::string& filename);
    void        load(const std::string& filename, const int& extno);
    void        load(const std::string& filename, const std::string& extname);
    std::string filename(void) const { return m_filename; }
    std::string extno(void) const { return m_extno; }
    std::string extname(void) const { return m_extname; }
    void        filename(const std::string& filename) { m_filename=filename; }
    void        extno(const std::string& extno) { m_extno=extno; }
    void        extname(const std::string& extname) { m_extname=extname; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLObservation& obs);
    void free_members(void);

    // Protected members
    std::string   m_instrument;   //!< Instrument name
    std::string   m_filename;     //!< Filename
    std::string   m_extno;        //!< Extension number
    std::string   m_extname;      //!< Extension name
    GMWLResponse* m_response;     //!< Pointer to response functions
    GMWLPointing* m_pointing;     //!< Pointer to pointing direction
};

#endif /* GMWLOBSERVATION_HPP */
