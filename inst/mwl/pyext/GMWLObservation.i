/***************************************************************************
 *         GMWLObservation.i - Multi-wavelength observation class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GMWLObservation.i
 * @brief Multi-wavelength observation class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLObservation.hpp"
%}


/***********************************************************************//**
 * @class GMWLObservation
 *
 * @brief Python interface class for multi-wavelength observations
 *
 * This class implements a multi-wavelength observation. A multi-wavelength
 * observation contains spectral points obtained with an unspecified
 * instrument. The spectral points are given in physical units.
 ***************************************************************************/
class GMWLObservation : public GObservation {

public:
    // Constructors and destructors
    GMWLObservation(void);
    explicit GMWLObservation(const GFilename& filename);
    GMWLObservation(const GMWLObservation& obs);
    virtual ~GMWLObservation(void);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GMWLObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GMWLResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    void             load(const GFilename& filename);
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
};


/***********************************************************************//**
 * @brief GMWLObservation class extension
 ***************************************************************************/
%extend GMWLObservation {
    GMWLObservation copy() {
        return (*self);
    }
};
