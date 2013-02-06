/***************************************************************************
 *              GXXXObservation.i  -  XXX observation class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXObservation.i
 * @brief XXX observation class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXObservation.hpp"
%}


/***********************************************************************//**
 * @class GXXXObservation
 *
 * @brief Interface class for XXX observations
 ***************************************************************************/
class GXXXObservation : public GObservation {
public:
    // Constructors and destructors
    GXXXObservation(void);
    GXXXObservation(const GXXXObservation& obs);
    virtual ~GXXXObservation(void);

    // Implement pure virtual methods
    virtual void             clear(void);
    virtual GXXXObservation* clone(void) const;
    virtual void             response(const GResponse& rsp);
    virtual GXXXResponse*    response(void) const;
    virtual GXXXPointing*    pointing(void) const;
    virtual std::string      instrument(void) const;
    virtual double           ontime(void) const;
    virtual double           livetime(void) const;
    virtual double           deadc(const GTime& time) const;
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GXXXObservation class extension
 ***************************************************************************/
%extend GXXXObservation {
    GXXXObservation copy() {
        return (*self);
    }
};
