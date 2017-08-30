/***************************************************************************
 *             GCOMObservation.i - COMPTEL observation class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file GCOMObservation.i
 * @brief COMPTEL observation class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMObservation.hpp"
%}


/***********************************************************************//**
 * @class GCOMObservation
 *
 * @brief Interface class for COMPTEL observations
 ***************************************************************************/
class GCOMObservation : public GObservation {
public:
    // Constructors and destructors
    GCOMObservation(void);
    GCOMObservation(const GFilename& drename,
                    const GFilename& drbname,
                    const GFilename& drgname,
                    const GFilename& drxname);
    GCOMObservation(const GCOMObservation& obs);
    virtual ~GCOMObservation(void);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GCOMObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GCOMResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    void           load(const GFilename& drename,
                        const GFilename& drbname,
                        const GFilename& drgname,
                        const GFilename& drxname);
    void           response(const GCaldb& caldb, const std::string& rspname);
    void           obs_id(const double& id);
    void           ontime(const double& ontime);
    void           livetime(const double& livetime);
    void           deadc(const double& deadc);
    void           ewidth(const double& ewidth);
    const double&  obs_id(void) const;
    const double&  ewidth(void) const;
    const GSkyMap& drb(void) const;
    const GSkyMap& drg(void) const;
    const GSkyMap& drx(void) const;
};


/***********************************************************************//**
 * @brief GCOMObservation class extension
 ***************************************************************************/
%extend GCOMObservation {
    GCOMObservation copy() {
        return (*self);
    }
};
