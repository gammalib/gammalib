/***************************************************************************
 *            GCOMObservation.i  -  COMPTEL observation class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
    explicit GCOMObservation(const std::string& drename,
                             const std::string& drbname,
                             const std::string& drgname,
                             const std::string& drxname);
    GCOMObservation(const GCOMObservation& obs);
    virtual ~GCOMObservation(void);

    // Implement pure virtual methods
    virtual void             clear(void);
    virtual GCOMObservation* clone(void) const;
    virtual void             response(const GResponse& rsp);
    virtual GCOMResponse*    response(void) const;
    virtual GCOMPointing*    pointing(void) const;
    virtual std::string      instrument(void) const;
    virtual double           ontime(void) const;
    virtual double           livetime(void) const;
    virtual double           deadc(const GTime& time) const;
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;

    // Other methods
    void           load(const std::string& drename,
                        const std::string& drbname,
                        const std::string& drgname,
                        const std::string& drxname);
    void           response(const std::string& iaqname,
                            const std::string& caldb = "");
    void           obs_id(const double& id);
    void           ontime(const double& ontime);
    void           livetime(const double& livetime);
    void           deadc(const double& deadc);
    double         obs_id(void) const;
    const GSkymap& drb(void) const;
    const GSkymap& drg(void) const;
    const GSkymap& drx(void) const;
};


/***********************************************************************//**
 * @brief GCOMObservation class extension
 ***************************************************************************/
%extend GCOMObservation {
    GCOMObservation copy() {
        return (*self);
    }
};
