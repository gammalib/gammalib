/***************************************************************************
 *         GCTAObservation.i  -  CTA Observation class interface           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GCTAObservation.i
 * @brief CTA observation class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAObservation.hpp"
%}


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief CTA observation class Python interface
 ***************************************************************************/
class GCTAObservation : public GObservation {
public:
    // Constructors and destructors
    GCTAObservation(void);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Implemented pure virtual base class methods
    virtual void             clear(void);
    virtual GCTAObservation* clone(void) const;
    virtual void             response(const GResponse& rsp);
    virtual GCTAResponse*    response(void) const;
    virtual GCTAPointing*    pointing(void) const;
    virtual std::string      instrument(void) const;
    virtual double           ontime(void) const;
    virtual double           livetime(void) const;
    virtual double           deadc(const GTime& time) const;
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;

    // Other methods
    void        load_unbinned(const std::string& filename);
    void        load_binned(const std::string& filename);
    void        save(const std::string& filename, bool clobber) const;
    void        response(const std::string& irfname, std::string caldb = "");
    void        pointing(const GCTAPointing& pointing);
    void        obs_id(const int& id);
    void        ra_obj(const double& ra);
    void        dec_obj(const double& dec);
    void        ontime(const double& ontime);
    void        livetime(const double& livetime);
    void        deadc(const double& deadc);
    int         obs_id(void) const;
    double      ra_obj(void) const;
    double      dec_obj(void) const;
    std::string eventfile(void) const;
    void        eventfile(const std::string& filename);
};


/***********************************************************************//**
 * @brief GCTAObservation class extension
 ***************************************************************************/
%extend GCTAObservation {
    GCTAObservation copy() {
        return (*self);
    }
};
