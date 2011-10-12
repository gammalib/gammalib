/***************************************************************************
 *         GCTAObservation.i  -  CTA Observation class interface           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @author J. Knodlseder
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
    virtual GCTAPointing*    pointing(const GTime& time) const;
    virtual std::string      instrument(void) const;

    // Other methods
    void   load_unbinned(const std::string& filename);
    void   load_binned(const std::string& filename);
    void   save(const std::string& filename, bool clobber) const;
    void   response(const std::string& irfname, std::string caldb = "");
    void   pointing(const GCTAPointing& pointing);
    void   obs_id(const int& id);
    void   ra_obj(const double& ra);
    void   dec_obj(const double& dec);
    int    obs_id(void) const;
    double livetime(void) const;
    double ra_obj(void) const;
    double dec_obj(void) const;
};


/***********************************************************************//**
 * @brief GCTAObservation class extension
 ***************************************************************************/
%extend GCTAObservation {
    GCTAObservation copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GCTAObservation type casts
 ***************************************************************************/
%inline %{
    GCTAObservation* cast_GCTAObservation(GObservation* obs) {
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs);
        if (cta == NULL) {
            throw GException::bad_type("cast_GCTAObservation(GObservation* obs)",
                                       "GObservation not of type GCTAObservation");
        }
        return cta;
    }
%}
