/***************************************************************************
 *              GLATObservation.i - Fermi/LAT observation class            *
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
 * @file GLATObservation.i
 * @brief Fermi/LAT observation class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATObservation.hpp"
%}


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Fermi/LAT observation class
 ***************************************************************************/
class GLATObservation : public GObservation {
public:
    // Constructors and destructors
    GLATObservation();
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation();

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GLATObservation*    clone(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GLATResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

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
};


/***********************************************************************//**
 * @brief GLATObservation class extension
 ***************************************************************************/
%extend GLATObservation {
    GLATObservation copy() {
        return (*self);
    }
};
