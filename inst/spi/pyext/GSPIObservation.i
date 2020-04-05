/***************************************************************************
 *            GSPIObservation.i - INTEGRAL/SPI observation class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIObservation.i
 * @brief INTEGRAL/SPI observation class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSPIObservation.hpp"
%}


/***********************************************************************//**
 * @class GSPIObservation
 *
 * @brief INTEGRAL/SPI observation class
 *
 * The INTEGRAL/SPI observation class defines an observation.
 ***************************************************************************/
class GSPIObservation : public GObservation {

public:
    // Constructors and destructors
    GSPIObservation(void);
    explicit GSPIObservation(const GXmlElement& xml);
    explicit GSPIObservation(const GFilename& filename);
    GSPIObservation(const GSPIObservation& obs);
    virtual ~GSPIObservation(void);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GSPIObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GSPIResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    void read(const GFits& fits);
    void load(const GFilename& filename);
    void ontime(const double& ontime);
    void livetime(const double& livetime);
    void deadc(const double& deadc);
};


/***********************************************************************//**
 * @brief GSPIObservation class extension
 ***************************************************************************/
%extend GSPIObservation {
    GSPIObservation copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.classname()) # TODO: Replace by appropriate class members
        return state
    def __setstate__(self, state):
        self.__init__()
}
};
