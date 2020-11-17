/***************************************************************************
 *         GCTAModelSpatialLookup.i - Spatial lookup table model           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAModelSpatialLookup.i
 * @brief Spatial lookup table model interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelSpatialLookup.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelSpatialLookup
 *
 * @brief Spatial lookup table model class
 ***************************************************************************/
class GCTAModelSpatialLookup  : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialLookup(void);
    explicit GCTAModelSpatialLookup(const GFilename& filename);
    explicit GCTAModelSpatialLookup(const GCTAResponseTable& table);
    explicit GCTAModelSpatialLookup(const GXmlElement& xml);
    GCTAModelSpatialLookup(const double&   radmax,
                           const double&   radbin,
                           const GEbounds& ebds);
    GCTAModelSpatialLookup(const GCTAModelSpatialLookup& model);
    virtual ~GCTAModelSpatialLookup(void);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelSpatialLookup* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const GCTAInstDir& dir,
                                         const GEnergy&     energy,
                                         const GTime&       time,
                                         const bool&        gradients = false) const;
    virtual double                  mc_max_value(const GCTAObservation& obs) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;

    // Other methods
    void                     fill(const GCTAObservation& obs);
    void                     fill(const GObservations& obs);
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     load(const GFilename& filename);
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;
    double                   norm(void) const;
    void                     norm(const double& norm);
};


/***********************************************************************//**
 * @brief GCTAModelSpatialLookup class extension
 ***************************************************************************/
%extend GCTAModelSpatialLookup {
    GCTAModelSpatialLookup copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        if state[0].elements('parameter') == 0:
            self.__init__()
        else:
            self.__init__(state[0])
}
};
