/***************************************************************************
 *                   GCTAPointing.i  -  CTA pointing class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Jurgen Knodlseder                           *
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
 * @file GCTAPointing.i
 * @brief CTA pointing class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GHorizDir.hpp"
#include "GCTAPointing.hpp"
%}


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief CTA pointing class
 ***************************************************************************/
class GCTAPointing : public GBase {
public:
    // Constructors and destructors
    GCTAPointing(void);
    explicit GCTAPointing(const GSkyDir& dir);
    explicit GCTAPointing(const GXmlElement& xml);
    explicit GCTAPointing(const GFilename& filename);
    GCTAPointing(const GCTAPointing& pnt);
    virtual ~GCTAPointing(void);

    // Methods
    void           clear(void);
    GCTAPointing*  clone(void) const;
    std::string    classname(void) const;
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);
    GCTAInstDir    instdir(const GSkyDir& skydir) const;
    GSkyDir        skydir(const GCTAInstDir& instdir) const;
    const GMatrix& rot(void) const;
    const double&  zenith(void) const;
    const double&  azimuth(void) const;
    void           zenith(const double& zenith);  
    void           azimuth(const double& azimuth); 
    GHorizDir      dir_horiz(const GTime& time) const;
    void           load(const GFilename& filename);
    void           read(const GFitsTable& table);
    void           read(const GXmlElement& xml);
    void           write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCTAPointing class extension
 ***************************************************************************/
%extend GCTAPointing {
    GCTAPointing copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.dir(), self.zenith(), self.azimuth())
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
        self.zenith(state[1])
        self.azimuth(state[2])
}
};
