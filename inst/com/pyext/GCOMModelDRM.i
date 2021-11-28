/***************************************************************************
 *                GCOMModelDRM.i - COMPTEL DRM model class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GCOMModelDRM.hpp
 * @brief COMPTEL DRM model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMModelDRM.hpp"
%}


/***********************************************************************//**
 * @class GCOMModelDRM
 *
 * @brief COMPTEL DRM model class
 ***************************************************************************/
class GCOMModelDRM : public GModelData {
public:
    // Constructors and destructors
    GCOMModelDRM(void);
    explicit GCOMModelDRM(const GXmlElement& xml);
    GCOMModelDRM(const GCOMModelDRM& model);
    virtual ~GCOMModelDRM(void);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GCOMModelDRM*  clone(void) const;
    virtual std::string    classname(void) const;
    virtual std::string    type(void) const;
    virtual bool           is_constant(void) const;
    virtual double         eval(const GEvent& event,
                                const GObservation& obs,
                                const bool& gradients = false) const;
    virtual double         npred(const GEnergy& obsEng, const GTime& obsTime,
                                 const GObservation& obs) const;
    virtual GCOMEventCube* mc(const GObservation& obs, GRan& ran) const;
    virtual void           read(const GXmlElement& xml);
    virtual void           write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCOMModelDRM class extension
 ***************************************************************************/
%extend GCOMModelDRM {
    GCOMModelDRM copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = {'xml': xml}
        return state
    def __setstate__(self, state):
        if state['xml'].elements('parameter') == 0:
            self.__init__()
        else:
            self.__init__(state['xml'])
}
};
