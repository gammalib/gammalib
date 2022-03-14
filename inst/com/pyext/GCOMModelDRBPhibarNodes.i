/***************************************************************************
 *        GCOMModelDRBPhibarNodes.i - COMPTEL DRB model fitting class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2022 by Juergen Knoedlseder                         *
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
 * @file GCOMModelDRBPhibarNodes.hpp
 * @brief COMPTEL DRB Phibar nodes model fitting class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMModelDRBPhibarNodes.hpp"
%}


/***********************************************************************//**
 * @class GCOMModelDRBPhibarNodes
 *
 * @brief COMPTEL DRB model fitting class
 ***************************************************************************/
class GCOMModelDRBPhibarNodes : public GModelData {
public:
    // Constructors and destructors
    GCOMModelDRBPhibarNodes(void);
    explicit GCOMModelDRBPhibarNodes(const GXmlElement& xml);
    GCOMModelDRBPhibarNodes(const GCOMModelDRBPhibarNodes& model);
    virtual ~GCOMModelDRBPhibarNodes(void);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCOMModelDRBPhibarNodes* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent& event,
                                          const GObservation& obs,
                                          const bool& gradients = false) const;
    virtual double                   npred(const GEnergy& obsEng, const GTime& obsTime,
                                           const GObservation& obs) const;
    virtual GCOMEventCube*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCOMModelDRBPhibarNodes class extension
 ***************************************************************************/
%extend GCOMModelDRBPhibarNodes {
    GCOMModelDRBPhibarNodes copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        if state[0][0].elements('node') == 0:
            self.__init__()
        else:
            self.__init__(state[0][0])
}
};
