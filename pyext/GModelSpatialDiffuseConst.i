/***************************************************************************
 *        GModelSpatialDiffuseConst.i - Spatial isotropic model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseConst.i
 * @brief Isotropic spatial model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialDiffuseConst.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialDiffuseConst
 *
 * @brief Isotropic spatial model
 ***************************************************************************/
class GModelSpatialDiffuseConst  : public GModelSpatialDiffuse {
public:
    // Constructors and destructors
    GModelSpatialDiffuseConst(void);
    explicit GModelSpatialDiffuseConst(const GXmlElement& xml);
    explicit GModelSpatialDiffuseConst(const double& value);
    GModelSpatialDiffuseConst(const GModelSpatialDiffuseConst& model);
    virtual ~GModelSpatialDiffuseConst(void);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GModelSpatialDiffuseConst* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GPhoton& photon,
                                            const bool& gradients = false) const;
    virtual GSkyDir                    mc(const GEnergy& energy,
                                          const GTime& time,
                                          GRan& ran) const;
    virtual double                     mc_norm(const GSkyDir& dir,
                                               const double&  radius) const;
    virtual bool                       contains(const GSkyDir& dir,
                                                const double&  margin = 0.0) const;
    virtual GSkyRegion*                region(void) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;

    // Other methods
    double value(void) const;
    void   value(const double& value);
};


/***********************************************************************//**
 * @brief GModelSpatialDiffuseConst class extension
 ***************************************************************************/
%extend GModelSpatialDiffuseConst {
    GModelSpatialDiffuseConst copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = xml,
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
