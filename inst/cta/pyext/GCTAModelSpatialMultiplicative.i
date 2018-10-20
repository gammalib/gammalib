/***************************************************************************
 *  GCTAModelSpatialMultiplicative.i - Multiplicative spatial model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialMultiplicative.i
 * @brief Multiplicative spatial model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelSpatialMultiplicative.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelSpatialMultiplicative
 *
 * @brief Multiplicative spatial model class
 ***************************************************************************/
class GCTAModelSpatialMultiplicative : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialMultiplicative(void);
    explicit GCTAModelSpatialMultiplicative(const GXmlElement& xml);
    GCTAModelSpatialMultiplicative(const GCTAModelSpatialMultiplicative& model);
    virtual ~GCTAModelSpatialMultiplicative(void);

    // Operators
    virtual GCTAModelSpatialMultiplicative& operator=(const GCTAModelSpatialMultiplicative& model);

    // Pure virtual methods
    virtual void                            clear(void);
    virtual GCTAModelSpatialMultiplicative* clone(void) const;
    virtual std::string                     classname(void) const;
    virtual std::string                     type(void) const;
    virtual double                          eval(const GCTAInstDir& dir,
                                                 const GEnergy&     energy,
                                                 const GTime&       time,
                                                 const bool&        gradients = false) const;
    virtual GCTAInstDir                     mc(const GEnergy& energy,
                                               const GTime&   time,
                                               GRan& ran) const;
    virtual void                            read(const GXmlElement& xml);
    virtual void                            write(GXmlElement& xml) const;

    // Other methods
    void                    append(const GCTAModelSpatial& spatial,
                                   const std::string&      name="");
    int                     components(void) const;
    const GCTAModelSpatial* component(const int& index) const;
    const GCTAModelSpatial* component(const std::string& name) const;
};


/***********************************************************************//**
 * @brief GCTAModelSpatialMultiplicative class extension
 ***************************************************************************/
%extend GCTAModelSpatialMultiplicative {
    GCTAModelSpatialMultiplicative copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
