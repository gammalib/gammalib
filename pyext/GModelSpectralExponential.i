/***************************************************************************
 *   GModelSpectralExponential.i - Exponential spectral model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2020 Luigi Tibaldo                                  *
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
 * @file GModelSpectralExponential.i
 * @brief Exponential spectral model class interface definition
 * @author Luigi Tibaldo
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralExponential.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralExponential
 *
 * @brief Exponential spectral model class
 ***************************************************************************/
class GModelSpectralExponential : public GModelSpectral {
  
public:
  
    // Constructors and destructors
    GModelSpectralExponential(void);
    explicit GModelSpectralExponential(const GXmlElement& xml);
    explicit GModelSpectralExponential(const GModelSpectral* spec);
    GModelSpectralExponential(const GModelSpectralExponential& model);
    virtual ~GModelSpectralExponential(void);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GModelSpectralExponential* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GEnergy& srcEng,
                                            const GTime&   srcTime = GTime(),
                                            const bool&    gradients = false) const;
    virtual double                     flux(const GEnergy& emin,
                                            const GEnergy& emax) const;
    virtual double                     eflux(const GEnergy& emin,
                                             const GEnergy& emax) const;
    virtual GEnergy                    mc(const GEnergy& emin,
                                          const GEnergy& emax,
                                          const GTime&   time,
                                          GRan&          ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;

    // Other methods
    void                  exponent(const GModelSpectral* spec);
    const GModelSpectral* exponent(void) const;
    
};


/***********************************************************************//**
 * @brief GModelSpectralExponential class extension
 ***************************************************************************/
%extend GModelSpectralExponential {
    GModelSpectralExponential copy() {
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
