/***************************************************************************
 *           GModelSpectralNodes.i - Spectral nodes model class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralNodes.i
 * @brief Spectral nodes model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralNodes.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralNodes
 *
 * @brief Spectral nodes model class
 ***************************************************************************/
class GModelSpectralNodes : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralNodes(void);
    GModelSpectralNodes(const GModelSpectral& model, const GEnergies& energies);
    explicit GModelSpectralNodes(const GXmlElement& xml);
    GModelSpectralNodes(const GModelSpectralNodes& model);
    virtual ~GModelSpectralNodes(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralNodes* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime& srcTime = GTime(),
                                      const bool& gradients = false) const;
    virtual double               flux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin,
                                       const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin,
                                    const GEnergy& emax,
                                    const GTime&   time,
                                    GRan&          ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;

    // Other methods
    int     nodes(void) const;
    void    append(const GEnergy& energy, const double& intensity);
    void    insert(const int& index, const GEnergy& energy,
                   const double& intensity);
    void    remove(const int& index);
    void    reserve(const int& num);
    void    extend(const GModelSpectralNodes& nodes);
    GEnergy energy(const int& index) const;
    void    energy(const int& index, const GEnergy& energy);
    double  intensity(const int& index) const;
    void    intensity(const int& index, const double& intensity);
};


/***********************************************************************//**
 * @brief GModelSpectralNodes class extension
 ***************************************************************************/
%extend GModelSpectralNodes {
    GModelSpectralNodes copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = xml,
        return state
    def __setstate__(self, state):
        if state[0].elements('node') == 0:
            self.__init__()
        else:
            self.__init__(state[0])
}
};
