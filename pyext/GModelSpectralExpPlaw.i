/***************************************************************************
 *      GModelSpectralExpPlaw.i - Exponential cut off power law model      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralExpPlaw.i
 * @brief Exponential cut off power law spectral class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralExpPlaw.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralExpPlaw
 *
 * @brief Exponential cut off power law spectral class
 ***************************************************************************/
class GModelSpectralExpPlaw : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralExpPlaw(void);
    GModelSpectralExpPlaw(const std::string& type,
                          const std::string& prefactor,
                          const std::string& index,
                          const std::string& pivot,
                          const std::string& cutoff);
    GModelSpectralExpPlaw(const double&  prefactor,
                          const double&  index,
                          const GEnergy& pivot,
                          const GEnergy& cutoff);
    explicit GModelSpectralExpPlaw(const GXmlElement& xml);
    GModelSpectralExpPlaw(const GModelSpectralExpPlaw& model);
    virtual ~GModelSpectralExpPlaw(void);

    // Implemented pure virtual methods
    virtual void                   clear(void);
    virtual GModelSpectralExpPlaw* clone(void) const;
    virtual std::string            classname(void) const;
    virtual std::string            type(void) const;
    virtual double                 eval(const GEnergy& srcEng,
                                        const GTime&   srcTime = GTime(),
                                        const bool&    gradients = false) const;
    virtual double                 flux(const GEnergy& emin,
                                        const GEnergy& emax) const;
    virtual double                 eflux(const GEnergy& emin,
                                         const GEnergy& emax) const;
    virtual GEnergy                mc(const GEnergy& emin,
                                      const GEnergy& emax,
                                      const GTime&   time,
                                      GRan&          ran) const;
    virtual void                   read(const GXmlElement& xml);
    virtual void                   write(GXmlElement& xml) const;

    // Other methods
    void    type(const std::string& type);
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index(void) const;
    void    index(const double& index);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);
    GEnergy cutoff(void) const;
    void    cutoff(const GEnergy& cutoff);
};


/***********************************************************************//**
 * @brief GModelSpectralExpPlaw class extension
 ***************************************************************************/
%extend GModelSpectralExpPlaw {
    GModelSpectralExpPlaw copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.type(), self[0], self[1], self[2], self[3])
        return state
    def __setstate__(self, state):
        self.__init__()
        self.type(state[0])
        self[0] = state[1]
        self[1] = state[2]
        self[2] = state[3]
        self[3] = state[4]
}
};
