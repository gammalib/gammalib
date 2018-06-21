/***************************************************************************
 *                  GModelSpectralSmoothBrokenPlaw.i                       *
 *               Smoothly broken power law spectrum class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2018 by Josh Cardenzana                             *
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
 * @file GModelSpectralSmoothBrokenPlaw.i
 * @brief Smoothly broken power law spectrum class definition
 * @author Josh Cardenzana
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralSmoothBrokenPlaw.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralSmoothBrokenPlaw
 *
 * @brief Smoothly broken Power law spectral model class
 ***************************************************************************/
class GModelSpectralSmoothBrokenPlaw : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralSmoothBrokenPlaw(void);
    GModelSpectralSmoothBrokenPlaw(const std::string& type,
                                   const std::string& prefactor,
                                   const std::string& index1,
                                   const std::string& pivot,
                                   const std::string& index2,
                                   const std::string& breakenergy,
                                   const std::string& beta);
    GModelSpectralSmoothBrokenPlaw(const double&  prefactor,
                                   const double&  index1,
                                   const GEnergy& pivot,
                                   const double&  index2,
                                   const GEnergy& breakenergy,
                                   const double&  beta);
    explicit GModelSpectralSmoothBrokenPlaw(const GXmlElement& xml);
    GModelSpectralSmoothBrokenPlaw(const GModelSpectralSmoothBrokenPlaw& model);
    virtual ~GModelSpectralSmoothBrokenPlaw(void);
    
    // Implemented pure virtual methods
    virtual void                            clear(void);
    virtual GModelSpectralSmoothBrokenPlaw* clone(void) const;
    virtual std::string                     classname(void) const;
    virtual std::string                     type(void) const;
    virtual double                          eval(const GEnergy& srcEng,
                                                 const GTime& srcTime = GTime(),
                                                 const bool& gradients = false) const;
    virtual double                          flux(const GEnergy& emin,
                                                 const GEnergy& emax) const;
    virtual double                          eflux(const GEnergy& emin,
                                                  const GEnergy& emax) const;
    virtual GEnergy                         mc(const GEnergy& emin,
                                               const GEnergy& emax,
                                               const GTime&   time,
                                               GRan&          ran) const;
    virtual void                            read(const GXmlElement& xml);
    virtual void                            write(GXmlElement& xml) const;
    
    // Other methods
    void    type(const std::string& type);
    double  prefactor(void) const;
    double  index1(void) const;
    double  index2(void) const;
    GEnergy pivot(void) const;
    GEnergy breakenergy(void) const;
    double  beta(void) const;
    void    prefactor(const double& prefactor);
    void    index1(const double& index1);
    void    index2(const double& index2);
    void    pivot(const GEnergy& pivot);
    void    breakenergy(const GEnergy& breakenergy);
    void    beta(const double& beta);
};


/***********************************************************************//**
 * @brief GModelSpectralSmoothBrokenPlaw class extension
 ***************************************************************************/
%extend GModelSpectralSmoothBrokenPlaw {
    GModelSpectralSmoothBrokenPlaw copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.type(), self[0], self[1], self[2], self[3], self[4], self[5])
        return state
    def __setstate__(self, state):
        self.__init__()
        self.type(state[0])
        self[0] = state[1]
        self[1] = state[2]
        self[2] = state[3]
        self[3] = state[4]
        self[4] = state[5]
        self[5] = state[6]
}
};
