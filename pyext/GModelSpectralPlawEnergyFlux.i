/***************************************************************************
 *    GModelSpectralPlawEnergyFlux.i - Spectral power law model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                                    *
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
 * @file GModelSpectralPlawEnergyFlux.i
 * @brief Energy flux normalized power law spectral model class Python interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlawEnergyFlux.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlawEnergyFlux
 *
 * @brief Flux normalized power law spectral model class
 ***************************************************************************/
class GModelSpectralPlawEnergyFlux : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlawEnergyFlux(void);
    GModelSpectralPlawEnergyFlux(const std::string& type,
                                 const std::string& eflux,
                                 const std::string& index,
                                 const std::string& emin,
                                 const std::string& emax);
    explicit GModelSpectralPlawEnergyFlux(const double&  eflux,
                                          const double&  index,
                                          const GEnergy& emin,
                                          const GEnergy& emax);
    explicit GModelSpectralPlawEnergyFlux(const GXmlElement& xml);
    GModelSpectralPlawEnergyFlux(const GModelSpectralPlawEnergyFlux& model);
    virtual ~GModelSpectralPlawEnergyFlux(void);

    // Implemented pure virtual base class methods
    virtual void                          clear(void);
    virtual GModelSpectralPlawEnergyFlux* clone(void) const;
    virtual std::string                   classname(void) const;
    virtual std::string                   type(void) const;
    virtual double                        eval(const GEnergy& srcEng,
                                               const GTime&   srcTime = GTime()) const;
    virtual double                        eval_gradients(const GEnergy& srcEng,
                                                         const GTime&   srcTime = GTime());
    virtual double                        flux(const GEnergy& emin,
                                               const GEnergy& emax) const;
    virtual double                        eflux(const GEnergy& emin,
                                                const GEnergy& emax) const;
    virtual GEnergy                       mc(const GEnergy& emin,
                                             const GEnergy& emax,
                                             const GTime&   time,
                                             GRan&          ran) const;
    virtual void                          read(const GXmlElement& xml);
    virtual void                          write(GXmlElement& xml) const;

    // Other methods
    double  eflux(void) const;
    void    eflux(const double& eflux);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(void) const;
    void    emin(const GEnergy& emin);
    GEnergy emax(void) const;
    void    emax(const GEnergy& emax);
};


/***********************************************************************//**
 * @brief GModelSpectralPlawEnergyFlux class extension
 ***************************************************************************/
%extend GModelSpectralPlawEnergyFlux {
    GModelSpectralPlawEnergyFlux copy() {
        return (*self);
    }
};
