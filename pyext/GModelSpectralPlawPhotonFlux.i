/***************************************************************************
 *    GModelSpectralPlawPhotonFlux.i - Spectral power law model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlawPhotonFlux.i
 * @brief Flux normalized power law spectral model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlawPhotonFlux.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlawPhotonFlux
 *
 * @brief Flux normalized power law spectral model class
 ***************************************************************************/
class GModelSpectralPlawPhotonFlux : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlawPhotonFlux(void);
    GModelSpectralPlawPhotonFlux(const std::string& type,
                        const std::string& photon_flux,
                        const std::string& index,
                        const std::string& emin,
                        const std::string& emax);
    GModelSpectralPlawPhotonFlux(const double&  photon_flux,
                        const double&  index,
                        const GEnergy& emin,
                        const GEnergy& emax);
    explicit GModelSpectralPlawPhotonFlux(const GXmlElement& xml);
    GModelSpectralPlawPhotonFlux(const GModelSpectralPlawPhotonFlux& model);
    virtual ~GModelSpectralPlawPhotonFlux(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralPlawPhotonFlux* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime = GTime()) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime&   srcTime = GTime());
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
    double  photon_flux(void) const;
    void    photon_flux(const double& photon_flux);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(void) const;
    void    emin(const GEnergy& emin);
    GEnergy emax(void) const;
    void    emax(const GEnergy& emax);
};


/***********************************************************************//**
 * @brief GModelSpectralPlawPhotonFlux class extension
 ***************************************************************************/
%extend GModelSpectralPlawPhotonFlux {
    GModelSpectralPlawPhotonFlux copy() {
        return (*self);
    }
};
