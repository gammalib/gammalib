/***************************************************************************
 *          GModelSpectralGauss.i - Spectral gaussian model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Christoph Deil & Ellis Owen                      *
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
 * @file GModelSpectralGauss.i
 * @brief Spectral constant model class interface definition
 * @author Ellis Owen, Christoph Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralGauss.hpp"
#include "GTools.hpp"
#include "GEnergy.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralGauss
 *
 * @brief Spectral gaussian model class
 ***************************************************************************/
class GModelSpectralGauss : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralGauss(void);
    explicit GModelSpectralGauss(const GXmlElement& xml);
    GModelSpectralGauss(const double&  norm,
                        const GEnergy& mean,
                        const GEnergy& sigma);
    GModelSpectralGauss(const GModelSpectralGauss& model);
    virtual ~GModelSpectralGauss(void);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralGauss* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime& srcTime) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime& srcTime);
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
    double  norm(void) const;
    void    norm(const double& norm);
    GEnergy mean(void) const;
    void    mean(const GEnergy& mean);
    GEnergy sigma(void) const;
    void    sigma(const GEnergy& sigma);
};


/***********************************************************************//**
 * @brief GModelSpectralGauss class extension
 ***************************************************************************/
%extend GModelSpectralGauss {
    GModelSpectralGauss copy() {
        return (*self);
    }
};
