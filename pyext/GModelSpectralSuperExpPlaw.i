/***************************************************************************
 *      GModelSpectralSuperExpPlaw.i - Super Exponential power law model   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralSuperExpPlaw.i
 * @brief Super Exponential cut off power law spectral class interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralSuperExpPlaw.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralSuperExpPlaw
 *
 * @brief Super Exponential cut off power law spectral class
 ***************************************************************************/
class GModelSpectralSuperExpPlaw : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralSuperExpPlaw(void);
    explicit GModelSpectralSuperExpPlaw(const double&  prefactor,
                                   const double&  index1,
                                   const GEnergy& pivot,
                                   const GEnergy& cutoff,
                                   const double index2);
    explicit GModelSpectralSuperExpPlaw(const GXmlElement& xml);
    GModelSpectralSuperExpPlaw(const GModelSpectralSuperExpPlaw& model);
    virtual ~GModelSpectralSuperExpPlaw(void);

    // Implemented pure virtual methods
    virtual void                   clear(void);
    virtual GModelSpectralSuperExpPlaw* clone(void) const;
    virtual std::string            type(void) const;
    virtual double                 eval(const GEnergy& srcEng,
                                        const GTime&   srcTime) const;
    virtual double                 eval_gradients(const GEnergy& srcEng,
                                                  const GTime&   srcTime);
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
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index1(void) const;
    void    index1(const double& index1);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);
    GEnergy cutoff(void) const;
    void    cutoff(const GEnergy& cutoff);
    double  index2(void) const;
    void    index2(const double& index2);
};


/***********************************************************************//**
 * @brief GModelSpectralExpPlaw class extension
 ***************************************************************************/
%extend GModelSpectralSuperExpPlaw {
    GModelSpectralSuperExpPlaw copy() {
        return (*self);
    }
};
