/***************************************************************************
 *          GModelSpectralPlaw.i - Spectral power law model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlaw.i
 * @brief Power law spectral model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlaw.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlaw
 *
 * @brief Power law spectral model class
 ***************************************************************************/
class GModelSpectralPlaw : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlaw(void);
    explicit GModelSpectralPlaw(const double&  prefactor,
                                const double&  index,
                                const GEnergy& pivot);
    explicit GModelSpectralPlaw(const GXmlElement& xml);
    GModelSpectralPlaw(const GModelSpectralPlaw& model);
    virtual ~GModelSpectralPlaw(void);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralPlaw* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime&   srcTime = GTime()) const;
    virtual double              eval_gradients(const GEnergy& srcEng,
                                               const GTime&   srcTime = GTime());
    virtual double              flux(const GEnergy& emin,
                                     const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin,
                                   const GEnergy& emax,
                                   const GTime&   time,
                                   GRan&          ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    double  prefactor(void) const;
    double  index(void) const;
    GEnergy pivot(void) const;
    void    prefactor(const double& prefactor);
    void    index(const double& index);
    void    pivot(const GEnergy& pivot);
};


/***********************************************************************//**
 * @brief GModelSpectralPlaw class extension
 ***************************************************************************/
%extend GModelSpectralPlaw {
    GModelSpectralPlaw copy() {
        return (*self);
    }
};
