/***************************************************************************
 *      GModelSpectralExpPlaw.i - Exponential cut off power law model      *
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
 * @file GModelSpectralExpPlaw2.i
 * @brief Exponential cut off power law spectral class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralExpPlaw2.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralExpPlaw2
 *
 * @brief Exponential cut off power law spectral class
 ***************************************************************************/
class GModelSpectralExpPlaw2 : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralExpPlaw2(void);
    explicit GModelSpectralExpPlaw2(const double&  prefactor,
                                    const double&  index,
                                    const GEnergy& pivot,
                                    const double& lambda);
    explicit GModelSpectralExpPlaw2(const GXmlElement& xml);
    GModelSpectralExpPlaw2(const GModelSpectralExpPlaw2& model);
    virtual ~GModelSpectralExpPlaw2(void);

    // Operators
    virtual GModelSpectralExpPlaw2& operator=(const GModelSpectralExpPlaw2& model);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GModelSpectralExpPlaw2* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const GEnergy& srcEng,
                                         const GTime&   srcTime = GTime()) const;
    virtual double                  eval_gradients(const GEnergy& srcEng,
                                                  const GTime&   srcTime = GTime());
    virtual double                  flux(const GEnergy& emin,
                                         const GEnergy& emax) const;
    virtual double                  eflux(const GEnergy& emin,
                                          const GEnergy& emax) const;
    virtual GEnergy                 mc(const GEnergy& emin,
                                       const GEnergy& emax,
                                       const GTime&   time,
                                       GRan&          ran) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;
    virtual std::string             print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index(void) const;
    void    index(const double& index);
    double  lambda(void) const;
    void    lambda(const double& lambda);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);
};


/***********************************************************************//**
 * @brief GModelSpectralExpPlaw2 class extension
 ***************************************************************************/
%extend GModelSpectralExpPlaw2 {
    GModelSpectralExpPlaw2 copy() {
        return (*self);
    }
};
