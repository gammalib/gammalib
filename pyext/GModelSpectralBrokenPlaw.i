/***************************************************************************
 *      GModelSpectralBrokenPlaw.i - Broken power law spectrum class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Anneli Schulz                                    *
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
 * @file GModelSpectralBrokenPlaw.i
 * @brief Broken power law spectrum class definition
 * @author Anneli Schulz
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralBrokenPlaw.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralBrokenPlaw
 *
 * @brief Broken Power law spectral model class
 ***************************************************************************/
class GModelSpectralBrokenPlaw : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralBrokenPlaw(void);
    GModelSpectralBrokenPlaw(const double&  prefactor,
                             const double&  index1,
                             const GEnergy& breakenergy,
                             const double&  index2);
    explicit GModelSpectralBrokenPlaw(const GXmlElement& xml);
    GModelSpectralBrokenPlaw(const GModelSpectralBrokenPlaw& model);
    virtual ~GModelSpectralBrokenPlaw(void);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpectralBrokenPlaw* clone(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GEnergy& srcEng,
                                           const GTime&   srcTime) const;
    virtual double                    eval_gradients(const GEnergy& srcEng,
                                                     const GTime&   srcTime);
    virtual double                    flux(const GEnergy& emin,
                                           const GEnergy& emax) const;
    virtual double                    eflux(const GEnergy& emin,
                                            const GEnergy& emax) const;
    virtual GEnergy                   mc(const GEnergy& emin,
                                         const GEnergy& emax,
                                         const GTime&   time,
                                         GRan&          ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    double  prefactor(void) const;
    double  index1(void) const;
    double  index2(void) const;
    GEnergy breakenergy(void) const;
    void    prefactor(const double& prefactor);
    void    index1(const double& index1);
    void    index2(const double& index2);
    void    breakenergy(const GEnergy& breakenergy);
};


/***********************************************************************//**
 * @brief GModelSpectralBrokenPlaw class extension
 ***************************************************************************/
%extend GModelSpectralBrokenPlaw {
    GModelSpectralBrokenPlaw copy() {
        return (*self);
    }
};
