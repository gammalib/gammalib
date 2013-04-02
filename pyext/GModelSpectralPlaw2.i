/***************************************************************************
 *         GModelSpectralPlaw2.i - Spectral power law model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlaw2.i
 * @brief Flux normalized power law spectral model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlaw2.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlaw2
 *
 * @brief Flux normalized power law spectral model class
 ***************************************************************************/
class GModelSpectralPlaw2 : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlaw2(void);
    explicit GModelSpectralPlaw2(const double&  integral,
                                 const double&  index,
                                 const GEnergy& emin,
                                 const GEnergy& emax);
    explicit GModelSpectralPlaw2(const GXmlElement& xml);
    GModelSpectralPlaw2(const GModelSpectralPlaw2& model);
    virtual ~GModelSpectralPlaw2(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralPlaw2* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime&   srcTime);
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
    double  integral(void) const;
    void    integral(const double& integral);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(void) const;
    void    emin(const GEnergy& emin);
    GEnergy emax(void) const;
    void    emax(const GEnergy& emax);
};


/***********************************************************************//**
 * @brief GModelSpectralPlaw2 class extension
 ***************************************************************************/
%extend GModelSpectralPlaw2 {
    GModelSpectralPlaw2 copy() {
        return (*self);
    }
};
