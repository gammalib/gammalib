/***************************************************************************
 *     GModelSpectralLogParabola.i - Log parabola spectral model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Michael Mayer                               *
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
 * @file GModelSpectralLogParabola.i
 * @brief Log parabola spectral model class definition
 * @author Michael Mayer 
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralLogParabola.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralLogParabola
 *
 * @brief LogParabola spectral model class
 ***************************************************************************/
class GModelSpectralLogParabola : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralLogParabola(void);
    explicit GModelSpectralLogParabola(const double& norm,
                                       const double& index,
                                       const double& curvature);
    explicit GModelSpectralLogParabola(const GXmlElement& xml);
    GModelSpectralLogParabola(const GModelSpectralLogParabola& model);
    virtual ~GModelSpectralLogParabola(void);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GModelSpectralLogParabola* clone(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GEnergy& srcEng,
                                            const GTime&   srcTime) const;
    virtual double                     eval_gradients(const GEnergy& srcEng,
                                                      const GTime&   srcTime);
    virtual double                     flux(const GEnergy& emin,
                                            const GEnergy& emax) const;
    virtual double                     eflux(const GEnergy& emin,
                                             const GEnergy& emax) const;
    virtual GEnergy                    mc(const GEnergy& emin,
                                          const GEnergy& emax,
                                          const GTime&   time,
                                          GRan&          ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;

    // Other methods
    double norm(void) const;
    double index(void) const;
    double curvature(void) const;
    double pivot(void) const;

};


/***********************************************************************//**
 * @brief GModelSpectralLogParabola class extension
 ***************************************************************************/
%extend GModelSpectralLogParabola {
    GModelSpectralLogParabola copy() {
        return (*self);
    }
};
