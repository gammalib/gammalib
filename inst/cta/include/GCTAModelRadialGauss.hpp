/***************************************************************************
 *      GCTAModelRadialGauss.hpp  -  Radial Gaussian CTA model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GCTAModelRadialGauss.hpp
 * @brief Radial Gaussian model class interface definition
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRADIALGAUSS_HPP
#define GCTAMODELRADIALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GXmlElement.hpp"
#include "GCTAModelRadial.hpp"
#include "GCTAInstDir.hpp"
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialGauss
 *
 * @brief Radial Gaussian CTA model class
 *
 * This class implements the radial function
 * \f[f(\theta) = \exp \left(-\frac{1}{2}
 *                     \left( \frac{\theta^2}{\sigma} \right)^2 \right)\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$\sigma\f$ is the width parameter (in degrees\f$^2\f$).
 *
 * This function represents a Gaussian in \f$\theta^2\f$.
 ***************************************************************************/
class GCTAModelRadialGauss  : public GCTAModelRadial {

public:
    // Constructors and destructors
    GCTAModelRadialGauss(void);
    explicit GCTAModelRadialGauss(const double& sigma);
    explicit GCTAModelRadialGauss(const GXmlElement& xml);
    GCTAModelRadialGauss(const GCTAModelRadialGauss& model);
    virtual ~GCTAModelRadialGauss(void);

    // Operators
    virtual GCTAModelRadialGauss& operator= (const GCTAModelRadialGauss& model);

    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GCTAModelRadialGauss* clone(void) const;
    virtual std::string           type(void) const { return "Gaussian"; }
    virtual double                eval(const double& offset) const;
    virtual double                eval_gradients(const double& offset) const;
    virtual GCTAInstDir           mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                omega(void) const;
    virtual void                  read(const GXmlElement& xml);
    virtual void                  write(GXmlElement& xml) const;
    virtual std::string           print(void) const;

    // Other methods
    double sigma(void) const { return m_sigma.real_value(); }
    void   sigma(const double& sigma) { m_sigma.real_value(sigma); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadialGauss& model);
    void free_members(void);

    // Radial integration class (used by omega() method). Note that the
    // sigma parameter is given in rad^2
    class integrand : public GIntegrand {
    public:
        integrand(double sigma) : m_sigma(sigma) { }
        double eval(double x) {
            double arg  = x * x / m_sigma;
            double arg2 = arg * arg;
            double f    = exp(-0.5 * arg2);
            return (f*sin(x));
        }
    private:
        double m_sigma;
    };

    // Protected members
    GModelPar m_sigma;        //!< Width parameter (degrees^2)
};

#endif /* GCTAMODELRADIALGAUSS_HPP */
