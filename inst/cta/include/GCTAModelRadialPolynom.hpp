/***************************************************************************
 *      GCTAModelRadialPolynom.hpp  -  Radial Polynom CTA model class      *
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
 * @file GCTAModelRadialPolynom.hpp
 * @brief Radial Polynom model class interface definition
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRADIALPOLYNOM_HPP
#define GCTAMODELRADIALPOLYNOM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelPar.hpp"
#include "GXmlElement.hpp"
#include "GCTAModelRadial.hpp"
#include "GCTAInstDir.hpp"
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialPolynom
 *
 * @brief Radial Polynom CTA model class
 *
 * This class implements the radial function
 * \f[f(\theta) = \sum_{i=0}^m c_i \theta^i\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$c_i\f$ are the polynomial coefficients.
 *
 * This function represents a Polynom in \f$\theta\f$.
 ***************************************************************************/
class GCTAModelRadialPolynom : public GCTAModelRadial {

public:
    // Constructors and destructors
    GCTAModelRadialPolynom(void);
    explicit GCTAModelRadialPolynom(const std::vector<double>& coeffs);
    explicit GCTAModelRadialPolynom(const GXmlElement& xml);
    GCTAModelRadialPolynom(const GCTAModelRadialPolynom& model);
    virtual ~GCTAModelRadialPolynom(void);

    // Operators
    virtual GCTAModelRadialPolynom& operator= (const GCTAModelRadialPolynom& model);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelRadialPolynom* clone(void) const;
    virtual std::string             type(void) const { return "Polynom"; }
    virtual double                  eval(const double& offset) const;
    virtual double                  eval_gradients(const double& offset) const;
    virtual GCTAInstDir             mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                  omega(void) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;
    virtual std::string             print(void) const;

    // Other methods
    int    size(void) const { return m_coeffs.size(); }
    //double coeff(void) const;
    //void   coeff(const double& value);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadialPolynom& model);
    void free_members(void);
    void update_pars(void);

    // Radial integration class (used by omega() method). Note that the
    // integration is done in radians
    class integrand : public GIntegrand {
    public:
        integrand(const GCTAModelRadialPolynom* model) : m_model(model) { }
        double eval(double x) {
            return (sin(x)*m_model->eval(x*rad2deg));
        }
    private:
        const GCTAModelRadialPolynom* m_model;
    };

    // Protected members
    std::vector<GModelPar> m_coeffs;        //!< Coefficients
};

#endif /* GCTAMODELRADIALPOLYNOM_HPP */
