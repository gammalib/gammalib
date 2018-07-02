/***************************************************************************
 *  GModelSpectralExponential.hpp - Exponential spectral model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Luigi Tibaldo                                    *
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
 * @file GModelSpectralExponential.hpp
 * @brief Exponential spectral model class interface definition
 * @author Luigi Tibaldo
 */

#ifndef GMODELSPECTRALEXPONENTIAL_HPP
#define GMODELSPECTRALEXPONENTIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GModelSpectralNodes.hpp"
#include "GFunction.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralExponential
 *
 * @brief Exponential spectral model class
 *
 * This class implements a spectral model that is the exponential of a
 * spectral model component. The spectral model component can be defined in
 * an XML file, or using the exponent() method.
 ***************************************************************************/
class GModelSpectralExponential : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralExponential(void);
    explicit GModelSpectralExponential(const GXmlElement& xml);
    explicit GModelSpectralExponential(const GModelSpectral* spec);
    GModelSpectralExponential(const GModelSpectralExponential& model);
    virtual ~GModelSpectralExponential(void);

    // Operators
    virtual GModelSpectralExponential& operator=(const GModelSpectralExponential& model);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GModelSpectralExponential* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GEnergy& srcEng,
                                            const GTime&   srcTime = GTime(),
                                            const bool&    gradients = false) const;
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
    virtual std::string                print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                  exponent(const GModelSpectral* spec);
    const GModelSpectral* exponent(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralExponential& model);
    void free_members(void);
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Class to determine the integral photon flux
    class flux_kern : public GFunction {
    public:
        // Constructor
        flux_kern(const GModelSpectral* spec) : m_exp(spec) {}

        // Method
        double eval(const double& x) {
            GEnergy energy(x, "MeV");
            double value = m_exp->eval(energy);
            value = exp(value);
            return value;
        }
    protected:
        const GModelSpectral* m_exp;
    };

    // Class to determine the integral energy flux, derived from flux_kern
    class eflux_kern : public flux_kern {
    public:
        // Constructor
        eflux_kern(const GModelSpectral* spec) : flux_kern(spec) {}

        // Method
        double eval(const double& x) {
            return x * flux_kern::eval(x);
        }
    };

    // Protected members
    std::string					 m_type;        //!< Model type
    GModelSpectral*				 m_exponent;	//!< Exponent model component

    // MC cache
    mutable GModelSpectralNodes  m_mc_spectrum; //!< MC spectrum cache
    mutable GEnergy              m_mc_emin;     //!< Last minimum energy
    mutable GEnergy              m_mc_emax;     //!< Last maximum energy
    mutable std::vector<double>  m_mc_values;   //!< Parameter values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralExponential").
 ***************************************************************************/
inline
std::string GModelSpectralExponential::classname(void) const
{
    return ("GModelSpectralExponential");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spectral Multiplicative model.
 ***************************************************************************/
inline
std::string GModelSpectralExponential::type(void) const
{
    return (m_type);
}

#endif /* GMODELSPECTRALEXPONENTIAL_HPP */
