/***************************************************************************
 *     GModelSpectralExpPlaw.hpp - Exponential cut off power law model     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralExpPlaw.hpp
 * @brief Exponential cut off power law spectral class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALEXPPLAW_HPP
#define GMODELSPECTRALEXPPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GModelSpectralExpPlaw
 *
 * @brief Exponential cut off power law spectral class
 *
 * This class implements a power law as the spectral component of the
 * gamma-ray sky model. The power law is defined as
 * \f[I(E)=norm (E/pivot)^{index} \exp(-E/ecut)\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy,
 * \f$index\f$ is the spectral index, and
 * \f$ecut\f$ is the cut off energy.
 ***************************************************************************/
class GModelSpectralExpPlaw : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralExpPlaw(void);
    explicit GModelSpectralExpPlaw(const double& norm, const double& index, const double& ecut);
    explicit GModelSpectralExpPlaw(const GXmlElement& xml);
    GModelSpectralExpPlaw(const GModelSpectralExpPlaw& model);
    virtual ~GModelSpectralExpPlaw(void);

    // Operators
    virtual GModelSpectralExpPlaw& operator= (const GModelSpectralExpPlaw& model);

    // Implemented pure virtual methods
    virtual void                   clear(void);
    virtual GModelSpectralExpPlaw* clone(void) const;
    virtual std::string            type(void) const { return "ExpCutoff"; }
    virtual double                 eval(const GEnergy& srcEng) const;
    virtual double                 eval_gradients(const GEnergy& srcEng) const;
    virtual double                 flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double                 eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy                mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                   read(const GXmlElement& xml);
    virtual void                   write(GXmlElement& xml) const;
    virtual std::string            print(void) const;

    // Other methods
    double norm(void) const { return m_norm.Value(); }
    double index(void) const { return m_index.Value(); }
    double ecut(void) const { return m_ecut.Value(); }
    double pivot(void) const { return m_pivot.Value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralExpPlaw& model);
    void free_members(void);
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Photon flux integration kernel
    class flux_kernel : public GFunction {
    public:
        flux_kernel(const double& norm,
                    const double& index,
                    const double& pivot,
                    const double& ecut) :
                    m_norm(norm),
                    m_index(index),
                    m_pivot(pivot),
                    m_ecut(ecut) {}
        double eval(double eng);
    protected:
        double m_norm;   //!< Normalization
        double m_index;  //!< Index
        double m_pivot;  //!< Pivot energy
        double m_ecut;   //!< Cut off energy
    };

    // Energy flux integration kernel
    class eflux_kernel : public GFunction {
    public:
        eflux_kernel(const double& norm,
                     const double& index,
                     const double& pivot,
                     const double& ecut) :
                     m_norm(norm),
                     m_index(index),
                     m_pivot(pivot),
                     m_ecut(ecut) {}
        double eval(double eng);
    protected:
        double m_norm;   //!< Normalization
        double m_index;  //!< Index
        double m_pivot;  //!< Pivot energy
        double m_ecut;   //!< Cut off energy
    };

    // Protected members
    GModelPar m_norm;               //!< Normalization factor
    GModelPar m_index;              //!< Spectral index
    GModelPar m_ecut;               //!< Exponential cut off energy
    GModelPar m_pivot;              //!< Pivot energy

    // Cached members used for pre-computations
    mutable double m_mc_emin;       //!< Minimum energy
    mutable double m_mc_emax;       //!< Maximum energy
    mutable double m_mc_exponent;   //!< Exponent (index+1)
    mutable double m_mc_pow_emin;   //!< Power of minimum energy
    mutable double m_mc_pow_ewidth; //!< Power of energy width
};

#endif /* GMODELSPECTRALEXPPLAW_HPP */
