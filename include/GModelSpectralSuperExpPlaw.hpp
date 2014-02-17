/***************************************************************************
 *   GModelSpectralSuperExpPlaw.hpp - Super exp. cut off power law model   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                                    *
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
 * @file GModelSpectralSuperExpPlaw.hpp
 * @brief Super exponential cut off power law spectral class interface definition
 * @author Michael Mayer
 */

#ifndef GMODELSPECTRALSUPEREXPPLAW_HPP
#define GMODELSPECTRALSUPEREXPPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GModelSpectralSuperExpPlaw
 *
 * @brief Super exponential cut off power law spectral class
 *
 * This class implements a power law spectrum with super exponential cut off.
 * The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index1}
 *    \exp \left( \frac{-E}{\tt m\_ecut}^{\tt m\_index2} \right)
 * \f]
 *
 * where
 *
 *      \f${\tt m\_norm}\f$ is the normalization or prefactor,
 *      \f${\tt m\_index1}\f$ is the spectral index,
 *      \f${\tt m\_ecut}\f$ is the cut off energy, and
 *      \f${\tt m\_pivot}\f$ is the pivot energy.
 *      \f${\tt m\_index2}\f$ is index defining the cutoff shape,
 ***************************************************************************/
class GModelSpectralSuperExpPlaw : public GModelSpectral {

public:
    // Constructors and destructors
	GModelSpectralSuperExpPlaw(void);
    GModelSpectralSuperExpPlaw(const double&  prefactor,
                               const double&  index1,
                               const GEnergy& pivot,
                               const GEnergy& cutoff,
                               const double&  index2);
    explicit GModelSpectralSuperExpPlaw(const GXmlElement& xml);
    GModelSpectralSuperExpPlaw(const GModelSpectralSuperExpPlaw& model);
    virtual ~GModelSpectralSuperExpPlaw(void);

    // Operators
    virtual GModelSpectralSuperExpPlaw& operator=(const GModelSpectralSuperExpPlaw& model);

    // Implemented pure virtual methods
    virtual void                        clear(void);
    virtual GModelSpectralSuperExpPlaw* clone(void) const;
    virtual std::string                 type(void) const;
    virtual double                      eval(const GEnergy& srcEng,
                                             const GTime&   srcTime) const;
    virtual double                      eval_gradients(const GEnergy& srcEng,
                                                       const GTime&   srcTime);
    virtual double                      flux(const GEnergy& emin,
                                             const GEnergy& emax) const;
    virtual double                      eflux(const GEnergy& emin,
                                              const GEnergy& emax) const;
    virtual GEnergy                     mc(const GEnergy& emin,
                                           const GEnergy& emax,
                                           const GTime&   time,
                                           GRan&          ran) const;
    virtual void                        read(const GXmlElement& xml);
    virtual void                        write(GXmlElement& xml) const;
    virtual std::string                 print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index1(void) const;
    void    index1(const double& index1);
    double  index2(void) const;
    void    index2(const double& index2);
    GEnergy cutoff(void) const;
    void    cutoff(const GEnergy& cutoff);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralSuperExpPlaw& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Photon flux integration kernel
    class flux_kernel : public GFunction {
    public:
        flux_kernel(const double& norm,
                    const double& index1,
                    const double& pivot,
                    const double& ecut,
                    const double& index2) :
                    m_norm(norm),
                    m_index1(index1),
                    m_inv_pivot(1.0/pivot),
                    m_inv_ecut(1.0/ecut),
                    m_index2(index2) {}
        double eval(const double& eng);
    protected:
        const double& m_norm;      //!< Normalization
        const double& m_index1;    //!< Index1
        double        m_inv_pivot; //!< 1 / Pivot energy
        double        m_inv_ecut;  //!< 1 / Cut off energy
        const double& m_index2;    //!< Index2
    };

    // Energy flux integration kernel
    class eflux_kernel : public GFunction {
    public:
        eflux_kernel(const double& norm,
                     const double& index1,
                     const double& pivot,
                     const double& ecut,
                     const double& index2) :
                     m_norm(norm),
                     m_index1(index1),
                     m_inv_pivot(1.0/pivot),
                     m_inv_ecut(1.0/ecut),
                     m_index2(index2) {}
        double eval(const double& eng);
    protected:
        const double& m_norm;      //!< Normalization
        const double& m_index1;    //!< Index1
        double        m_inv_pivot; //!< 1 / Pivot energy
        double        m_inv_ecut;  //!< 1 / Cut off energy
        const double& m_index2;    //!< Index2
    };

    // Protected members
    GModelPar m_norm;    //!< Normalization factor
    GModelPar m_index1;  //!< Spectral index
    GModelPar m_index2;  //!< Index of cutoff
    GModelPar m_ecut;    //!< Exponential cut off energy
    GModelPar m_pivot;   //!< Pivot energy

    // Cached members used for pre-computations
    mutable GEnergy m_last_energy;   //!< Last energy value
    mutable double  m_last_index1;   //!< Last index1 parameter
    mutable double  m_last_ecut;     //!< Last energy cut-off parameter
    mutable double  m_last_pivot;    //!< Last pivot parameter
    mutable double  m_last_index2;   //!< Last index2 parameter
    mutable double  m_last_e_norm;   //!< Last E/Epivot value
    mutable double  m_last_e_cut;    //!< Last E/Ecut value
    mutable double  m_last_exponent; //!< last pow(E/Ecut,index2) value
    mutable double  m_last_power;    //!< Last power value
    mutable double  m_mc_emin;       //!< Minimum energy
    mutable double  m_mc_emax;       //!< Maximum energy
    mutable double  m_mc_exponent;   //!< Exponent (index+1)
    mutable double  m_mc_pow_emin;   //!< Power of minimum energy
    mutable double  m_mc_pow_ewidth; //!< Power of energy width
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "SuperExpCutoff".
 *
 * Returns the type of the exponentially cut off power law model.
 ***************************************************************************/
inline
std::string GModelSpectralSuperExpPlaw::type(void) const
{
    return "PLSuperExpCutoff";
}


/***********************************************************************//**
 * @brief Return pre factor
 *
 * @return Pre factor (ph/cm2/s/MeV).
 *
 * Returns the pre factor.
 ***************************************************************************/
inline
double GModelSpectralSuperExpPlaw::prefactor(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set pre factor 
 *
 * @param[in] prefactor Pre factor (ph/cm2/s/MeV).
 *
 * Sets the pre factor.
 ***************************************************************************/
inline
void GModelSpectralSuperExpPlaw::prefactor(const double& prefactor)
{
    m_norm.value(prefactor);
    return;
}


/***********************************************************************//**
 * @brief Return power law index
 *
 * @return Power law index.
 *
 * Returns the power law index.
 ***************************************************************************/
inline
double GModelSpectralSuperExpPlaw::index1(void) const
{
    return (m_index1.value());
}


/***********************************************************************//**
 * @brief Set power law index 
 *
 * @param[in] index Power law index.
 *
 * Sets the power law index.
 ***************************************************************************/
inline
void GModelSpectralSuperExpPlaw::index1(const double& index1)
{
    m_index1.value(index1);
    return;
}

/***********************************************************************//**
 * @brief Return cut off index
 *
 * @return cut off index.
 *
 * Returns the cut off index.
 ***************************************************************************/
inline
double GModelSpectralSuperExpPlaw::index2(void) const
{
    return (m_index2.value());
}


/***********************************************************************//**
 * @brief Set cut off index
 *
 * @param[in] index Cut off index.
 *
 * Sets the cut off index.
 ***************************************************************************/
inline
void GModelSpectralSuperExpPlaw::index2(const double& index2)
{
    m_index2.value(index2);
    return;
}


/***********************************************************************//**
 * @brief Return pivot energy
 *
 * @return Pivot energy.
 *
 * Returns the pivot energy.
 ***************************************************************************/
inline
GEnergy GModelSpectralSuperExpPlaw::pivot(void) const
{
    GEnergy energy;
    energy.MeV(m_pivot.value());
    return energy;
}


/***********************************************************************//**
 * @brief Set pivot energy
 *
 * @param[in] pivot Pivot energy.
 *
 * Sets the pivot energy.
 ***************************************************************************/
inline
void GModelSpectralSuperExpPlaw::pivot(const GEnergy& pivot)
{
    m_pivot.value(pivot.MeV());
    return;
}


/***********************************************************************//**
 * @brief Return exponential cut-off energy
 *
 * @return Exponential cut-off energy.
 *
 * Returns the exponential cut-off energy.
 ***************************************************************************/
inline
GEnergy GModelSpectralSuperExpPlaw::cutoff(void) const
{
    GEnergy energy;
    energy.MeV(m_ecut.value());
    return energy;
}


/***********************************************************************//**
 * @brief Set exponential cut-off energy
 *
 * @param[in] cutoff Exponential cut-off energy.
 *
 * Sets the exponential cut-off energy.
 ***************************************************************************/
inline
void GModelSpectralSuperExpPlaw::cutoff(const GEnergy& cutoff)
{
    m_ecut.value(cutoff.MeV());
    return;
}

#endif /* GMODELSPECTRALSUPEREXPPLAW_HPP */
