/***************************************************************************
 *     GModelSpectralBrokenPlaw.hpp - Broken power law spectrum class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Anneli Schulz                               *
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
 * @file GModelSpectralBrokenPlaw.hpp
 * @brief Broken power law spectrum class definition
 * @author Anneli Schulz
 */

#ifndef GMODELSPECTRALBROKENPLAW_HPP
#define GMODELSPECTRALBROKENPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralBrokenPlaw
 *
 * @brief Broken power law spectral model class
 *
 * This class implements a broken power law spectrum. The model is defined
 * by
 *
 * \f[
 *    S_{\rm E}(E | t) = k_0 \times \left \{
 *    \begin{eqnarray}
 *     \left( \frac{E}{E_b} \right)^{\gamma_1} & {\rm if\,\,} E < E_b \\
 *     \left( \frac{E}{E_b} \right)^{\gamma_2} & {\rm otherwise}
 *    \end{eqnarray}
 *    \right .
 * \f]
 *
 * where
 * \f$k_0\f$ is the normalization or prefactor,
 * \f$\gamma_1\f$ is the spectral index before the break,
 * \f$\gamma_2\f$ is the spectral index after the break, and
 * \f$E_b\f$ is the break energy.
 ***************************************************************************/
class GModelSpectralBrokenPlaw : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralBrokenPlaw(void);
    GModelSpectralBrokenPlaw(const std::string& type,
                             const std::string& prefactor,
                             const std::string& index1,
                             const std::string& breakenergy,
                             const std::string& index2);
    GModelSpectralBrokenPlaw(const double&  prefactor,
                             const double&  index1,
                             const GEnergy& breakenergy,
                             const double&  index2);
    explicit GModelSpectralBrokenPlaw(const GXmlElement& xml);
    GModelSpectralBrokenPlaw(const GModelSpectralBrokenPlaw& model);
    virtual ~GModelSpectralBrokenPlaw(void);

    // Operators
    virtual GModelSpectralBrokenPlaw& operator=(const GModelSpectralBrokenPlaw& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpectralBrokenPlaw* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GEnergy& srcEng,
                                           const GTime& srcTime = GTime(),
                                           const bool& gradients = false) const;
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
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  prefactor(void) const;
    double  index1(void) const;
    double  index2(void) const;
    GEnergy breakenergy(void) const;
    void    prefactor(const double& prefactor);
    void    index1(const double& index1);
    void    index2(const double& index2);
    void    breakenergy(const GEnergy& breakenergy);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralBrokenPlaw& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    std::string m_type;                   //!< Model type
    GModelPar   m_norm;                   //!< Normalization factor
    GModelPar   m_index1;                 //!< Spectral index1
    GModelPar   m_index2;                 //!< Spectral index2
    GModelPar   m_breakenergy;            //!< Energy of spectral break

    // Cached members used for pre-computations
    mutable GEnergy m_last_energy;        //!< Last energy value
    mutable double  m_last_index1;        //!< Last index1 parameter
    mutable double  m_last_index2;        //!< Last index1 parameter
    mutable double  m_last_breakenergy;   //!< Last breakenergy parameter
    mutable double  m_last_e_norm;        //!< Last E/Ebreakenergy value
    mutable double  m_last_log_e_norm;    //!< Last ln(E/Ebreakenergy) value
    mutable double  m_last_power;         //!< Last power value
    mutable double  m_mc_emin;            //!< Minimum energy
    mutable double  m_mc_emax;            //!< Maximum energy
    mutable double  m_mc_exponent1;       //!< Exponent (index1+1)
    mutable double  m_mc_exponent2;       //!< Exponent (index2+1)
    mutable double  m_mc_pow_emin;        //!< Power of minimum energy
    mutable double  m_mc_pow_ewidth;      //!< Power of energy width
    mutable std::vector<double> m_mc_cum; //!< Cumulative distribution
    mutable std::vector<double> m_mc_min; //!< Lower boundary for MC
    mutable std::vector<double> m_mc_max; //!< Upper boundary for MC
    mutable std::vector<double> m_mc_exp; //!< Exponent for MC
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralBrokenPlaw").
 ***************************************************************************/
inline
std::string GModelSpectralBrokenPlaw::classname(void) const
{
    return ("GModelSpectralBrokenPlaw");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "PowerLaw".
 *
 * Returns the type of the spectral broken power law model.
 ***************************************************************************/
inline
std::string GModelSpectralBrokenPlaw::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return pre factor
 *
 * @return Pre factor (ph/cm2/s/MeV).
 *
 * Returns the pre factor.
 ***************************************************************************/
inline
double GModelSpectralBrokenPlaw::prefactor(void) const
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
void GModelSpectralBrokenPlaw::prefactor(const double& prefactor)
{
    m_norm.value(prefactor);
    return;
}


/***********************************************************************//**
 * @brief Return power law index1
 *
 * @return Power law index1.
 *
 * Returns the power law index1.
 ***************************************************************************/
inline
double GModelSpectralBrokenPlaw::index1(void) const
{
    return (m_index1.value());
}


/***********************************************************************//**
 * @brief Set power law index1
 *
 * @param[in] index1 Power law index1.
 *
 * Sets the power law index1.
 ***************************************************************************/
inline
void GModelSpectralBrokenPlaw::index1(const double& index1)
{
    m_index1.value(index1);
    return;
}


/***********************************************************************//**
 * @brief Return power law index2
 *
 * @return Power law index2.
 *
 * Returns the power law index2.
 ***************************************************************************/
inline
double GModelSpectralBrokenPlaw::index2(void) const
{
    return (m_index2.value());
}


/***********************************************************************//**
 * @brief Set power law index2
 *
 * @param[in] index2 Power law index2.
 *
 * Sets the power law index2.
 ***************************************************************************/
inline
void GModelSpectralBrokenPlaw::index2(const double& index2)
{
    m_index2.value(index2);
    return;
}


/***********************************************************************//**
 * @brief Return breakenergy energy
 *
 * @return breakenergy energy.
 *
 * Returns the breakenergy energy.
 ***************************************************************************/
inline
GEnergy GModelSpectralBrokenPlaw::breakenergy(void) const
{
    GEnergy energy;
    energy.MeV(m_breakenergy.value());
    return energy;
}


/***********************************************************************//**
 * @brief Set breakenergy energy
 *
 * @param[in] breakenergy breakenergy energy.
 *
 * Sets the breakenergy energy.
 ***************************************************************************/
inline
void GModelSpectralBrokenPlaw::breakenergy(const GEnergy& breakenergy)
{
    m_breakenergy.value(breakenergy.MeV());
    return;
}

#endif /* GMODELSPECTRALPLAW_HPP */
