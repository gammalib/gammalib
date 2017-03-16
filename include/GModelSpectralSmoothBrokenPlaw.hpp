/***************************************************************************
 *               GModelSpectralSmoothBrokenPlaw.cpp -                      *
 *             Smoothly broken power law spectrum class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Josh Cardenzana                                  *
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
 * @file GModelSpectralSmoothBrokenPlaw.hpp
 * @brief Smoothly broken power law spectrum class definition
 * @author Josh Cardenzana
 */

#ifndef GMODELSPECTRALSMOOTHBROKENPLAW_HPP
#define GMODELSPECTRALSMOOTHBROKENPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpectral.hpp"
#include "GFunction.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;

// TODO: Updated class documentation
/***********************************************************************//**
 * @class GModelSpectralSmoothBrokenPlaw
 *
 * @brief Smoothly broken power law spectral model class
 *
 * This class implements a smoothly broken power law spectrum. The model is 
 * defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = k_0 \times \left \{
 *     \left( \frac{E}{E_b} \right)^{\gamma_2} & {\rm otherwise}
 *    \right .
 * \f]
 *
 * where
 * \f$k_0\f$ is the normalization or prefactor,
 * \f$\gamma_1\f$ is the spectral index before the break,
 * \f$\gamma_2\f$ is the spectral index after the break, and
 * \f$E_b\f$ is the break energy.
 * \f$\beta\f$ defines the break smoothness
 ***************************************************************************/
class GModelSpectralSmoothBrokenPlaw : public GModelSpectral {
    
public:
    // Constructors and destructors
    GModelSpectralSmoothBrokenPlaw(void);
    GModelSpectralSmoothBrokenPlaw(const std::string& type,
                                   const std::string& prefactor,
                                   const std::string& index1,
                                   const std::string& pivot,
                                   const std::string& index2,
                                   const std::string& breakenergy,
                                   const std::string& beta);
    GModelSpectralSmoothBrokenPlaw(const double&  prefactor,
                                   const double&  index1,
                                   const GEnergy& pivot,
                                   const double&  index2,
                                   const GEnergy& breakenergy,
                                   const double&  beta);
    explicit GModelSpectralSmoothBrokenPlaw(const GXmlElement& xml);
    GModelSpectralSmoothBrokenPlaw(const GModelSpectralSmoothBrokenPlaw& model);
    virtual ~GModelSpectralSmoothBrokenPlaw(void);
    
    // Operators
    virtual GModelSpectralSmoothBrokenPlaw& operator=(const GModelSpectralSmoothBrokenPlaw& model);
    
    // Implemented pure virtual methods
    virtual void                            clear(void);
    virtual GModelSpectralSmoothBrokenPlaw* clone(void) const;
    virtual std::string                     classname(void) const;
    virtual std::string                     type(void) const;
    virtual double                          eval(const GEnergy& srcEng,
                                                 const GTime& srcTime = GTime(),
                                                 const bool& gradients = false) const;
    virtual double                          flux(const GEnergy& emin,
                                                 const GEnergy& emax) const;
    virtual double                          eflux(const GEnergy& emin,
                                                  const GEnergy& emax) const;
    virtual GEnergy                         mc(const GEnergy& emin,
                                               const GEnergy& emax,
                                               const GTime&   time,
                                               GRan&          ran) const;
    virtual void                            read(const GXmlElement& xml);
    virtual void                            write(GXmlElement& xml) const;
    virtual std::string                     print(const GChatter& chatter = NORMAL) const;
    
    // Other methods
    // Methods for getting the current values
    double  prefactor(void) const;
    double  index1(void) const;
    double  index2(void) const;
    GEnergy pivot(void) const;
    GEnergy breakenergy(void) const;
    double  beta(void) const;
    // Methods for setting the parameter values
    void    prefactor(const double& prefactor);
    void    index1(const double& index1);
    void    index2(const double& index2);
    void    pivot(const GEnergy& pivot);
    void    breakenergy(const GEnergy& breakenergy);
    void    beta(const double& beta);
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralSmoothBrokenPlaw& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;
    
    // Class to determine to the integral photon flux
    class flux_kern : public GFunction {
    public:
        // Constructor
        flux_kern(const double&  prefactor,
                  const double&  index1,
                  const GEnergy& pivot,
                  const double&  index2,
                  const GEnergy& breakenergy,
                  const double&  beta) :
            m_prefactor(prefactor),
            m_index1(index1),
            m_index2(index2),
            m_pivot(pivot),
            m_breakenergy(breakenergy),
            m_beta(beta)
        {};
        
        // Method
        double eval(const double& energy) {
            double epivot = energy / m_pivot.MeV();
            double ebreak = energy / m_breakenergy.MeV();
            return m_prefactor * std::pow(epivot, m_index1) *
            std::pow(1.0 + std::pow(ebreak,(m_index1-m_index2)/m_beta),-m_beta);
        }
        
        // Members
    protected:
        double  m_prefactor; //!< Normalization
        double  m_index1;	 //!< Spectral index1
        double  m_index2;    //!< Spectral index2
        GEnergy m_pivot;	 //!< Pivot energy
        GEnergy m_breakenergy; //!< Break energy
        double  m_beta;      //!< Break smoothness parameter
    };
    
    // Class to determine the integral energy flux, derivation of flux_kern
    class eflux_kern : public flux_kern {
    public:
        // Constructor
        eflux_kern(const double&  prefactor,
                   const double&  index1,
                   const GEnergy& pivot,
                   const double&  index2,
                   const GEnergy& breakenergy,
                   const double&  beta):
            flux_kern(prefactor, index1, pivot, index2, breakenergy, beta)
        {};
        
        // Method
        double eval(const double& energy) {
            return energy * flux_kern::eval(energy);
        }
    };
    
    // Protected members
    std::string m_type;                   //!< Model type
    GModelPar   m_norm;                   //!< Normalization factor
    GModelPar   m_index1;                 //!< Spectral index1
    GModelPar   m_index2;                 //!< Spectral index2
    GModelPar   m_pivot;                  //!< Pivot energy
    GModelPar   m_breakenergy;            //!< Energy of spectral break
    GModelPar   m_beta;                   //!< Break smoothness
    
    // Cached members used for pre-computations
    mutable GEnergy m_last_energy;        //!< Last energy value
    mutable double  m_last_index1;        //!< Last index1 parameter
    mutable double  m_last_index2;        //!< Last index2 parameter
    mutable double  m_last_pivot;         //!< Last pivot parameter
    mutable double  m_last_breakenergy;   //!< Last breakenergy parameter
    mutable double  m_last_beta;          //!< Last beta parameter
    mutable double  m_last_epivot_norm;   //!< Last E/Epivot value
    mutable double  m_last_ebreak_norm;   //!< Last E/Ebreakenergy value
    mutable double  m_last_log_epivot_norm; //!< Last ln(E/Epivot) value
    mutable double  m_last_log_ebreak_norm; //!< Last ln(E/Ebreakenergy) value
    mutable double  m_last_epivot_pow;    //!< Last pow(E/Epivot,index1) value
    mutable double  m_last_ebreak_pow;    //!< Last pow(E/Ebreakenergy,(index1-index2)/beta)

    mutable double  m_mc_emin;            //!< Minimum energy
    mutable double  m_mc_emax;            //!< Maximum energy
    mutable double  m_mc_plaw_prefactor;  //!< Prefactor for comparison p-laws
    mutable double  m_mc_exponentS;       //!< Exponent (index+1) for softer index
    mutable double  m_mc_exponentH;       //!< Exponent (index+1) for harder index
    mutable double  m_mc_pow_ewidth_low;  //!< width of energy range below break
    mutable double  m_mc_norm;            //!< Normalization term for energy generation
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralSmoothBrokenPlaw").
 ***************************************************************************/
inline
std::string GModelSpectralSmoothBrokenPlaw::classname(void) const
{
    return ("GModelSpectralSmoothBrokenPlaw");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "PowerLaw".
 *
 * Returns the type of the spectral smoothly broken power law model.
 ***************************************************************************/
inline
std::string GModelSpectralSmoothBrokenPlaw::type(void) const
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
double GModelSpectralSmoothBrokenPlaw::prefactor(void) const
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
void GModelSpectralSmoothBrokenPlaw::prefactor(const double& prefactor)
{
    m_norm.value(prefactor);
    return;
}


/***********************************************************************//**
 * @brief Return smoothly broken power law index1
 *
 * @return Power law index1.
 *
 * Returns the power law index1.
 ***************************************************************************/
inline
double GModelSpectralSmoothBrokenPlaw::index1(void) const
{
    return (m_index1.value());
}


/***********************************************************************//**
 * @brief Set smoothly broken power law index1
 *
 * @param[in] index1 Power law index1.
 *
 * Sets the power law index1.
 ***************************************************************************/
inline
void GModelSpectralSmoothBrokenPlaw::index1(const double& index1)
{
    m_index1.value(index1);
    return;
}


/***********************************************************************//**
 * @brief Return smoothly broken power law index2
 *
 * @return Power law index2.
 *
 * Returns the power law index2.
 ***************************************************************************/
inline
double GModelSpectralSmoothBrokenPlaw::index2(void) const
{
    return (m_index2.value());
}


/***********************************************************************//**
 * @brief Set smoothly broken power law index2
 *
 * @param[in] index2 Power law index2.
 *
 * Sets the power law index2.
 ***************************************************************************/
inline
void GModelSpectralSmoothBrokenPlaw::index2(const double& index2)
{
    m_index2.value(index2);
    return;
}


/***********************************************************************//**
* @brief Return pivot energy
*
* @return Smoothly broken power law pivot energy.
*
* Returns the smoothly broken power law scale energy.
***************************************************************************/
inline
GEnergy GModelSpectralSmoothBrokenPlaw::pivot(void) const
{
    GEnergy energy;
    energy.MeV(m_pivot.value());
    return energy;
}


/***********************************************************************//**
* @brief Set pivot energy
*
* @param[in] pivot Smoothly broken power law pivot energy.
*
* Sets the power law index2.
***************************************************************************/
inline
void GModelSpectralSmoothBrokenPlaw::pivot(const GEnergy& pivot)
{
    m_pivot.value(pivot.MeV());
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
GEnergy GModelSpectralSmoothBrokenPlaw::breakenergy(void) const
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
void GModelSpectralSmoothBrokenPlaw::breakenergy(const GEnergy& breakenergy)
{
    m_breakenergy.value(breakenergy.MeV());
    return;
}


/***********************************************************************//**
* @brief Returns break smoothness parameter beta
*
* @return break smoothness parameter.
*
* Returns the break smoothness parameter 'beta'.
***************************************************************************/
inline
double GModelSpectralSmoothBrokenPlaw::beta(void) const
{
    return (m_beta.value());
}


/***********************************************************************//**
* @brief Set break smoothness
*
* @param[in] beta break smoothness parameter.
*
* Sets the beta break smoothness parameter.
***************************************************************************/
inline
void GModelSpectralSmoothBrokenPlaw::beta(const double& beta)
{
    m_beta.value(beta);
    return;
}
#endif /* GMODELSPECTRALSMOOTHBROKENPLAW_HPP */
