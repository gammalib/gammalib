/***************************************************************************
 *     GModelSpectralExpPlaw2.hpp - Exponential cut off power law model    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Alexander Ziegler                                *
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
 * @file GModelSpectralExpPlaw2.hpp
 * @brief Exponential cut off power law spectral class interface definition
 * @author Alexander Ziegler
 */

#ifndef GMODELSPECTRALEXPPLAW2_HPP
#define GMODELSPECTRALEXPPLAW2_HPP

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


/***********************************************************************//**
 * @class GModelSpectralExpPlaw2
 *
 * @brief Exponential cut off power law spectral class
 *
 * This class implements a power law spectrum with exponential cut off. The
 * model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index}
 *    \exp \left( -{\tt m\_lamda}*E \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_pivot}\f$ is the pivot energy,
 * - \f${\tt m\_index}\f$ is the spectral index, and
 * - \f${\tt m\_lambda}\f$ is the cut-off parameter (1/cut-off energy).
 ***************************************************************************/
class GModelSpectralExpPlaw2 : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralExpPlaw2(void);
    explicit GModelSpectralExpPlaw2(const double&  prefactor,
                                    const double&  index,
                                    const GEnergy& pivot,
                                    const double& lambda);
    explicit GModelSpectralExpPlaw2(const GXmlElement& xml);
    GModelSpectralExpPlaw2(const GModelSpectralExpPlaw2& model);
    virtual ~GModelSpectralExpPlaw2(void);

    // Operators
    virtual GModelSpectralExpPlaw2& operator=(const GModelSpectralExpPlaw2& model);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GModelSpectralExpPlaw2* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const GEnergy& srcEng,
                                         const GTime&   srcTime = GTime()) const;
    virtual double                  eval_gradients(const GEnergy& srcEng,
                                                  const GTime&   srcTime = GTime());
    virtual double                  flux(const GEnergy& emin,
                                         const GEnergy& emax) const;
    virtual double                  eflux(const GEnergy& emin,
                                          const GEnergy& emax) const;
    virtual GEnergy                 mc(const GEnergy& emin,
                                       const GEnergy& emax,
                                       const GTime&   time,
                                       GRan&          ran) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;
    virtual std::string             print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index(void) const;
    void    index(const double& index);
    double  lambda(void) const;
    void    lambda(const double& lambda);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralExpPlaw2& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Photon flux integration kernel
    class flux_kernel : public GFunction {
    public:
        flux_kernel(const double& norm,
                    const double& index,
                    const double& pivot,
                    const double& lambda) :
                    m_norm(norm),
                    m_index(index),
                    m_inv_pivot(1.0/pivot),
                    m_lambda(lambda) {}
        double eval(const double& eng);
    protected:
        double m_norm;      //!< Normalization
        double m_index;     //!< Index
        double m_inv_pivot; //!< 1 / Pivot energy
        double m_lambda;    //!< Cut-off parameter
    };

    // Energy flux integration kernel
    class eflux_kernel : public GFunction {
    public:
        eflux_kernel(const double& norm,
                     const double& index,
                     const double& pivot,
                     const double& lambda) :
                     m_norm(norm),
                     m_index(index),
                     m_inv_pivot(1.0/pivot),
                     m_lambda(lambda) {}
        double eval(const double& eng);
    protected:
        double m_norm;      //!< Normalization
        double m_index;     //!< Index
        double m_inv_pivot; //!< 1 / Pivot energy
        double m_lambda;    //!< Cut-off parameter
    };

    // Protected members
    GModelPar m_norm;               //!< Normalization factor
    GModelPar m_index;              //!< Spectral index
    GModelPar m_lambda;             //!< Cut-off parameter
    GModelPar m_pivot;              //!< Pivot energy

    // Cached members used for pre-computations
    mutable GEnergy m_last_energy;   //!< Last energy value
    mutable double  m_last_index;    //!< Last index parameter
    mutable double  m_last_lambda;   //!< Last cut-off parameter
    mutable double  m_last_pivot;    //!< Last pivot parameter
    mutable double  m_last_e_norm;   //!< Last E/Epivot value
    mutable double  m_last_e_lambda; //!< Last E*lambda value
    mutable double  m_last_power;    //!< Last power value
    mutable double  m_mc_emin;       //!< Minimum energy
    mutable double  m_mc_emax;       //!< Maximum energy
    mutable double  m_mc_exponent;   //!< Exponent (index+1)
    mutable double  m_mc_pow_emin;   //!< Power of minimum energy
    mutable double  m_mc_pow_ewidth; //!< Power of energy width
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralExpPlaw2").
 ***************************************************************************/
inline
std::string GModelSpectralExpPlaw2::classname(void) const
{
    return ("GModelSpectralExpPlaw2");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "ExpCutoff2".
 *
 * Returns the type of the exponentially cut off power law model.
 ***************************************************************************/
inline
std::string GModelSpectralExpPlaw2::type(void) const
{
    return "ExpCutoff2";
}


/***********************************************************************//**
 * @brief Return pre factor
 *
 * @return Pre factor (ph/cm2/s/MeV).
 *
 * Returns the pre factor.
 ***************************************************************************/
inline
double GModelSpectralExpPlaw2::prefactor(void) const
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
void GModelSpectralExpPlaw2::prefactor(const double& prefactor)
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
double GModelSpectralExpPlaw2::index(void) const
{
    return (m_index.value());
}


/***********************************************************************//**
 * @brief Set power law index 
 *
 * @param[in] index Power law index.
 *
 * Sets the power law index.
 ***************************************************************************/
inline
void GModelSpectralExpPlaw2::index(const double& index)
{
    m_index.value(index);
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
GEnergy GModelSpectralExpPlaw2::pivot(void) const
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
void GModelSpectralExpPlaw2::pivot(const GEnergy& pivot)
{
    m_pivot.value(pivot.MeV());
    return;
}


/***********************************************************************//**
 * @brief Return exponential cut-off parameter
 *
 * @return Exponential cut-off parameter.
 *
 * Returns the exponential cut-off parameter.
 ***************************************************************************/
inline
double GModelSpectralExpPlaw2::lambda(void) const
{
	return (m_lambda.value());
}


/***********************************************************************//**
 * @brief Set exponential cut-off parameter
 *
 * @param[in] lambda Exponential cut-off parameter.
 *
 * Sets the exponential cut-off parameter.
 ***************************************************************************/
inline
void GModelSpectralExpPlaw2::lambda(const double& lambda)
{
    m_lambda.value(lambda);
    return;
}

#endif /* GMODELSPECTRALEXPPLAW2_HPP */
