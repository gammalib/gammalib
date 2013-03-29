/***************************************************************************
 *         GModelSpectralPlaw.hpp - Spectral power law model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlaw.hpp
 * @brief Power law spectral model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALPLAW_HPP
#define GMODELSPECTRALPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralPlaw
 *
 * @brief Power law spectral model class
 *
 * This class implements a power law spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_index}\f$ is the spectral index, and
 * - \f${\tt m\_pivot}\f$ is the pivot energy.
 ***************************************************************************/
class GModelSpectralPlaw : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlaw(void);
    explicit GModelSpectralPlaw(const double& norm, const double& index,
                                const GEnergy& pivot);
    explicit GModelSpectralPlaw(const GXmlElement& xml);
    GModelSpectralPlaw(const GModelSpectralPlaw& model);
    virtual ~GModelSpectralPlaw(void);

    // Operators
    virtual GModelSpectralPlaw& operator=(const GModelSpectralPlaw& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralPlaw* clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime&   srcTime) const;
    virtual double              eval_gradients(const GEnergy& srcEng,
                                               const GTime&   srcTime);
    virtual double              flux(const GEnergy& emin,
                                     const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin,
                                   const GEnergy& emax,
                                   const GTime&   time,
                                   GRan&          ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double  norm(void) const;
    double  index(void) const;
    GEnergy pivot(void) const;
    void    norm(const double& norm);
    void    index(const double& index);
    void    pivot(const GEnergy& pivot);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralPlaw& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    GModelPar m_norm;                //!< Normalization factor
    GModelPar m_index;               //!< Spectral index
    GModelPar m_pivot;               //!< Pivot energy

    // Cached members used for pre-computations
    mutable GEnergy m_last_energy;     //!< Last energy value
    mutable double  m_last_index;      //!< Last index parameter
    mutable double  m_last_pivot;      //!< Last pivot parameter
    mutable double  m_last_e_norm;     //!< Last E/Epivot value
    mutable double  m_last_log_e_norm; //!< Last ln(E/Epivot) value
    mutable double  m_last_power;      //!< Last power value
    mutable double  m_mc_emin;         //!< Minimum energy
    mutable double  m_mc_emax;         //!< Maximum energy
    mutable double  m_mc_exponent;     //!< Exponent (index+1)
    mutable double  m_mc_pow_emin;     //!< Power of minimum energy
    mutable double  m_mc_pow_ewidth;   //!< Power of energy width
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "PowerLaw".
 *
 * Returns the type of the spectral power law model.
 ***************************************************************************/
inline
std::string GModelSpectralPlaw::type(void) const
{
    return "PowerLaw";
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor (ph/cm2/s/MeV).
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelSpectralPlaw::norm(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set normalization factor 
 *
 * @param[in] norm Normalization factor (ph/cm2/s/MeV).
 *
 * Sets the normalization factor.
 ***************************************************************************/
inline
void GModelSpectralPlaw::norm(const double& norm)
{
    m_norm.value(norm);
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
double GModelSpectralPlaw::index(void) const
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
void GModelSpectralPlaw::index(const double& index)
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
GEnergy GModelSpectralPlaw::pivot(void) const
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
void GModelSpectralPlaw::pivot(const GEnergy& pivot)
{
    m_pivot.value(pivot.MeV());
    return;
}

#endif /* GMODELSPECTRALPLAW_HPP */
