/***************************************************************************
 *         GModelSpectralGauss.hpp - Spectral Gaussian model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Christoph Deil & Ellis Owen                      *
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
 * @file GModelSpectralGauss.hpp
 * @brief Gaussian spectral model class interface definition
 * @author Christoph Deil & Ellis Owen
 */

#ifndef GMODELSPECTRALGAUSS_HPP
#define GMODELSPECTRALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GModelSpectralGauss
 *
 * @brief Gaussian spectral model class
 *
 * This class implements a Gaussian spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) =
 *    \frac{\tt m\_norm}{\sqrt{2 \pi} m\_sigma^2} \exp
 *    \left(\frac{-(E - m\_mean)^2}{2 m_\sigma^2} \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization,
 * - \f${\tt m\_mean}\f$ is the mean energy,
 * - \f${\tt m\_sigma}\f$ is the energy width.
 ***************************************************************************/
class GModelSpectralGauss : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralGauss(void);
    explicit GModelSpectralGauss(const GXmlElement& xml);
    GModelSpectralGauss(const double&  norm,
                        const GEnergy& mean,
                        const GEnergy& sigma);
    GModelSpectralGauss(const GModelSpectralGauss& model);
    virtual ~GModelSpectralGauss(void);

    // Operators
    virtual GModelSpectralGauss& operator=(const GModelSpectralGauss& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralGauss* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime&   srcTime);
    virtual double               flux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin,
                                       const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin,
                                    const GEnergy& emax,
                                    const GTime&   time,
                                    GRan&          ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  norm(void) const;
    void    norm(const double& norm);
    GEnergy mean(void) const;
    void    mean(const GEnergy& mean);
    GEnergy sigma(void) const;
    void    sigma(const GEnergy& sigma);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralGauss& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;

    // Energy flux integration kernel
    class eflux_kernel : public GFunction {
    public:
        eflux_kernel(const double& norm,
                     const double& mean,
                     const double& sigma) :
                     m_norm(norm),
                     m_mean(mean),
                     m_sigma(sigma) {}
        double eval(const double& eng);
    protected:
        const double& m_norm;      //!< Normalization
        const double& m_mean;      //!< Mean
        const double& m_sigma;     //!< Sigma
    };

    // Protected members
    GModelPar m_norm;  //!< Normalization factor
    GModelPar m_mean;  //!< Gaussian mean energy
    GModelPar m_sigma; //!< Gaussian energy width
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "ConstantValue".
 *
 * Returns the type of the Gaussian spectral model.
 ***************************************************************************/
inline
std::string GModelSpectralGauss::type(void) const
{
    return "Gaussian";
}


/***********************************************************************//**
 * @brief Return normalization
 *
 * @return Normalization (ph/cm2/s).
 ***************************************************************************/
inline
double GModelSpectralGauss::norm(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set normalization
 *
 * @param[in] norm Normalization (ph/cm2/s).
 ***************************************************************************/
inline
void GModelSpectralGauss::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}

/***********************************************************************//**
 * @brief Return mean energy
 *
 * @return Mean energy.
 ***************************************************************************/
inline
GEnergy GModelSpectralGauss::mean(void) const
{
  GEnergy energy;
  energy.MeV(m_mean.value());
  return energy;
}


/***********************************************************************//**
 * @brief Set mean energy
 *
 * @param[in] Mean energy.
 ***************************************************************************/
inline
void GModelSpectralGauss::mean(const GEnergy& mean)
{
    m_mean.value(mean.MeV());
    return;
}

/***********************************************************************//**
 * @brief Return energy width
 *
 * @return Energy width
 ***************************************************************************/
inline
GEnergy GModelSpectralGauss::sigma(void) const
{
  GEnergy energy;
  energy.MeV(m_sigma.value());
  return energy;
}


/***********************************************************************//**
 * @brief Set energy width
 *
 * @param[in] Energy width
 ***************************************************************************/
inline
void GModelSpectralGauss::sigma(const GEnergy& sigma)
{
    m_sigma.value(sigma.MeV());
    return;
}

#endif /* GMODELSPECTRALGAUSS_HPP */
