/***************************************************************************
 *    GModelSpectralLogParabola.hpp - Log parabola spectral model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Michael Mayer                               *
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
 * @file GModelSpectralLogParabola.hpp
 * @brief Log parabola spectral model class definition
 * @author Michael Mayer 
 */

#ifndef GMODELSPECTRALLOGPARABOLA_HPP
#define GMODELSPECTRALLOGPARABOLA_HPP


/* __ Includes ___________________________________________________________ */
#include <string>
#include <cmath>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GIntegral.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GModelSpectralLogParabola
 *
 * @brief LogParabola spectral model class
 *
 * This class implements a log parabola spectrum. The spectrum is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{{\tt m\_index} +
 *    {\tt m\_curvature} \, \ln \frac{E}{\tt m\_pivot}}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_index}\f$ is the spectral index,
 * - \f${\tt m\_curvature}\f$ is the spectral curvature, and
 * - \f${\tt m\_pivot}\f$ is the pivot energy.
 ***************************************************************************/
class GModelSpectralLogParabola : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralLogParabola(void);
    explicit GModelSpectralLogParabola(const double&  prefactor,
                                       const double&  index,
                                       const GEnergy& pivot,
                                       const double&  curvature);
    explicit GModelSpectralLogParabola(const GXmlElement& xml);
    GModelSpectralLogParabola(const GModelSpectralLogParabola& model);
    virtual ~GModelSpectralLogParabola(void);

    // Operators
    virtual GModelSpectralLogParabola& operator=(const GModelSpectralLogParabola& model);

    // Implemented pure virtual base class methods
    virtual void                       clear(void);
    virtual GModelSpectralLogParabola* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GEnergy& srcEng,
                                            const GTime&   srcTime) const;
    virtual double                     eval_gradients(const GEnergy& srcEng,
                                                      const GTime&   srcTime);
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
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index(void) const;
    void    index(const double& index);
    double  curvature(void) const;
    void    curvature(const double& curvature);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralLogParabola& model);
    void free_members(void);
    void update_eval_cache(const GEnergy& energy) const;
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax,
                         const GTime&   time) const;

    // Class to determine to the integral photon flux
    class flux_kern : public GFunction {
    public:
    	// Constructor
    	flux_kern(const double&  norm,
                  const double&  index,
                  const double&  curvature,
                  const GEnergy& pivot) :
    		      m_norm(norm),
                  m_index(index),
                  m_curvature(curvature),
                  m_pivot(pivot) {};

    	// Method
    	double eval(const double& x) {
            double xrel = x/m_pivot.MeV();
            return m_norm*std::pow(xrel, m_index + m_curvature * std::log(xrel));
        }

    	// Members
    protected:
    	const double&  m_norm;       //!< Normalization
    	const double&  m_index;	     //!< Spectral index at pivot
    	const double&  m_curvature;  //!< Curvature
    	const GEnergy& m_pivot;	     //!< Pivot energy
    };

    // Class to determine the integral energyflux, derviation of flux_kern
    class eflux_kern : public flux_kern {
    public:
    	// Constructor
    	eflux_kern(const double&  norm,
                   const double&  index,
                   const double&  curvature,
                   const GEnergy& pivot) :
                   flux_kern(norm, index, curvature, pivot) {};

    	// Method
    	double eval(const double& x) {
            return x * flux_kern::eval(x);
        }
    };


    // Protected members
    GModelPar m_norm;         //!< Normalization factor
    GModelPar m_index;        //!< Spectral index
    GModelPar m_curvature;    //!< Curvature
    GModelPar m_pivot;        //!< Pivot energy

    // Cached members used for pre-computations
    mutable GEnergy m_last_energy;     //!< Last energy value
    mutable double  m_last_index;      //!< Last index parameter
    mutable double  m_last_curvature;  //!< Last curvature parameters
    mutable double  m_last_pivot;      //!< Last pivot parameter
    mutable double  m_last_e_norm;     //!< Last E/Epivot value
    mutable double  m_last_log_e_norm; //!< Last ln(E/Epivot) value
    mutable double  m_last_exponent;   //!< Last exponent
    mutable double  m_last_power;      //!< Last power value

    // Cached members used for pre-computations
    mutable double m_mc_emin;       //!< Minimum energy
    mutable double m_mc_emax;       //!< Maximum energy
    mutable double m_mc_exponent;   //!< Exponent (index+1)
    mutable double m_mc_pow_emin;   //!< Power of minimum energy
    mutable double m_mc_pow_ewidth; //!< Power of energy width
    mutable double m_mc_norm; 	    //!< Norm of powerlaw model at logparabola pivot energy
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralLogParabola").
 ***************************************************************************/
inline
std::string GModelSpectralLogParabola::classname(void) const
{
    return ("GModelSpectralLogParabola");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "LogParabola".
 *
 * Returns the type of the log parabola spectral model.
 ***************************************************************************/
inline
std::string GModelSpectralLogParabola::type(void) const
{
    return "LogParabola";
}


/***********************************************************************//**
 * @brief Return pre factor
 *
 * @return Pre factor (ph/cm2/s/MeV).
 *
 * Returns the pre factor.
 ***************************************************************************/
inline
double GModelSpectralLogParabola::prefactor(void) const
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
void GModelSpectralLogParabola::prefactor(const double& prefactor)
{
    m_norm.value(prefactor);
    return;
}


/***********************************************************************//**
 * @brief Return spectral index
 *
 * @return Spectral index.
 *
 * Returns the spectral index.
 ***************************************************************************/
inline
double GModelSpectralLogParabola::index(void) const
{
    return (m_index.value());
}


/***********************************************************************//**
 * @brief Set spectral index 
 *
 * @param[in] index Spectral index.
 *
 * Sets the spectral index.
 ***************************************************************************/
inline
void GModelSpectralLogParabola::index(const double& index)
{
    m_index.value(index);
    return;
}


/***********************************************************************//**
 * @brief Return spectral curvature
 *
 * @return Spectral curvature.
 *
 * Returns the spectral curvature.
 ***************************************************************************/
inline
double GModelSpectralLogParabola::curvature(void) const
{
    return (m_curvature.value());
}


/***********************************************************************//**
 * @brief Set spectral curvature 
 *
 * @param[in] curvature Spectral curvature.
 *
 * Sets the spectral curvature.
 ***************************************************************************/
inline
void GModelSpectralLogParabola::curvature(const double& curvature)
{
    m_curvature.value(curvature);
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
GEnergy GModelSpectralLogParabola::pivot(void) const
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
void GModelSpectralLogParabola::pivot(const GEnergy& pivot)
{
    m_pivot.value(pivot.MeV());
    return;
}

#endif /* GMODELSPECTRALLOGPARABOLA_HPP */
