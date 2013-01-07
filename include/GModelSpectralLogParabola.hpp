/***************************************************************************
 *    GModelSpectralLogParabola.hpp - Log parabola spectral model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Michael Mayer                                    *
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
 * @brief LogParabola spectrall model class
 *
 * This class implements a logparabole as the spectral component of the
 * gamma-ray sky model. The power law is defined as
 * \f[I(E)=norm (E/pivot)^{index-curvature\ln(E/pivot)}\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy, and
 * \f$index\f$ is the spectral index.
 * \f$curvature\f$ is the curvature
 ***************************************************************************/
class GModelSpectralLogParabola : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralLogParabola(void);
    explicit GModelSpectralLogParabola(const double& norm, const double& index, const double &curvature);
    explicit GModelSpectralLogParabola(const GXmlElement& xml);
    GModelSpectralLogParabola(const GModelSpectralLogParabola& model);
    virtual ~GModelSpectralLogParabola(void);

    // Operators
    virtual GModelSpectralLogParabola& operator=(const GModelSpectralLogParabola& model);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GModelSpectralLogParabola* clone(void) const;
    virtual std::string                type(void) const { return "LogParabola"; }
    virtual double                     eval(const GEnergy& srcEng) const;
    virtual double                     eval_gradients(const GEnergy& srcEng) const;
    virtual double                     flux(const GEnergy& emin,
                                            const GEnergy& emax) const;
    virtual double                     eflux(const GEnergy& emin,
                                             const GEnergy& emax) const;
    virtual GEnergy                    mc(const GEnergy& emin,
                                          const GEnergy& emax, GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;
    virtual std::string                print(void) const;

    // Other methods
    void   autoscale(void);
    double norm(void) const { return m_norm.real_value(); }
    double index(void) const { return m_index.real_value(); }
    double curvature(void) const { return m_curvature.real_value(); }
    double pivot(void) const { return m_pivot.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralLogParabola& model);
    void free_members(void);
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Class to determine to the integral photon flux
    class flux_kern : public GFunction {
    public:
    	// Constructor
    	flux_kern(double norm, double index, double curvature, double pivot) :
    		m_norm(norm), m_index(index), m_curvature(curvature), m_pivot(pivot) {};

    	// Method
    	double eval(double x) {
            return m_norm*std::pow(x/m_pivot,m_index+m_curvature*std::log(x/m_pivot));
        }

    	// Members
    protected:
    	double m_norm;        //!< Normalization
    	double m_index;	      //!< Spectral index at pivot
    	double m_curvature;   //!< Curvature
    	double m_pivot;	      //!< Pivot energy
    };

    // Class to determine the integral energyflux, derviation of flux_kern
    class eflux_kern : public flux_kern {
    public:
    	//Constructor
    	eflux_kern(double norm, double index, double curvature, double pivot) :
            flux_kern(norm,index,curvature,pivot) {};

    	// Method
    	double eval(double x) {
            return x*flux_kern::eval(x);
        }
    };


    // Protected members
    GModelPar m_norm;         //!< Normalization factor
    GModelPar m_index;        //!< Spectral index
    GModelPar m_curvature;    //!< Curvature
    GModelPar m_pivot;        //!< Pivot energy

    // Cached members used for pre-computations
    mutable double m_mc_emin;       //!< Minimum energy
    mutable double m_mc_emax;       //!< Maximum energy
    mutable double m_mc_exponent;   //!< Exponent (index+1)
    mutable double m_mc_pow_emin;   //!< Power of minimum energy
    mutable double m_mc_pow_ewidth; //!< Power of energy width
};

#endif /* GMODELSPECTRALLOGPARABOLA_HPP */
