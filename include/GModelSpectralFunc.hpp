/***************************************************************************
 *         GModelSpectralFunc.hpp - Spectral function model class          *
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
 * @file GModelSpectralFunc.hpp
 * @brief Spectral function model class definition
 * @author J. Knoedlseder
 */

#ifndef GMODELSPECTRALFUNC_HPP
#define GMODELSPECTRALFUNC_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GModelSpectralFunc
 *
 * @brief Spectral function model class
 *
 * This class implements an arbitrary function as the spectral component of
 * the model. The function is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 ***************************************************************************/
class GModelSpectralFunc : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralFunc(void);
    explicit GModelSpectralFunc(const std::string& filename);
    explicit GModelSpectralFunc(const GXmlElement& xml);
    GModelSpectralFunc(const GModelSpectralFunc& model);
    virtual ~GModelSpectralFunc(void);

    // Operators
    virtual GModelSpectralFunc& operator=(const GModelSpectralFunc& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralFunc* clone(void) const;
    virtual std::string         type(void) const { return "FileFunction"; }
    virtual double              eval(const GEnergy& srcEng) const;
    virtual double              eval_gradients(const GEnergy& srcEng) const;
    virtual double              flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double norm(void) const { return m_norm.Value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralFunc& model);
    void free_members(void);
    void load_nodes(const std::string& filename);
    void set_cache(void) const;
    void mc_update(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    GModelPar           m_norm;       //!< Normalization factor
    mutable GNodeArray  m_lin_nodes;  //!< Energy nodes of function
    mutable GNodeArray  m_log_nodes;  //!< lof10(Energy) nodes of function
    std::vector<double> m_lin_values; //!< Function values at nodes
    std::vector<double> m_log_values; //!< log10(Function) values at nodes
    std::string         m_filename;   //!< Name of file function

    // Cached members used for pre-computations
    mutable std::vector<double> m_prefactor; //!< Power-law normalisations
    mutable std::vector<double> m_gamma;     //!< Power-law indices
    mutable std::vector<double> m_epivot;    //!< Power-law pivot energies
    mutable std::vector<double> m_flux;      //!< Photon fluxes
    mutable std::vector<double> m_eflux;     //!< Energy fluxes
    
    // Cached members for MC
    mutable GEnergy             m_mc_emin;   //!< Minimum energy
    mutable GEnergy             m_mc_emax;   //!< Maximum energy
    mutable std::vector<double> m_mc_cum;    //!< Cumulative distribution
    mutable std::vector<double> m_mc_min;    //!< Lower boundary for MC
    mutable std::vector<double> m_mc_max;    //!< Upper boundary for MC
    mutable std::vector<double> m_mc_exp;    //!< Exponent for MC
};

#endif /* GMODELSPECTRALFUNC_HPP */
