/***************************************************************************
 *         GModelSpectralNodes.hpp  -  Spectral nodes model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GModelSpectralNodes.hpp
 * @brief Spectral nodes model class definition
 * @author J. Knoedlseder
 */

#ifndef GMODELSPECTRALNODES_HPP
#define GMODELSPECTRALNODES_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GModelSpectralNodes
 *
 * @brief Spectral nodes model class
 *
 * This class implements a piecewise power law between spectral nodes
 * \f$(E_i, I_i)\f$,
 * where 
 * \f$E_i\f$ is the energy and \f$I_i\f$ is the intensity of node \f$i\f$.
 * For a given energy \f$E\f$, the piecewise powerlaw is computing by
 * finding the bracketing energies \f$E_1 <= E <= E_2\f$ and computing
 * \f[I(E)=10^{(\log v_1 + \log s_1) w_1 + (\log v_2 + \log s_2) w_2}\f]
 * where
 * \f$I_1 = v_1 s_1\f$ is the intensity of node 1 (\f$E_1\f$),
 * \f$I_2 = v_2 s_2\f$ is the intensity of node 2 (\f$E_2\f$),
 * \f$w_1\f$ is the weighting of node 1, and
 * \f$w_2\f$ is the weighting of node 2.
 * The weightings \f$w_1\f$ and \f$w_2\f$ are computed by linear
 * interpolation (in the log-log plane) between the nodes
 * \f$(\log E_1, \log I_1)\f$
 * and
 * \f$(\log E_2, \log I_2)\f$
 * to the requested energy \f$\log E\f$.
 ***************************************************************************/
class GModelSpectralNodes : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralNodes(void);
    explicit GModelSpectralNodes(const GXmlElement& xml);
    GModelSpectralNodes(const GModelSpectralNodes& model);
    virtual ~GModelSpectralNodes(void);

    // Operators
    virtual GModelSpectralNodes& operator=(const GModelSpectralNodes& model);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralNodes* clone(void) const;
    virtual std::string          type(void) const { return "NodeFunction"; }
    virtual double               eval(const GEnergy& srcEng) const;
    virtual double               eval_gradients(const GEnergy& srcEng) const;
    virtual double               flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralNodes& model);
    void free_members(void);
    void update_pars(void);
    void load_nodes(const std::string& filename);
    void set_cache(void) const;
    void set_eval_cache(void) const;
    void set_flux_cache(void) const;
    void update_eval_cache(void) const;
    void update_flux_cache(void) const;
    void mc_update(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    std::vector<GModelPar>      m_energies;     //!< Node energies
    std::vector<GModelPar>      m_values;       //!< Node values

    // Evaluation cache
    mutable std::vector<double> m_old_energies; //!< Old energies
    mutable std::vector<double> m_old_values;   //!< Old values
    mutable GNodeArray          m_log_energies; //!< log10(energy) of nodes
    mutable std::vector<double> m_log_values;   //!< log10(value) of nodes

    // Flux computation cache
    mutable GNodeArray          m_lin_energies; //!< Energy of nodes
    mutable std::vector<double> m_lin_values;   //!< Values of nodes
    mutable std::vector<double> m_prefactor;    //!< Power-law normalisations
    mutable std::vector<double> m_gamma;        //!< Power-law indices
    mutable std::vector<double> m_epivot;       //!< Power-law pivot energies
    mutable std::vector<double> m_flux;         //!< Photon fluxes
    mutable std::vector<double> m_eflux;        //!< Energy fluxes
    
    // Cached members for MC
    mutable GEnergy             m_mc_emin;      //!< Minimum energy
    mutable GEnergy             m_mc_emax;      //!< Maximum energy
    mutable std::vector<double> m_mc_cum;       //!< Cumulative distribution
    mutable std::vector<double> m_mc_min;       //!< Lower boundary for MC
    mutable std::vector<double> m_mc_max;       //!< Upper boundary for MC
    mutable std::vector<double> m_mc_exp;       //!< Exponent for MC
};

#endif /* GMODELSPECTRALNODES_HPP */
