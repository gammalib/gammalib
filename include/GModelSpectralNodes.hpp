/***************************************************************************
 *          GModelSpectralNodes.hpp - Spectral nodes model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALNODES_HPP
#define GMODELSPECTRALNODES_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GNodeArray.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GEnergies;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralNodes
 *
 * @brief Spectral nodes model class
 *
 * This class implements a piecewise power law between spectral nodes
 *
 * \f[
 *    ({\tt m\_energies[i]}, {\tt m\_values[i]})
 * \f]
 *
 * where 
 * - \f${\tt m\_energies[i]}\f$ is the energy, and
 * - \f${\tt m\_values[i]}\f$ is the intensity (in photons/cm2/s/MeV)
 *   of node \f$i\f$.
 *
 * For a given energy \f$E\f$, the piecewise power law is computing by
 * finding the bracketing energies 
 * \f${\tt m\_energies[i]} <= E <= {\tt m\_energies[i+1]}\f$ and computing
 *
 * \f[
 *    S_{\rm E}(E | t) =
 *    10^{(\log {\tt m\_values[i]}) w_{i} + 
 *        (\log {\tt m\_values[i+1]}) w_{i+1}}
 * \f]
 *
 * where
 * - \f${\tt m\_values[i]}\f$ is the intensity of node \f$i\f$,
 * - \f${\tt m\_values[i+1]}\f$ is the intensity of node \f$i+1\f$,
 * - \f$w_{i}\f$ is the weighting of node \f$i\f$, and
 * - \f$w_{i+1}\f$ is the weighting of node \f$i+1\f$.
 *
 * The weightings \f$w_{i}\f$ and \f$w_{i+1}\f$ are computed by linear
 * interpolation (in the log-log plane) between the nodes
 * \f$(\log {\tt m\_energies[i]}, \log{\tt m\_values[i]})\f$
 * and
 * \f$(\log {\tt m\_energies[i+1]}, \log{\tt m\_values[i+1]})\f$
 * to the requested energy \f$\log E\f$.
 ***************************************************************************/
class GModelSpectralNodes : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralNodes(void);
    GModelSpectralNodes(const GModelSpectral& model, const GEnergies& energies);
    explicit GModelSpectralNodes(const GXmlElement& xml);
    GModelSpectralNodes(const GModelSpectralNodes& model);
    virtual ~GModelSpectralNodes(void);

    // Operators
    virtual GModelSpectralNodes& operator=(const GModelSpectralNodes& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralNodes* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime& srcTime = GTime(),
                                      const bool& gradients = false) const;
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
    int     nodes(void) const;
    void    append(const GEnergy& energy, const double& intensity);
    void    insert(const int& index, const GEnergy& energy,
                   const double& intensity);
    void    remove(const int& index);
    void    reserve(const int& num);
    void    extend(const GModelSpectralNodes& nodes);
    GEnergy energy(const int& index) const;
    void    energy(const int& index, const GEnergy& energy);
    double  intensity(const int& index) const;
    void    intensity(const int& index, const double& intensity);
    double  error(const int& index) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralNodes& model);
    void free_members(void);
    void update_pars(void);
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


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralNodes").
 ***************************************************************************/
inline
std::string GModelSpectralNodes::classname(void) const
{
    return ("GModelSpectralNodes");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "NodeFunction".
 *
 * Returns the type of the spectral node function model.
 ***************************************************************************/
inline
std::string GModelSpectralNodes::type(void) const
{
    return "NodeFunction";
}


/***********************************************************************//**
 * @brief Return number of nodes
 *
 * @return Number of nodes.
 *
 * Returns the number of nodes in the node function model.
 ***************************************************************************/
inline
int GModelSpectralNodes::nodes(void) const
{
    return (int)m_energies.size();
}

#endif /* GMODELSPECTRALNODES_HPP */
