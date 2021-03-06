/***************************************************************************
 *         GModelSpectralFunc.hpp - Spectral function model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALFUNC_HPP
#define GMODELSPECTRALFUNC_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"
#include "GNodeArray.hpp"
#include "GFilename.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GEnergies;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralFunc
 *
 * @brief Spectral function model class
 *
 * This class implements an arbitrary function  as the spectral model
 * component. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization of the function.
 ***************************************************************************/
class GModelSpectralFunc : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralFunc(void);
    GModelSpectralFunc(const GFilename& filename, const double& norm);
    GModelSpectralFunc(const GModelSpectral& model, const GEnergies& energies);
    explicit GModelSpectralFunc(const GXmlElement& xml);
    GModelSpectralFunc(const GModelSpectralFunc& model);
    virtual ~GModelSpectralFunc(void);

    // Operators
    virtual GModelSpectralFunc& operator=(const GModelSpectralFunc& model);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GModelSpectralFunc* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime&   srcTime = GTime(),
                                     const bool&    gradients = false) const;
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
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int              nodes(void) const;
    bool             is_empty(void) const;
    void             append(const GEnergy& energy, const double& intensity);
    void             insert(const GEnergy& energy, const double& intensity);
    void             remove(const int& index);
    void             reserve(const int& num);
    void             extend(const GModelSpectralFunc& filefct);
    GEnergy          energy(const int& index) const;
    void             energy(const int& index, const GEnergy& energy);
    double           intensity(const int& index) const;
    void             intensity(const int& index, const double& intensity);
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
    double           norm(void) const;
    void             norm(const double& norm);
    void             save(const GFilename& filename,
                          const bool&      clobber = false) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralFunc& model);
    void free_members(void);
    void load_nodes(const GFilename& filename);
    void set_cache(void) const;
    void mc_update(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    GModelPar           m_norm;       //!< Normalization factor
    mutable GNodeArray  m_lin_nodes;  //!< Energy nodes of function
    mutable GNodeArray  m_log_nodes;  //!< log10(Energy) nodes of function
    std::vector<double> m_lin_values; //!< Function values at nodes
    std::vector<double> m_log_values; //!< log10(Function) values at nodes
    mutable GFilename   m_filename;   //!< Name of file function

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


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralFunc").
 ***************************************************************************/
inline
std::string GModelSpectralFunc::classname(void) const
{
    return ("GModelSpectralFunc");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "FileFunction".
 *
 * Returns the type of the spectral function model.
 ***************************************************************************/
inline
std::string GModelSpectralFunc::type(void) const
{
    return "FileFunction";
}


/***********************************************************************//**
 * @brief Return number of nodes in file function
 *
 * @return Number of nodes in file function.
 *
 * Returns the number of nodes in the file function model.
 ***************************************************************************/
inline
int GModelSpectralFunc::nodes(void) const
{
    return (int)m_lin_nodes.size();
}


/***********************************************************************//**
 * @brief Signals if there are nodes in file function
 *
 * @return True if file function is empty, false otherwise.
 *
 * Signals if the file function does not contain any node.
 ***************************************************************************/
inline
bool GModelSpectralFunc::is_empty(void) const
{
    return (m_lin_nodes.is_empty());
}


/***********************************************************************//**
 * @brief Reserves space for nodes in file function
 *
 * @param[in] num Number of nodes
 *
 * Reserves space for @p num nodes in file function.
 ***************************************************************************/
inline
void GModelSpectralFunc::reserve(const int& num)
{
    m_lin_nodes.reserve(num);
    m_log_nodes.reserve(num);
    m_lin_values.reserve(num);
    m_log_values.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor.
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelSpectralFunc::norm(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set normalization factor 
 *
 * @param[in] norm Normalization factor.
 *
 * Sets the normalization factor.
 ***************************************************************************/
inline
void GModelSpectralFunc::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return Name of node file.
 *
 * Returns the name of the file function node file.
 ***************************************************************************/
inline
const GFilename& GModelSpectralFunc::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Loads nodes from node file and set filename
 *
 * @param[in] filename Node file name.
 *
 * Loads the nodes from a file function node file and sets the filename.
 ***************************************************************************/
inline
void GModelSpectralFunc::filename(const GFilename& filename)
{
    load_nodes(filename);
    return;
}

#endif /* GMODELSPECTRALFUNC_HPP */
