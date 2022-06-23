/***************************************************************************
 *          GModelSpectralTable.hpp - Spectral table model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTable.hpp
 * @brief Spectral table model class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALTABLE_HPP
#define GMODELSPECTRALTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelSpectral.hpp"
#include "GEbounds.hpp"
#include "GModelPar.hpp"
#include "GNdarray.hpp"
#include "GNodeArray.hpp"
#include "GFilename.hpp"
#include "GModelSpectralTablePars.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GEnergy;
class GTime;
class GFits;
class GXmlElement;
class GFitsBinTable;


/***********************************************************************//**
 * @class GModelSpectralTable
 *
 * @brief Spectral table model class
 ***************************************************************************/
class GModelSpectralTable : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralTable(void);
    GModelSpectralTable(const GFilename& filename, const double& norm);
    GModelSpectralTable(const GEbounds&                ebounds,
                        const GModelSpectralTablePars& pars,
                        const GNdarray&                spectra);
    explicit GModelSpectralTable(const GXmlElement& xml);
    GModelSpectralTable(const GModelSpectralTable& model);
    virtual ~GModelSpectralTable(void);

    // Operators
    virtual GModelSpectralTable& operator=(const GModelSpectralTable& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralTable* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime = GTime(),
                                      const bool&    gradients = false) const;
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
    int                           nspectra(void) const;
    GModelSpectralTablePar&       table_par(const int& index);
    const GModelSpectralTablePar& table_par(const int& index) const;
    GModelSpectralTablePar&       table_par(const std::string& name);
    const GModelSpectralTablePar& table_par(const std::string& name) const;
    double                        norm(void) const;
    void                          norm(const double& norm);
    const GEbounds&               ebounds(void) const;
    void                          load(const GFilename& filename);
    void                          save(const GFilename& filename,
                                       const bool&      clobber = false) const;
    const GFilename&              filename(void) const;
    void                          energy_scale(const std::string& name);

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModelSpectralTable& model);
    void          free_members(void);
    void          set_par_pointers(void);
    void          set_energy_nodes(void);
    void          scale_energy(void);
    bool          has_energy_scale(void) const;
    GFitsBinTable create_par_table(void) const;
    GFitsBinTable create_eng_table(void) const;
    GFitsBinTable create_spec_table(void) const;
    void          load_par(const GFits& fits);
    void          load_eng(const GFits& fits);
    void          load_spec(const GFits& fits);
    int           par_index(const std::string& name) const;
    void          update(void) const;
    void          update_flux(void) const;
    void          update_mc(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    GModelPar                   m_norm;        //!< Normalization factor
    GModelSpectralTablePars     m_table_pars;  //!< Table model parameters
    GNdarray                    m_spectra;     //!< Spectra
    GEbounds                    m_ebounds;     //!< Energy boundaries
    mutable GFilename           m_filename;    //!< Filename of table
    std::string                 m_escale_par;  //!< Energy scaling parameter
    double                      m_escale;      //!< Energy scale
    double                      m_log10escale; //!< log10 of energy scale

    // Cached members used for pre-computations
    mutable int                 m_npars;       //!< Number of parameters
    mutable int                 m_nebins;      //!< Number of energy bins
    mutable std::vector<double> m_last_values; //!< Last parameter values
    mutable GNodeArray          m_lin_nodes;   //!< Energy nodes of function
    mutable GNodeArray          m_log_nodes;   //!< log10(Energy) nodes of function
    mutable GNdarray            m_lin_values;  //!< Function values and grad's
    mutable GNdarray            m_log_values;  //!< log10(Function) values and grad's
    mutable std::vector<double> m_prefactor;   //!< Power-law normalisations
    mutable std::vector<double> m_gamma;       //!< Power-law indices
    mutable std::vector<double> m_epivot;      //!< Power-law pivot energies
    mutable std::vector<double> m_flux;        //!< Photon fluxes
    mutable std::vector<double> m_eflux;       //!< Energy fluxes
    mutable GEnergy             m_mc_emin;     //!< Minimum energy
    mutable GEnergy             m_mc_emax;     //!< Maximum energy
    mutable std::vector<double> m_mc_cum;      //!< Cumulative distribution
    mutable std::vector<double> m_mc_min;      //!< Lower boundary for MC
    mutable std::vector<double> m_mc_max;      //!< Upper boundary for MC
    mutable std::vector<double> m_mc_exp;      //!< Exponent for MC
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralTable").
 ***************************************************************************/
inline
std::string GModelSpectralTable::classname(void) const
{
    return ("GModelSpectralTable");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "TableModel".
 *
 * Returns the type of the spectral model.
 ***************************************************************************/
inline
std::string GModelSpectralTable::type(void) const
{
    return ("TableModel");
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor.
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelSpectralTable::norm(void) const
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
void GModelSpectralTable::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}


/***********************************************************************//**
 * @brief Return reference to energy boundaries
 *
 * @return Reference to energy boundaries.
 *
 * Returns a reference to energy boundaries.
 ***************************************************************************/
inline
const GEbounds& GModelSpectralTable::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return Name of table file.
 *
 * Returns the name of the table file.
 ***************************************************************************/
inline
const GFilename& GModelSpectralTable::filename(void) const
{
    return (m_filename);
}

#endif /* GMODELSPECTRALTABLE_HPP */
