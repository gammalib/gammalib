/***************************************************************************
 *           GModelSpectralBins.hpp - Spectral bins model class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GModelSpectralBins.hpp
 * @brief Spectral bins model class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALBINS_HPP
#define GMODELSPECTRALBINS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GEbounds;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralBins
 *
 * @brief Spectral bins model class
 *
 * This class implements spectral bins that have a constant intensity within
 * their boundaries. A spectral bin is defined by a lower and upper energy
 * limit, with the lower limit included and the upper limit excluded from an
 * energy bin.
 ***************************************************************************/
class GModelSpectralBins : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralBins(void);
    GModelSpectralBins(const GModelSpectral& model,
                       const GEbounds&       ebounds,
                       const double&         index = -2.0);
    explicit GModelSpectralBins(const GXmlElement& xml);
    GModelSpectralBins(const GModelSpectralBins& model);
    virtual ~GModelSpectralBins(void);

    // Operators
    virtual GModelSpectralBins& operator=(const GModelSpectralBins& model);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GModelSpectralBins* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime& srcTime = GTime(),
                                     const bool& gradients = false) const;
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
    int     bins(void) const;
    void    append(const GEnergy& emin,
                   const GEnergy& emax,
                   const double&  intensity);
    void    insert(const int&     index,
                   const GEnergy& emin,
                   const GEnergy& emax,
                   const double&  intensity);
    void    remove(const int& index);
    void    reserve(const int& num);
    void    extend(const GModelSpectralBins& bins);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(const int& index) const;
    GEnergy emax(const int& index) const;
    void    emin(const int& index, const GEnergy& emin);
    void    emax(const int& index, const GEnergy& emax);
    double  intensity(const int& index) const;
    void    intensity(const int& index, const double& intensity);
    double  error(const int& index) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralBins& model);
    void free_members(void);
    int  bin_index(const GEnergy& energy) const;
    void update_pars(void);
    void set_cache(void) const;
    void mc_update(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    GModelPar                   m_index;     //!< Spectral index of all bins
    std::vector<GModelPar>      m_emin;      //!< Lower energy limits
    std::vector<GModelPar>      m_emax;      //!< Upper energy limits
    std::vector<GModelPar>      m_values;    //!< Bin values

    // Evaluation cache
    mutable std::vector<double> m_epivot;    //!< Power-law pivot energies in MeV

    // MC cache
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
 * @return String containing the class name ("GModelSpectralBins").
 ***************************************************************************/
inline
std::string GModelSpectralBins::classname(void) const
{
    return ("GModelSpectralBins");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "BinFunction".
 *
 * Returns the type of the spectral bin function model.
 ***************************************************************************/
inline
std::string GModelSpectralBins::type(void) const
{
    return "BinFunction";
}


/***********************************************************************//**
 * @brief Return number of bins
 *
 * @return Number of bins.
 *
 * Returns the number of bins in the bin function model.
 ***************************************************************************/
inline
int GModelSpectralBins::bins(void) const
{
    return (int)m_values.size();
}

#endif /* GMODELSPECTRALBINS_HPP */
