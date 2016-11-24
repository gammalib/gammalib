/***************************************************************************
 *    GModelSpectralPlawEnergyFlux.hpp - Spectral power law model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                                    *
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
 * @file GModelSpectralPlawEnergyFlux.hpp
 * @brief Energy flux normalized power law spectral model class interface definition
 * @author Michael Mayer
 */

#ifndef GMODELSPECTRALPLAWENERGYFLUX_HPP
#define GMODELSPECTRALPLAWENERGYFLUX_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralPlawEnergyFlux
 *
 * @brief Energy flux normalized power law spectral model class
 *
 * This class implements a power law spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_eflux}
 *    \frac{{\tt m\_index}+1}
 *         {{\tt e\_max}^{{\tt m\_index}+2} -
 *          {\tt e\_min}^{{\tt m\_index}+2}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} \ne -2\f$ and
 *
 * \f[
 *    S_{\rm E}(E | t) = 
 *    \frac{{\tt m\_eflux}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -2\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_eflux}\f$ is the energy flux between
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 ***************************************************************************/
class GModelSpectralPlawEnergyFlux : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlawEnergyFlux(void);
    GModelSpectralPlawEnergyFlux(const std::string& type,
                                 const std::string& eflux,
                                 const std::string& index,
                                 const std::string& emin,
                                 const std::string& emax);
    explicit GModelSpectralPlawEnergyFlux(const double&  eflux,
                                          const double&  index,
                                          const GEnergy& emin,
                                          const GEnergy& emax);
    explicit GModelSpectralPlawEnergyFlux(const GXmlElement& xml);
    GModelSpectralPlawEnergyFlux(const GModelSpectralPlawEnergyFlux& model);
    virtual ~GModelSpectralPlawEnergyFlux(void);

    // Operators
    virtual GModelSpectralPlawEnergyFlux& operator=(const GModelSpectralPlawEnergyFlux& model);

    // Implemented pure virtual base class methods
    virtual void                          clear(void);
    virtual GModelSpectralPlawEnergyFlux* clone(void) const;
    virtual std::string                   classname(void) const;
    virtual std::string                   type(void) const;
    virtual double                        eval(const GEnergy& srcEng,
                                               const GTime&   srcTime = GTime(),
                                               const bool&    gradients = false) const;
    virtual double                        flux(const GEnergy& emin,
                                               const GEnergy& emax) const;
    virtual double                        eflux(const GEnergy& emin,
                                                const GEnergy& emax) const;
    virtual GEnergy                       mc(const GEnergy& emin,
                                             const GEnergy& emax,
                                             const GTime&   time,
                                             GRan&          ran) const;
    virtual void                          read(const GXmlElement& xml);
    virtual void                          write(GXmlElement& xml) const;
    virtual std::string                   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  eflux(void) const;
    void    eflux(const double& eflux);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(void) const;
    void    emin(const GEnergy& emin);
    GEnergy emax(void) const;
    void    emax(const GEnergy& emax);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralPlawEnergyFlux& model);
    void free_members(void);
    void update(const GEnergy& srcEng) const;

    // Protected members
	std::string     m_type;            //!< Model type
	GModelPar       m_eflux;           //!< Energy flux (erg/cm2/s)
	GModelPar       m_index;           //!< Spectral index
	GModelPar       m_emin;            //!< Lower energy limit (MeV)
	GModelPar       m_emax;            //!< Upper energy limit (MeV)

	// Cached members used for pre-computations
	mutable double  m_log_emin;        //!< Log(emin)
	mutable double  m_log_emax;        //!< Log(emax)
	mutable double  m_pow_emin;        //!< emin^(index+1)
	mutable double  m_pow_emax;        //!< emax^(index+1)
	mutable double  m_norm;            //!< Power-law normalization (for pivot energy 1 MeV)
	mutable double  m_g_norm;          //!< Power-law normalization gradient
	mutable double  m_power;           //!< Power-law factor
	mutable double  m_last_index;      //!< Last spectral index (MeV)
	mutable GEnergy m_last_emin;       //!< Last lower energy limit
	mutable GEnergy m_last_emax;       //!< Last upper energy limit
	mutable GEnergy m_last_energy;     //!< Last source energy
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralPlawEnergyFlux").
 ***************************************************************************/
inline
std::string GModelSpectralPlawEnergyFlux::classname(void) const
{
    return ("GModelSpectralPlawEnergyFlux");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spectral power law model.
 ***************************************************************************/
inline
std::string GModelSpectralPlawEnergyFlux::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return energy flux
 *
 * @return Energy flux (erg/cm2/s).
 *
 * Returns the energy flux.
 ***************************************************************************/
inline
double GModelSpectralPlawEnergyFlux::eflux(void) const
{
    return (m_eflux.value());
}


/***********************************************************************//**
 * @brief Set energy flux
 *
 * @param[in] eflux Energy flux (erg/cm2/s).
 *
 * Sets the energy flux.
 ***************************************************************************/
inline
void GModelSpectralPlawEnergyFlux::eflux(const double& eflux)
{
    m_eflux.value(eflux);
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
double GModelSpectralPlawEnergyFlux::index(void) const
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
void GModelSpectralPlawEnergyFlux::index(const double& index)
{
    m_index.value(index);
    return;
}


/***********************************************************************//**
 * @brief Return minimum energy
 *
 * @return Minimum energy.
 *
 * Returns the minimum energy.
 ***************************************************************************/
inline
GEnergy GModelSpectralPlawEnergyFlux::emin(void) const
{
    GEnergy energy;
    energy.MeV(m_emin.value());
    return energy;
}


/***********************************************************************//**
 * @brief Set minimum energy
 *
 * @param[in] emin Minimum energy.
 *
 * Sets the minimum energy.
 ***************************************************************************/
inline
void GModelSpectralPlawEnergyFlux::emin(const GEnergy& emin)
{
    m_emin.value(emin.MeV());
    return;
}


/***********************************************************************//**
 * @brief Return maximum energy
 *
 * @return Maximum energy.
 *
 * Returns the maximum energy.
 ***************************************************************************/
inline
GEnergy GModelSpectralPlawEnergyFlux::emax(void) const
{
    GEnergy energy;
    energy.MeV(m_emax.value());
    return energy;
}


/***********************************************************************//**
 * @brief Set maximum energy
 *
 * @param[in] emax Maximum energy.
 *
 * Sets the maximum energy.
 ***************************************************************************/
inline
void GModelSpectralPlawEnergyFlux::emax(const GEnergy& emax)
{
    m_emax.value(emax.MeV());
    return;
}

#endif /* GMODELSPECTRALPLAWENERGYFLUX_HPP */
