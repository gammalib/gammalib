/***************************************************************************
 *   GModelSpectralPlawPhotonFlux.hpp - Spectral power law model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlawPhotonFlux.hpp
 * @brief Flux normalized power law spectral model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALPLAWPHOTONFLUX_HPP
#define GMODELSPECTRALPLAWPHOTONFLUX_HPP

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
 * @class GModelSpectralPlawPhotonFlux
 *
 * @brief Photon flux normalized power law spectral model class
 *
 * This class implements a power law spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_flux}
 *    \frac{{\tt m\_index}+1}
 *         {{\tt e\_max}^{{\tt m\_index}+1} -
 *          {\tt e\_min}^{{\tt m\_index}+1}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} \ne -1\f$ and
 *
 * \f[
 *    S_{\rm E}(E | t) = 
 *    \frac{{\tt m\_flux}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_flux}\f$ is the photon flux between
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 ***************************************************************************/
class GModelSpectralPlawPhotonFlux : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlawPhotonFlux(void);
    GModelSpectralPlawPhotonFlux(const std::string& type,
                                 const std::string& flux,
                                 const std::string& index,
                                 const std::string& emin,
                                 const std::string& emax);
    GModelSpectralPlawPhotonFlux(const double&  flux,
                                 const double&  index,
                                 const GEnergy& emin,
                                 const GEnergy& emax);
    explicit GModelSpectralPlawPhotonFlux(const GXmlElement& xml);
    GModelSpectralPlawPhotonFlux(const GModelSpectralPlawPhotonFlux& model);
    virtual ~GModelSpectralPlawPhotonFlux(void);

    // Operators
    virtual GModelSpectralPlawPhotonFlux& operator=(const GModelSpectralPlawPhotonFlux& model);

    // Implemented pure virtual base class methods
    virtual void                          clear(void);
    virtual GModelSpectralPlawPhotonFlux* clone(void) const;
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
    void    type(const std::string& type);
    double  flux(void) const;
    void    flux(const double& flux);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(void) const;
    void    emin(const GEnergy& emin);
    GEnergy emax(void) const;
    void    emax(const GEnergy& emax);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralPlawPhotonFlux& model);
    void free_members(void);
    void update(const GEnergy& srcEng) const;

    // Protected members
    std::string     m_type;            //!< Model type
    GModelPar       m_flux;            //!< Photon flux (ph/cm2/s)
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
 * @return String containing the class name ("GModelSpectralPlawPhotonFlux").
 ***************************************************************************/
inline
std::string GModelSpectralPlawPhotonFlux::classname(void) const
{
    return ("GModelSpectralPlawPhotonFlux");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spectral power law model.
 ***************************************************************************/
inline
std::string GModelSpectralPlawPhotonFlux::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Set model type
 *
 * @param[in] type Model type.
 *
 * Set the type of the spectral power law model.
 ***************************************************************************/
inline
void GModelSpectralPlawPhotonFlux::type(const std::string& type)
{
    m_type = type;
    return;
}


/***********************************************************************//**
 * @brief Return photon flux
 *
 * @return Photon flux (ph/cm2/s).
 *
 * Returns the photon flux.
 ***************************************************************************/
inline
double GModelSpectralPlawPhotonFlux::flux(void) const
{
    return (m_flux.value());
}


/***********************************************************************//**
 * @brief Set photon flux
 *
 * @param[in] flux Photon flux (ph/cm2/s).
 *
 * Sets the photon flux.
 ***************************************************************************/
inline
void GModelSpectralPlawPhotonFlux::flux(const double& flux)
{
    m_flux.value(flux);
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
double GModelSpectralPlawPhotonFlux::index(void) const
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
void GModelSpectralPlawPhotonFlux::index(const double& index)
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
GEnergy GModelSpectralPlawPhotonFlux::emin(void) const
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
void GModelSpectralPlawPhotonFlux::emin(const GEnergy& emin)
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
GEnergy GModelSpectralPlawPhotonFlux::emax(void) const
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
void GModelSpectralPlawPhotonFlux::emax(const GEnergy& emax)
{
    m_emax.value(emax.MeV());
    return;
}

#endif /* GMODELSPECTRALPLAWPHOTONFLUX_HPP */
