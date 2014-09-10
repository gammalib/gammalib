/***************************************************************************
 *        GModelSpectralPlaw2.hpp - Spectral power law model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlaw2.hpp
 * @brief Flux normalized power law spectral model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALPLAW2_HPP
#define GMODELSPECTRALPLAW2_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralPlaw2
 *
 * @brief Flux normalized power law spectral model class
 *
 * This class implements a power law spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_integral}
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
 *    \frac{{\tt m\_integral}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_integral}\f$ is the integral flux between 
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 ***************************************************************************/
class GModelSpectralPlaw2 : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlaw2(void);
    explicit GModelSpectralPlaw2(const double&  integral,
                                 const double&  index,
                                 const GEnergy& emin,
                                 const GEnergy& emax);
    explicit GModelSpectralPlaw2(const GXmlElement& xml);
    GModelSpectralPlaw2(const GModelSpectralPlaw2& model);
    virtual ~GModelSpectralPlaw2(void);

    // Operators
    virtual GModelSpectralPlaw2& operator=(const GModelSpectralPlaw2& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralPlaw2* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime&   srcTime);
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
    double  integral(void) const;
    void    integral(const double& integral);
    double  index(void) const;
    void    index(const double& index);
    GEnergy emin(void) const;
    void    emin(const GEnergy& emin);
    GEnergy emax(void) const;
    void    emax(const GEnergy& emax);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralPlaw2& model);
    void free_members(void);
    void update(const GEnergy& srcEng) const;

    // Protected members
    GModelPar       m_integral;        //!< Integral flux
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
    mutable double  m_last_integral;   //!< Last integral flux
    mutable double  m_last_index;      //!< Last spectral index (MeV)
    mutable GEnergy m_last_emin;       //!< Last lower energy limit
    mutable GEnergy m_last_emax;       //!< Last upper energy limit
    mutable GEnergy m_last_energy;     //!< Last source energy
    mutable double  m_last_value;      //!< Last function value
    mutable double  m_last_g_integral; //!< Last integral flux gradient
    mutable double  m_last_g_index;    //!< Last spectral index gradient
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralPlaw2").
 ***************************************************************************/
inline
std::string GModelSpectralPlaw2::classname(void) const
{
    return ("GModelSpectralPlaw2");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "PowerLaw2".
 *
 * Returns the type of the spectral power law model.
 ***************************************************************************/
inline
std::string GModelSpectralPlaw2::type(void) const
{
    return "PowerLaw2";
}


/***********************************************************************//**
 * @brief Return integral flux
 *
 * @return Integral flux (ph/cm2/s).
 *
 * Returns the integral flux.
 ***************************************************************************/
inline
double GModelSpectralPlaw2::integral(void) const
{
    return (m_integral.value());
}


/***********************************************************************//**
 * @brief Set integral flux
 *
 * @param[in] integral Integral flux (ph/cm2/s).
 *
 * Sets the integral flux.
 ***************************************************************************/
inline
void GModelSpectralPlaw2::integral(const double& integral)
{
    m_integral.value(integral);
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
double GModelSpectralPlaw2::index(void) const
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
void GModelSpectralPlaw2::index(const double& index)
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
GEnergy GModelSpectralPlaw2::emin(void) const
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
void GModelSpectralPlaw2::emin(const GEnergy& emin)
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
GEnergy GModelSpectralPlaw2::emax(void) const
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
void GModelSpectralPlaw2::emax(const GEnergy& emax)
{
    m_emax.value(emax.MeV());
    return;
}

#endif /* GMODELSPECTRALPLAW2_HPP */
