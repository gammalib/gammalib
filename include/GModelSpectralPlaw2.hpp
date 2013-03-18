/***************************************************************************
 *        GModelSpectralPlaw2.hpp - Spectral power law model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * This class implements a power law as the spectral component of the
 * gamma-ray sky model.
 * The power law function is defined as
 * \f[I(E)=integral (index+1)/(emax^{index+1}-emin^{index+1}) E^{index}\f]
 * for \f$index \ne -1\f$ and
 * \f[I(E)=integral / (\log(emax)-\log(emin)) E^{index}\f]
 * for \f$index = -1\f$, where
 * \f$integral\f$ is the integral flux between \f$emin\f$ and \f$emax\f$,
 * \f$index\f$ is the spectral index,
 * \f$emin\f$ is the lower energy limit, and
 * \f$emax\f$ is the upper energy limit.
 ***************************************************************************/
class GModelSpectralPlaw2 : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlaw2(void);
    explicit GModelSpectralPlaw2(const double& integral, const double& index);
    explicit GModelSpectralPlaw2(const GXmlElement& xml);
    GModelSpectralPlaw2(const GModelSpectralPlaw2& model);
    virtual ~GModelSpectralPlaw2(void);

    // Operators
    virtual GModelSpectralPlaw2& operator= (const GModelSpectralPlaw2& model);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralPlaw2* clone(void) const;
    virtual std::string          type(void) const { return "PowerLaw2"; }
    virtual double               eval(const GEnergy& srcEng) const;
    virtual double               eval_gradients(const GEnergy& srcEng) const;
    virtual double               flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(void) const;

    // Other methods
    double integral(void) const { return m_integral.Value(); }
    double index(void) const { return m_index.Value(); }
    double emin(void) const { return m_emin.Value(); }
    double emax(void) const { return m_emax.Value(); }

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
    mutable double  m_last_emin;       //!< Last lower energy limit (MeV)
    mutable double  m_last_emax;       //!< Last upper energy limit (MeV)
    mutable GEnergy m_last_energy;     //!< Last source energy
    mutable double  m_last_value;      //!< Last function value
    mutable double  m_last_g_integral; //!< Last integral flux gradient
    mutable double  m_last_g_index;    //!< Last spectral index gradient
};

#endif /* GMODELSPECTRALPLAW2_HPP */
