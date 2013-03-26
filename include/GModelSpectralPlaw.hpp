/***************************************************************************
 *         GModelSpectralPlaw.hpp - Spectral power law model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlaw.hpp
 * @brief Power law spectral model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALPLAW_HPP
#define GMODELSPECTRALPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralPlaw
 *
 * @brief Power law spectral model class
 *
 * This class implements a power law as the spectral component of the
 * gamma-ray sky model. The power law is defined as
 * \f[I(E)=norm (E/pivot)^{index}\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy, and
 * \f$index\f$ is the spectral index.
 ***************************************************************************/
class GModelSpectralPlaw : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlaw(void);
    explicit GModelSpectralPlaw(const double& norm, const double& index,
                                const double& pivot);
    explicit GModelSpectralPlaw(const GXmlElement& xml);
    GModelSpectralPlaw(const GModelSpectralPlaw& model);
    virtual ~GModelSpectralPlaw(void);

    // Operators
    virtual GModelSpectralPlaw& operator=(const GModelSpectralPlaw& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralPlaw* clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime&   srcTime) const;
    virtual double              eval_gradients(const GEnergy& srcEng,
                                               const GTime&   srcTime);
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
    virtual std::string         print(void) const;

    // Other methods
    double norm(void) const { return m_norm.value(); }
    double index(void) const { return m_index.value(); }
    double pivot(void) const { return m_pivot.value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralPlaw& model);
    void free_members(void);

    // Protected members
    GModelPar m_norm;            //!< Normalization factor
    GModelPar m_index;           //!< Spectral index
    GModelPar m_pivot;           //!< Pivot energy
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "PowerLaw".
 *
 * Returns the type of the spectral power law model.
 ***************************************************************************/
inline
std::string GModelSpectralPlaw::type(void) const
{
    return "PowerLaw";
}

#endif /* GMODELSPECTRALPLAW_HPP */
