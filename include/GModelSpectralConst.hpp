/***************************************************************************
 *         GModelSpectralConst.hpp - Spectral constant model class         *
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
 * @file GModelSpectralConst.hpp
 * @brief Constant spectral model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALCONST_HPP
#define GMODELSPECTRALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralConst
 *
 * @brief Constant spectral model class
 *
 * This class implements a constant spectrum. The model is defined by
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * \f${\tt m\_norm}\f$ is the normalization constant in units of
 * ph/cm2/s/MeV.
 ***************************************************************************/
class GModelSpectralConst : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralConst(void);
    explicit GModelSpectralConst(const GXmlElement& xml);
    explicit GModelSpectralConst(const double& norm);
    GModelSpectralConst(const GModelSpectralConst& model);
    virtual ~GModelSpectralConst(void);

    // Operators
    virtual GModelSpectralConst& operator=(const GModelSpectralConst& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralConst* clone(void) const;
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
    virtual std::string          print(void) const;

    // Other methods
    double norm(void) const;
    void   norm(const double& norm);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_norm;  //!< Normalization factor
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "ConstantValue".
 *
 * Returns the type of the constant spectral model.
 ***************************************************************************/
inline
std::string GModelSpectralConst::type(void) const
{
    return "ConstantValue";
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor (ph/cm2/s/MeV).
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelSpectralConst::norm(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set normalization factor 
 *
 * @param[in] norm Normalization factor (ph/cm2/s/MeV).
 *
 * Sets the normalization factor.
 ***************************************************************************/
inline
void GModelSpectralConst::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}

#endif /* GMODELSPECTRALCONST_HPP */
