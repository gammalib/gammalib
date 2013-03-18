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
 * This class implements a constant as the spectral component of the
 * gamma-ray sky model. The function is defined as
 * \f[I(E)=norm\f]
 * where
 * \f$norm\f$ is the normalization constant.
 ***************************************************************************/
class GModelSpectralConst : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralConst(void);
    explicit GModelSpectralConst(const GXmlElement& xml);
    GModelSpectralConst(const GModelSpectralConst& model);
    virtual ~GModelSpectralConst(void);

    // Operators
    virtual GModelSpectralConst& operator=(const GModelSpectralConst& model);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralConst* clone(void) const;
    virtual std::string          type(void) const { return "ConstantValue"; }
    virtual double               eval(const GEnergy& srcEng) const;
    virtual double               eval_gradients(const GEnergy& srcEng) const;
    virtual double               flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(void) const;

    // Other methods
    double norm(void) const { return m_norm.Value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_norm;      //!< Normalization factor
};

#endif /* GMODELSPECTRALCONST_HPP */
