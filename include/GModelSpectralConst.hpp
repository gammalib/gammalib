/***************************************************************************
 *         GModelSpectralConst.hpp - Spectral constant model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;


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
    GModelSpectralConst(const std::string& type, const std::string& value);
    explicit GModelSpectralConst(const GXmlElement& xml);
    explicit GModelSpectralConst(const double& value);
    GModelSpectralConst(const GModelSpectralConst& model);
    virtual ~GModelSpectralConst(void);

    // Operators
    virtual GModelSpectralConst& operator=(const GModelSpectralConst& model);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GModelSpectralConst* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime = GTime()) const;
    virtual double               eval_gradients(const GEnergy& srcEng,
                                                const GTime&   srcTime = GTime());
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
    double value(void) const;
    void   value(const double& value);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralConst& model);
    void free_members(void);

    // Protected members
    std::string m_type;     //!< Model type
    GModelPar   m_norm;     //!< Normalization factor
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralConst").
 ***************************************************************************/
inline
std::string GModelSpectralConst::classname(void) const
{
    return ("GModelSpectralConst");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the constant spectral model.
 ***************************************************************************/
inline
std::string GModelSpectralConst::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return model value
 *
 * @return Model value (ph/cm2/s/MeV).
 *
 * Returns the model value.
 ***************************************************************************/
inline
double GModelSpectralConst::value(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set model value
 *
 * @param[in] value Model value (ph/cm2/s/MeV).
 *
 * Sets the model value.
 ***************************************************************************/
inline
void GModelSpectralConst::value(const double& value)
{
    m_norm.value(value);
    return;
}

#endif /* GMODELSPECTRALCONST_HPP */
