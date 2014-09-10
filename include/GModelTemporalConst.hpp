/***************************************************************************
 *         GModelTemporalConst.hpp - Temporal constant model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GModelTemporalConst.hpp
 * @brief Constant temporal model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELTEMPORALCONST_HPP
#define GMODELTEMPORALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelTemporal.hpp"
#include "GTime.hpp"
#include "GTimes.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief Constant temporal model class
 *
 * This class implements a constant light curve. The model is defined by
 *
 * \f[
 *    S_{\rm t}(t) = {\tt m\_norm}
 * \f]
 *
 * where
 * \f${\tt m\_norm}\f$ is the normalization constant which by default is
 * set to unity. The parameter \f${\tt m\_norm}\f$ is not supposed to be
 * fitted.
 ***************************************************************************/
class GModelTemporalConst : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalConst(void);
    explicit GModelTemporalConst(const double& norm);
    GModelTemporalConst(const GModelTemporalConst& model);
    virtual ~GModelTemporalConst(void);

    // Operators
    virtual GModelTemporalConst& operator= (const GModelTemporalConst& model);

    // Implemented virtual base class methods
    virtual void                 clear(void);
    virtual GModelTemporalConst* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GTime& srcTime) const;
    virtual double               eval_gradients(const GTime& srcTime);
    virtual GTimes               mc(const double& rate, const GTime& tmin,
                                    const GTime& tmax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double norm(void) const;
    void   norm(const double& norm);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_norm;    //!< Constant
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelTemporalConst").
 ***************************************************************************/
inline
std::string GModelTemporalConst::classname(void) const
{
    return ("GModelTemporalConst");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "ConstantValue".
 *
 * Returns the type of the constant temporal model.
 ***************************************************************************/
inline
std::string GModelTemporalConst::type(void) const
{
    return "Constant";
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor.
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelTemporalConst::norm(void) const
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
void GModelTemporalConst::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}

#endif /* GMODELTEMPORALCONST_HPP */
