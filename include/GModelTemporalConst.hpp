/***************************************************************************
 *        GModelTemporalConst.hpp  -  Temporal constant model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
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
 * @author J. Knodlseder
 */

#ifndef GMODELTEMPORALCONST_HPP
#define GMODELTEMPORALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelTemporal.hpp"
#include "GTime.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief Constant temporal model class
 *
 * This class implements the temporal component of the factorised model for
 * a model that is constant in time. It has a single parameter norm that is
 * used for scaling, and that normally is set to 1. The parameter is not
 * supposed to be fitted.
 ***************************************************************************/
class GModelTemporalConst : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalConst(void);
    GModelTemporalConst(const GModelTemporalConst& model);
    virtual ~GModelTemporalConst(void);

    // Operators
    virtual GModelTemporalConst& operator= (const GModelTemporalConst& model);

    // Implemented virtual methods
    virtual void                 clear(void);
    virtual GModelTemporalConst* clone(void) const;
    virtual std::string          type(void) const { return "Constant"; }
    virtual double               eval(const GTime& srcTime) const;
    virtual double               eval_gradients(const GTime& srcTime) const;
    virtual GTimes               mc(const double& rate, const GTime& tmin,
                                    const GTime& tmax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(void) const;

    // Other methods
    double norm(void) const { return m_norm.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_norm;    //!< Constant
};

#endif /* GMODELTEMPORALCONST_HPP */
