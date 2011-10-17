/***************************************************************************
 *         GModelTemporalConst.i  -  Temporal constant model class         *
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
 * @file GModelTemporalConst.i
 * @brief Constant temporal model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporalConst.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief Constant temporal model class
 ***************************************************************************/
class GModelTemporalConst  : public GModelTemporal {
public:
    // Constructors and destructors
    GModelTemporalConst(void);
    GModelTemporalConst(const GModelTemporalConst& model);
    virtual ~GModelTemporalConst(void);

    // Implemented virtual methods
    virtual void                 clear(void);
    virtual GModelTemporalConst* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GTime& srcTime) const;
    virtual double               eval_gradients(const GTime& srcTime) const;
    virtual GTimes               mc(const double& rate, const GTime& tmin,
                                    const GTime& tmax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;

    // Other methods
    double norm(void) const;
};


/***********************************************************************//**
 * @brief GModelTemporalConst class extension
 ***************************************************************************/
%extend GModelTemporalConst {
    GModelTemporalConst copy() {
        return (*self);
    }
};
