/***************************************************************************
 *        GModelTemporalFunc.i - Temporal file function model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GModelTemporalFunc.i
 * @brief File function temporal model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporalFunc.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporalFunc
 *
 * @brief File function temporal model class
 ***************************************************************************/
class GModelTemporalFunc : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalFunc(void);
    explicit GModelTemporalFunc(const GXmlElement& xml);
    GModelTemporalFunc(const GFilename& filename,
                       const double&    norm = 1.0);
    GModelTemporalFunc(const GModelTemporalFunc& model);
    virtual ~GModelTemporalFunc(void);

    // Implemented virtual base class methods
    virtual void                clear(void);
    virtual GModelTemporalFunc* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GTime& srcTime,
                                     const bool& gradients = false) const;
    virtual GTimes              mc(const double& rate, const GTime& tmin,
                                   const GTime& tmax, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
    double           norm(void) const;
    void             norm(const double& norm);
};


/***********************************************************************//**
 * @brief GModelTemporalFunc class extension
 ***************************************************************************/
%extend GModelTemporalFunc {
    GModelTemporalFunc copy() {
        return (*self);
    }
};
