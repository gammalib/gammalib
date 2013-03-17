/***************************************************************************
 *          GModelTemporal.i - Abstract temporal model base class          *
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
 * @file GModelTemporal.i
 * @brief Abstract temporal model base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporal.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporal
 *
 * @brief Abstract temporal model base class
 ***************************************************************************/
class GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporal(void);
    GModelTemporal(const GModelTemporal& model);
    virtual ~GModelTemporal(void);

    // Virtual methods
    virtual void            clear(void) = 0;
    virtual GModelTemporal* clone(void) const = 0;
    virtual std::string     type(void) const = 0;
    virtual double          eval(const GTime& srcTime) const = 0;
    virtual double          eval_gradients(const GTime& srcTime) const = 0;
    virtual GTimes          mc(const double& rate, const GTime& tmin,
                               const GTime& tmax, GRan& ran) const = 0;
    virtual void            read(const GXmlElement& xml) = 0;
    virtual void            write(GXmlElement& xml) const = 0;

    // Methods
    int  size(void) const;
    void autoscale(void);
};


/***********************************************************************//**
 * @brief GModelTemporal class extension
 ***************************************************************************/
%extend GModelTemporal {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar& __getitem__(const int& index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    GModelPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GModelPar& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    void __setitem__(const std::string& name, const GModelPar& val) {
        (*self)[name] = val;
        return;
    }
};
