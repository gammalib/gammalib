/***************************************************************************
 *          GCTAModelSpatial.i - Spatial model abstract base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatial.i
 * @brief Abstract spatial model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelSpatial.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelSpatial
 *
 * @brief Abstract spatial model class
 ***************************************************************************/
class GCTAModelSpatial : public GBase {

public:
    // Constructors and destructors
    GCTAModelSpatial(void);
    GCTAModelSpatial(const GCTAModelSpatial& model);
    virtual ~GCTAModelSpatial(void);

    // Operators
    virtual GCTAModelSpatial& operator=(const GCTAModelSpatial& model);
    virtual GModelPar&        operator[](const int& index);
    virtual GModelPar&        operator[](const std::string& name);

    // Pure virtual methods
    virtual void              clear(void) = 0;
    virtual GCTAModelSpatial* clone(void) const = 0;
    virtual std::string       classname(void) const = 0;
    virtual std::string       type(void) const = 0;
    virtual double            eval(const GCTAInstDir& dir,
                                   const GEnergy&     energy,
                                   const GTime&       time,
                                   const bool&        gradients = false) const = 0;
    virtual double            mc_max_value(const GCTAObservation& obs) const = 0;
    virtual void              read(const GXmlElement& xml) = 0;
    virtual void              write(GXmlElement& xml) const = 0;

    // Implemented virtual methods
    virtual GCTAInstDir       mc(const GEnergy&         energy,
                                 const GTime&           time,
                                 const GCTAObservation& obs,
                                 GRan&                  ran) const;

    // Methods
    int                       size(void) const;
    virtual double            npred(const GEnergy&      energy,
                                    const GTime&        time,
                                    const GObservation& obs) const;
};


/***********************************************************************//**
 * @brief GCTAModelSpatial class extension
 ***************************************************************************/
%extend GCTAModelSpatial {
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
