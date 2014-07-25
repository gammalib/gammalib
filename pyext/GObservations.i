/***************************************************************************
 *              GObservations.i - Observations container class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @file GObservations.i
 * @brief Observations container class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GObservations.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GObservations
 *
 * @brief Observation container class
 ***************************************************************************/
class GObservations : public GContainer {
public:
    // Constructors and destructors
    GObservations(void);
    GObservations(const GObservations& obs);
    explicit GObservations(const std::string& filename);
    virtual ~GObservations(void);

    // Methods
    void           clear(void);
    GObservations* clone(void) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GObservation*  set(const int& index, const GObservation& obs);
    GObservation*  append(const GObservation& obs);
    GObservation*  insert(const int& index, const GObservation& obs);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GObservations& obs);
    bool           contains(const std::string& instrument,
                            const std::string& id) const;
    void           load(const std::string& filename);
    void           save(const std::string& filename) const;
    void           read(const GXml& xml);
    void           write(GXml& xml) const;
    void           models(const GModels& models);
    void           models(const std::string& filename);
    const GModels& models(void);
    void           optimize(GOptimizer& opt);
    void           eval(void);
    double         logL(void) const;
    double         npred(void) const;

    // Optimizer function access method
    const GObservations::likelihood& function(void) const;
};
    
// Likelihood function
class likelihood : public GOptimizerFunction {
public:
    // Constructors and destructors
    likelihood(void);
    likelihood(GObservations* obs);
    likelihood(const likelihood& fct);
    ~likelihood(void);

    // Implemented pure virtual base class methods
    virtual void           eval(const GOptimizerPars& pars);
    virtual double         value(void) const;
    virtual GVector*       gradient(void);
    virtual GMatrixSparse* curvature(void);

    // Other methods
    void   set(GObservations* obs);
    double npred(void) const;
};
%nestedworkaround GObservations::likelihood;
%{
// SWIG thinks that likelihood is a global class, so we need to trick the C++
// compiler into understanding this so called global type.
typedef GObservations::likelihood likelihood;
%}


/***********************************************************************//**
 * @brief GObservations class extension
 ***************************************************************************/
%extend GObservations {
    GObservation* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GObservation& val) {
        if (index>=0 && index < self->size()) {
            self->set(index, val);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GObservations copy() {
        return (*self);
    }
};
