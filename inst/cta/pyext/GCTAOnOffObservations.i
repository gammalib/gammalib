/***************************************************************************
 *     GCTAOnOffObservations.i - CTA on-off observation container class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Pierrick Martin                                  *
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
 * @file GCTAOnOffObservations.i
 * @brief CTA on-off observation container class definition
 * @author Pierrick Martin
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAOnOffObservations.hpp"
%}


/***********************************************************************//**
 * @class GCTAOnOffObservations
 *
 * @brief ON/OFF Observation container class
 ***************************************************************************/
class GCTAOnOffObservations : public GContainer {
public:
    // Constructors and destructors
    GCTAOnOffObservations(void);
    GCTAOnOffObservations(const GCTAOnOffObservations& obs);
    explicit GCTAOnOffObservations(const std::string& filename);
    virtual ~GCTAOnOffObservations(void);

    // Methods
    void                        clear(void);
    GCTAOnOffObservations*      clone(void) const;
    GCTAOnOffObservation*       at(const int& index);
    const GCTAOnOffObservation* at(const int& index) const;
    int                         size(void) const;
    bool                        is_empty(void) const;
    GCTAOnOffObservation*       set(const int& index, const GCTAOnOffObservation& obs);
    GCTAOnOffObservation*       append(const GCTAOnOffObservation& obs);
    GCTAOnOffObservation*       insert(const int& index, const GCTAOnOffObservation& obs);
    void                        remove(const int& index);
    void                        reserve(const int& num);
    void                        extend(const GCTAOnOffObservations& obs);
    bool                        contains(const std::string& instrument,
                                         const std::string& id) const;
	void                        load(const std::string& filename);
    void                        save(const std::string& filename) const;
    void                        read(const GXml& xml);
    void                        write(GXml& xml) const;
	void                        models(const GModels& models);
    void                        models(const std::string& filename);
    const GModels&              models(void) const;	
    std::string                 print(const GChatter& chatter = NORMAL) const;
    void                        optimize(GOptimizer& opt);
	double                      npred(void) const;
	
	// Optimizer function access method
    const GCTAOnOffObservations::likelihood& function(void) const;
	
};

// Likelihood function
class likelihood : public GOptimizerFunction {
public:
    // Constructors and destructors
    likelihood(void);
    likelihood(GCTAOnOffObservations* obs);
    likelihood(const likelihood& fct);
    ~likelihood(void);
	
    // Implemented pure virtual base class methods
    double         value(void);
    double         npred(void) const;
    GVector*       gradient(void);
    GMatrixSparse* curvature(void);
	
    // Other methods
    void set(GCTAOnOffObservations* obs);
    void eval(const GOptimizerPars& pars);
};
%nestedworkaround GCTAOnOffObservations::likelihood;
%{
	// SWIG thinks that likelihood is a global class, so we need to trick the C++
	// compiler into understanding this so called global type.
	typedef GCTAOnOffObservations::likelihood likelihood;
	%}


/***********************************************************************//**
 * @brief GCTAOnOffObservations class extension
 ***************************************************************************/
%extend GCTAOnOffObservations {
    GCTAOnOffObservation* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GCTAOnOffObservation& val) {
        if (index>=0 && index < self->size()) {
            self->set(index, val);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GCTAOnOffObservations copy() {
        return (*self);
    }
};
