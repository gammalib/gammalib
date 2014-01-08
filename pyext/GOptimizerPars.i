/***************************************************************************
 *     GOptimizerPars.i - Abstract optimizer parameter container class     *
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
 * @file GOptimizerPars.i
 * @brief Abstract optimizer parameters base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizerPars.hpp"
%}


/***********************************************************************//**
 * @class GOptimizerPars
 *
 * @brief Optimizer parameter container class
 ***************************************************************************/
class GOptimizerPars : public GContainer {
public:
    // Constructors and destructors
    GOptimizerPars(void);
    explicit GOptimizerPars(const int& number);
    GOptimizerPars(const GOptimizerPars& pars);
    virtual ~GOptimizerPars(void);

    // Methods
    void            clear(void);
    GOptimizerPars* clone(void) const;
    GOptimizerPar*  at(const int& index);
    int             size(void) const;
    bool            is_empty(void) const;
    int             nfree(void) const;
    GOptimizerPar*  set(const int& index, const GOptimizerPar& par);
    GOptimizerPar*  set(const std::string& name, const GOptimizerPar& par);
    void            attach(GOptimizerPar *par);
    void            attach(const int& index, GOptimizerPar* par);
    void            attach(const std::string& name, GOptimizerPar* par);
    GOptimizerPar*  insert(const int& index, const GOptimizerPar& par);
    GOptimizerPar*  insert(const std::string& name, const GOptimizerPar& par);
    void            remove(const int& index);
    void            remove(const std::string& name);
    void            reserve(const int& num);
    void            extend(const GOptimizerPars& pars);
    bool            contains(const std::string& name) const;
};


/***********************************************************************//**
 * @brief GOptimizerPars class extension
 ***************************************************************************/
%extend GOptimizerPars {
    GOptimizerPar* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Parameter index",
                                           index, self->size());
        }
    }
    GOptimizerPar* __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, GOptimizerPar* par) {
        if (index >= 0 && index < self->size()) {
            self->attach(index, par);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "Parameter index",
                                           index, self->size());
        }
    }
    void __setitem__(const std::string& name, GOptimizerPar* par) {
        self->attach(name, par);
        return;
    }
    GOptimizerPars copy() {
        return (*self);
    }
};
