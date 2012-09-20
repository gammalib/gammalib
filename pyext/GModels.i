/***************************************************************************
 *                   GModels.i  -  Model container class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
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
 * @file GModels.i
 * @brief GModels class SWIG interface.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModels.hpp"
#include "GModelSky.hpp"
#include "GModelDiffuseSource.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelPointSource.hpp"
#include "GModelData.hpp"
#include "GTools.hpp"
%}

/* __ Inform about base classes ___________________________________________*/
%import(module="opt") "GOptimizerPars.i";

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GModels
 *
 * @brief GModels class Python bindings.
 ***************************************************************************/
class GModels : public GOptimizerPars {
public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    explicit GModels(const std::string& filename);
    virtual ~GModels(void);
 
    // Methods
    void     clear(void);
    GModels* clone(void) const;
    int      size(void) const;
    void     append(const GModel& model);
    void     load(const std::string& filename);
    void     save(const std::string& filename) const;
    void     read(const GXml& xml);
    void     write(GXml& xml) const;
    double   eval(const GEvent& event, const GObservation& obs) const;
    double   eval_gradients(const GEvent& event, const GObservation& obs) const;
};


/***********************************************************************//**
 * @brief GModels class extension
 ***************************************************************************/
%extend GModels {
    char *__str__() {
        return tochar(self->print());
    }
    GModel& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    GModel& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GModel& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    void __setitem__(const std::string& name, const GModel& val) {
        (*self)[name] = val;
        return;
    }
    GModels copy() {
        return (*self);
    }
};
