/***************************************************************************
 *                    GModels.i - Model container class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2017 by Juergen Knoedlseder                         *
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
 * @brief Model container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModels.hpp"
#include "GModel.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"
%}


/***********************************************************************//**
 * @class GModels
 *
 * @brief Model container class
 ***************************************************************************/
class GModels : public GContainer {
public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    explicit GModels(const GFilename& filename);
    virtual ~GModels(void);
 
    // Methods
    void           clear(void);
    GModels*       clone(void) const;
    std::string    classname(void) const;
    GModel*        at(const int& index);
    int            size(void) const;
    bool           is_empty(void) const;
    GModel*        set(const int& index, const GModel& model);
    GModel*        set(const std::string& name, const GModel& model);
    GModel*        append(const GModel& model);
    GModel*        insert(const int& index, const GModel& model);
    GModel*        insert(const std::string& name, const GModel& model);
    void           remove(const int& index);
    void           remove(const std::string& name);
    void           reserve(const int& num);
    void           extend(const GModels& models);
    bool           contains(const std::string& name) const;
    void           load(const GFilename& filename);
    void           save(const GFilename& filename) const;
    void           read(const GXml& xml);
    void           write(GXml& xml) const;
    int            npars(void) const;
    GOptimizerPars pars(void) const;
    double         eval(const GEvent& event,
                        const GObservation& obs,
                        const bool& gradients = false) const;
};


/***********************************************************************//**
 * @brief GModels class extension
 ***************************************************************************/
%extend GModels {
    GModel* __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Model index",
                                           index, self->size());
        }
    }
    GModel* __getitem__(const std::string& name) {
        return (*self)[name];
    }
    GModels* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PySlice_GetIndices((PySliceObject*)param, len, &start, &stop, &step) == 0) {
                GModels* models = new GModels;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        models->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        models->append(*(*self)[i]);
                    }
                }
                return models;
            }
            else {
                throw GException::invalid_argument("__getitem__(PyObject)",
                                                   "Invalid slice indices");
            }
        }
        else {
            throw GException::invalid_argument("__getitem__(PyObject)","");
        }
    }
    void __setitem__(const int& index, const GModel& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            self->set(index, val);
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            self->set(self->size()+index, val);
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "Model index",
                                           index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GModel& val) {
        self->set(name, val);
        return;
    }
    GModels copy() {
        return (*self);
    }
};
