/***************************************************************************
 *     GObservations.i  -  Observations container class SWIG interface     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @brief Observation container class Python interface definition
 * @author J. Knodlseder
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
class GObservations {
public:
    // Constructors and destructors
    GObservations(void);
    GObservations(const GObservations& obs);
    virtual ~GObservations(void);

    // Methods
    void     clear(void);
    int      size(void) const;
    void     append(GObservation &obs);
    void     models(const GModels& models);
    void     models(const std::string& filename);
    GModels& models(void);
    void     optimize(GOptimizer& opt);
    double   npred(void) const;
};


/***********************************************************************//**
 * @brief GObservations class extension
 ***************************************************************************/
%extend GObservations {
    char *__str__() {
        return tochar(self->print());
    }
    GObservation& __getitem__(const int& index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(const int& index, const GObservation& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    GObservations copy() {
        return (*self);
    }
};
