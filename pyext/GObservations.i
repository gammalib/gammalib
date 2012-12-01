/***************************************************************************
 *             GObservations.i  -  Observations container class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
class GObservations : public GBase {
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
    void           append(GObservation& obs);
    void           load(const std::string& filename);
    void           save(const std::string& filename) const;
    void           read(const GXml& xml);
    void           write(GXml& xml) const;
    void           models(const GModels& models);
    void           models(const std::string& filename);
    GModels&       models(void) { return m_models; }
    void           optimize(GOptimizer& opt);
    double         npred(void) const;
};


/***********************************************************************//**
 * @brief GObservations class extension
 ***************************************************************************/
%extend GObservations {
    char *__str__() {
        return tochar(self->print());
    }
    GObservation& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GObservation& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
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
