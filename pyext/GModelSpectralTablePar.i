/***************************************************************************
 *     GModelSpectralTablePar.i - Spectral table model parameter class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTable.i
 * @brief Spectral table model parameter class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralTablePar.hpp"
%}
%include "std_vector.i"
%template(DoubleVector) std::vector<double>;


/***********************************************************************//**
 * @class GModelSpectralTablePar
 *
 * @brief Spectral table model parameter class
 ***************************************************************************/
class GModelSpectralTablePar : public GBase {

public:
    // Constructors and destructors
    GModelSpectralTablePar(void);
    GModelSpectralTablePar(const GModelPar&           par,
                           const std::vector<double>& values);
    GModelSpectralTablePar(const GModelSpectralTablePar& par);
    virtual ~GModelSpectralTablePar(void);

    // Methods
    void                    clear(void);
    GModelSpectralTablePar* clone(void) const;
    std::string             classname(void) const;
    int                     size(void) const;
    bool                    is_empty(void) const;
    GModelPar&              par(void);
    const GNodeArray&       values(void) const;
    const int&              method(void) const;
    void                    method(const int& method);
};


/***********************************************************************//**
 * @brief GModelSpectralTablePar class extension
 ***************************************************************************/
%extend GModelSpectralTablePar {
    GModelSpectralTablePar copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.par(), self.value)
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1])
}
};
