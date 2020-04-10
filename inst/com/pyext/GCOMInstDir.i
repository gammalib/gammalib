/***************************************************************************
 *           GCOMInstDir.i - COMPTEL instrument direction class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2020 by Juergen Knoedlseder                         *
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
 * @file GCOMInstDir.hpp
 * @brief COMPTEL instrument direction class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMInstDir.hpp"
%}


/***********************************************************************//**
 * @class GCOMInstDir
 *
 * @brief Interface for the COMPTEL instrument direction class
 ***************************************************************************/
class GCOMInstDir : public GInstDir {
public:
    // Constructors and destructors
    GCOMInstDir(void);
    GCOMInstDir(const GCOMInstDir& dir);
    GCOMInstDir(const GSkyDir& dir, const double& phibar);
    virtual ~GCOMInstDir(void);

    // Methods
    virtual void         clear(void);
    virtual GCOMInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual uint64_t     hash(void) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    const GSkyDir& dir(void) const;
    void           phibar(const double& phibar);
    const double&  phibar(void) const;
};


/***********************************************************************//**
 * @brief GCOMInstDir class extension
 ***************************************************************************/
%extend GCOMInstDir {
    GCOMInstDir copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.dir(), self.phibar())
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1])
}
};
