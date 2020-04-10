/***************************************************************************
 *         GSPIInstDir.i - INTEGRAL/SPI instrument direction class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIInstDir.i
 * @brief INTEGRAL/SPI instrument direction class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSPIInstDir.hpp"
%}


/***********************************************************************//**
 * @class GSPIInstDir
 *
 * @brief INTEGRAL/SPI instrument direction class
 ***************************************************************************/
class GSPIInstDir : public GInstDir {

public:
    // Constructors and destructors
    GSPIInstDir(void);
    GSPIInstDir(const GSkyDir& dir, const int& detid);
    GSPIInstDir(const GSPIInstDir& dir);
    virtual ~GSPIInstDir(void);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GSPIInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual u_int64_t    hash(void) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    const GSkyDir& dir(void) const;
    void           detid(const int& detid);
    const int&     detid(void) const;
};


/***********************************************************************//**
 * @brief GSPIInstDir class extension
 ***************************************************************************/
%extend GSPIInstDir {
    GSPIInstDir copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.dir(), self.detid())
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1])
}
};
