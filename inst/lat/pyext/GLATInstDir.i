/***************************************************************************
 *           GLATInstDir.i - Fermi/LAT instrument direction class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GLATInstDir.i
 * @brief Fermi/LAT instrument direction class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATInstDir.hpp"
%}


/***********************************************************************//**
 * @class GLATInstDir
 *
 * @brief Python bindings for the LAT instrument direction class
 ***************************************************************************/
class GLATInstDir : public GInstDir {
public:
    // Constructors and destructors
    GLATInstDir(void);
    explicit GLATInstDir(const GSkyDir& dir);
    GLATInstDir(const GLATInstDir& dir);
    virtual ~GLATInstDir(void);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GLATInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual double       hash(void) const;

    // Other methods
    void     dir(const GSkyDir& dir);
    GSkyDir& dir(void);
};


/***********************************************************************//**
 * @brief GLATInstDir class extension
 ***************************************************************************/
%extend GLATInstDir {
    GLATInstDir copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.dir(),)
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
