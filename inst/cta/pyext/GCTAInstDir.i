/***************************************************************************
 *              GCTAInstDir.i - CTA instrument direction class             *
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
 * @file GCTAInstDir.i
 * @brief CTA instrument direction class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAInstDir.hpp"
%}


/***********************************************************************//**
 * @class GCTAInstDir
 *
 * @brief CTA instrument direction class
 ***************************************************************************/
class GCTAInstDir : public GInstDir {

public:
    // Constructors and destructors
    GCTAInstDir(void);
    explicit GCTAInstDir(const GSkyDir& dir);
    GCTAInstDir(const double& detx, const double& dety);
    GCTAInstDir(const GSkyDir& dir, const double& detx, const double& dety);
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GCTAInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual uint64_t     hash(void) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    const GSkyDir& dir(void);
    void           detx(const double &x);
    void           dety(const double &y);
    const double&  detx(void) const;
    const double&  dety(void) const;
    double         theta(void) const;
    double         phi(void) const;
    const bool&    has_dir(void) const;
    const bool&    has_detx(void) const;
    const bool&    has_dety(void) const;
};


/***********************************************************************//**
 * @brief GCTAInstDir class extension
 ***************************************************************************/
%extend GCTAInstDir {
    GCTAInstDir copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        if self.has_dir():
            dir = self.dir()
        else:
            dir = None
        if self.has_detx():
            detx = self.detx()
        else:
            detx = None
        if self.has_dety():
            dety = self.dety()
        else:
            dety = None
        state = (dir, detx, dety)
        return state
    def __setstate__(self, state):
        self.__init__()
        if state[0] is not None:
            self.dir(state[0])
        if state[1] is not None:
            self.detx(state[1])
        if state[2] is not None:
            self.dety(state[2])
}
};
