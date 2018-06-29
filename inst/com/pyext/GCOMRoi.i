/***************************************************************************
 *               GCOMRoi.i - COMPTEL region of interest class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2018 by Juergen Knoedlseder                         *
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
 * @file GCOMRoi.i
 * @brief COMPTEL region of interest class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMRoi.hpp"
%}


/***********************************************************************//**
 * @class GCOMRoi
 *
 * @brief COMPTEL region of interest class
 ***************************************************************************/
class GCOMRoi : public GRoi {

public:
    // Constructors and destructors
    GCOMRoi(void);
    GCOMRoi(const GCOMRoi& roi);
    GCOMRoi(const GCOMInstDir& centre, const double& radius,
            const double& phibar_min, const double& phibar_max);
    virtual ~GCOMRoi(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMRoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;

    // Other methods
    const GCOMInstDir& centre(void) const;
    const double&      radius(void) const;
    const double&      phibar_min(void) const;
    const double&      phibar_max(void) const;
    void               centre(const GCOMInstDir& centre);
    void               radius(const double& radius);
    void               phibar_min(const double& phibar_min);
    void               phibar_max(const double& phibar_max);
};


/***********************************************************************//**
 * @brief GCOMRoi class extension
 ***************************************************************************/
%extend GCOMRoi {
    GCOMRoi copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.centre(), self.radius(), self.phibar_min(), self.phibar_max())
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1], state[2], state[3])
}
};
