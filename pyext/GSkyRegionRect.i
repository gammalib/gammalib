/***************************************************************************
 *                      GSkyRegionRect.i - Sky region class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Andreas Specovius                           *
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
 * @file GSkyRegionRect.i
 * @brief Rectangular sky region class interface definition
 * @author Andreas Specovius
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegionRect.hpp"
%}


/***********************************************************************//**
 * @class GSkyRegionRect
 *
 * @brief Interface for the rectangular sky region class
 ***************************************************************************/
class GSkyRegionRect : public GSkyRegion {
public:
    // Constructors and destructors
    GSkyRegionRect(void);
    GSkyRegionRect(const GSkyDir& centre
                  ,const double& w, const double& h
                  ,const double& posang_deg
                  );
    GSkyRegionRect(const double& ra, const double& dec
                  ,const double& w, const double& h
                  ,const double& posang_deg
                  );
    explicit GSkyRegionRect(const std::string& line);
    GSkyRegionRect(const GSkyRegionRect& region);
    virtual ~GSkyRegionRect(void);

    // Implemented methods
    void              clear(void);
    GSkyRegionRect*   clone(void) const;
    std::string       classname(void) const;
    double            posang(void) const;
    void              posang(const double& posang);
    double            posang_deg(void) const;
    void              posang_deg(const double& posang);
    double            width(void) const;
    void              width(const double& width);
    double            height(void) const;
    void              height(const double& height);
    const GSkyDir&    centre(void) const;
    void              centre(const GSkyDir& centre);
    void              centre(const double& ra,const double& dec);
    double            ra(void) const;
    double            dec(void) const;
    void              read(const std::string& line);
    std::string       write(void) const;
    bool              contains(const GSkyDir& dir) const;
    bool              contains(const GSkyRegion& reg) const;
    bool              overlaps(const GSkyRegion& reg) const;

    GSkyDir           transform_to_local(const GSkyDir& skydir) const;
};


/***********************************************************************//**
 * @brief GSkyRegionRect class extension
 ***************************************************************************/
%extend GSkyRegionRect {
    GSkyRegionRect copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.centre(), self.width(), self.height(), self.posang(), self.name())
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1], state[2], state[3])
        self.name(state[4])
}
};
