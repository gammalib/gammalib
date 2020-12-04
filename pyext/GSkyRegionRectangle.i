/***************************************************************************
 *          GSkyRegionRectangle.i - Rectangular sky region class           *
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
 * @file GSkyRegionRectangle.i
 * @brief Rectangular sky region class interface definition
 * @author Andreas Specovius
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegionRectangle.hpp"
%}


/***********************************************************************//**
 * @class GSkyRegionRectangle
 *
 * @brief Interface for the rectangular sky region class
 ***************************************************************************/
class GSkyRegionRectangle : public GSkyRegion {

public:
    // Constructors and destructors
    GSkyRegionRectangle(void);
    GSkyRegionRectangle(const GSkyDir& centre,
                        const double&  width,
                        const double&  height,
                        const double&  posang);
    GSkyRegionRectangle(const double& ra,
                        const double& dec,
                        const double& width,
                        const double& height,
                        const double& posang);
    explicit GSkyRegionRectangle(const std::string& line);
    GSkyRegionRectangle(const GSkyRegionRectangle& region);
    virtual ~GSkyRegionRectangle(void);

    // Implemented methods
    void                 clear(void);
    GSkyRegionRectangle* clone(void) const;
    std::string          classname(void) const;
    const GSkyDir&       centre(void) const;
    void                 centre(const GSkyDir& centre);
    void                 centre(const double& ra,const double& dec);
    double               ra(void) const;
    double               dec(void) const;
    const double&        width(void) const;
    void                 width(const double& width);
    const double&        height(void) const;
    void                 height(const double& height);
    const double&        posang(void) const;
    void                 posang(const double& posang);
    GSkyDir              corner(const int& index) const;
    void                 read(const std::string& line);
    std::string          write(void) const;
    bool                 contains(const GSkyDir& dir) const;
    bool                 contains(const GSkyRegion& reg) const;
    bool                 overlaps(const GSkyRegion& reg) const;
};


/***********************************************************************//**
 * @brief GSkyRegionRectangle class extension
 ***************************************************************************/
%extend GSkyRegionRectangle {
    GSkyRegionRectangle copy() {
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
