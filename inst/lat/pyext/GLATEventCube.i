/***************************************************************************
 *              GLATEventCube.i - Fermi/LAT event cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @file GLATEventCube.i
 * @brief Fermi/LAT event cube class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventCube.hpp"
%}


/***********************************************************************//**
 * @class GLATEventCube
 *
 * @brief LAT event bin container class Python interface
 ***************************************************************************/
class GLATEventCube : public GEventCube {

public:
    // Constructors and destructors
    GLATEventCube(void);
    GLATEventCube(const GLATEventCube& cube);
    virtual ~GLATEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GLATEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename,
                                const bool& clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;

    // Other methods
    void              time(const GTime& time);
    void              map(const GSkymap& map);
    void              enodes(const GNodeArray& enodes);
    void              ontime(const double& ontime);
    const GTime&      time(void) const;
    const GSkymap&    map(void) const;
    const GNodeArray& enodes(void);
    const double&     ontime(void) const;
    int               nx(void) const;
    int               ny(void) const;
    int               npix(void) const;
    int               ebins(void) const;
    int               ndiffrsp(void) const;
    std::string       diffname(const int& index) const;
    GSkymap*          diffrsp(const int& index) const;
    double            maxrad(const GSkyDir& dir) const;
};


/***********************************************************************//**
 * @brief GLATEventCube class extension
 ***************************************************************************/
%extend GLATEventCube {
    GLATEventCube copy() {
        return (*self);
    }
    GLATEventCube(GEventCube* cube) {
        GLATEventCube* ptr = dynamic_cast<GLATEventCube*>(cube);
        if (ptr != NULL) {
            return (ptr->clone());
        }
        else {
            throw GException::bad_type("GLATEventCube(GEventCube*)",
                                       "GEventCube not of type GLATEventCube");
        }
    }
};
