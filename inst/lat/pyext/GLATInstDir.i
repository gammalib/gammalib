/***************************************************************************
 *          GLATInstDir.i - Fermi-LAT instrument direction class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @brief Fermi-LAT instrument direction class interface definition
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
    GLATInstDir(const GLATInstDir& dir);
    virtual ~GLATInstDir(void);

    // Methods
    void         clear(void);
    GLATInstDir* clone(void) const;
    void         skydir(const GSkyDir& dir);
    void         radec(const double& ra, const double& dec);
    void         radec_deg(const double& ra, const double& dec);
    void         lb(const double& l, const double& b);
    void         lb_deg(const double& l, const double& b);
    GSkyDir      skydir(void) const;
    double       l(void) const;
    double       l_deg(void) const;
    double       b(void) const;
    double       b_deg(void) const;
    double       ra(void) const;
    double       ra_deg(void) const;
    double       dec(void) const;
    double       dec_deg(void) const;
    double       dist(GSkyDir& dir) const;
    double       dist_deg(GSkyDir& dir) const;
    double       dist(GLATInstDir& dir) const;
    double       dist_deg(GLATInstDir& dir) const;
};


/***********************************************************************//**
 * @brief GLATInstDir class extension
 ***************************************************************************/
%extend GLATInstDir {
    GLATInstDir copy() {
        return (*self);
    }
};
