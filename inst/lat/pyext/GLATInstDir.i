/***************************************************************************
 *           GLATInstDir.i  -  LAT instrument direction class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @brief GLATInstDir class python bindings
 * @author J. Knodlseder
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
    void         skydir(const GSkyDir& dir) { m_dir=dir; }
    void         radec(const double& ra, const double& dec) { m_dir.radec(ra,dec); }
    void         radec_deg(const double& ra, const double& dec) { m_dir.radec_deg(ra,dec); }
    void         lb(const double& l, const double& b) { m_dir.lb(l,b); }
    void         lb_deg(const double& l, const double& b) { m_dir.lb_deg(l,b); }
    GSkyDir      skydir(void) const { return m_dir; }
    double       l(void) const { return m_dir.l(); }
    double       l_deg(void) const { return m_dir.l_deg(); }
    double       b(void) const { return m_dir.b(); }
    double       b_deg(void) const { return m_dir.b_deg(); }
    double       ra(void) const { return m_dir.ra(); }
    double       ra_deg(void) const { return m_dir.ra_deg(); }
    double       dec(void) const { return m_dir.dec(); }
    double       dec_deg(void) const { return m_dir.dec_deg(); }
    double       dist(GSkyDir& dir) const { return m_dir.dist(dir); }
    double       dist_deg(GSkyDir& dir) const { return m_dir.dist_deg(dir); }
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


/***********************************************************************//**
 * @brief GLATInstDir type casts
 ***************************************************************************/
%inline %{
    GLATInstDir* cast_GLATInstDir(GInstDir* dir) {
        return dynamic_cast<GLATInstDir*>(dir);
    }
%}
