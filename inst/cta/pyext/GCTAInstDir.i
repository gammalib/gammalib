/***************************************************************************
 *           GCTAInstDir.i  -  CTA instrument direction class              *
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
 * @file GCTAInstDir.i
 * @brief CTA instrument direction class Python interface definition
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
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Methods
    void         clear(void);
    GCTAInstDir* clone(void) const;
    void         dir(const GSkyDir& dir);
    void         radec(const double& ra, const double& dec);
    void         radec_deg(const double& ra, const double& dec);
    void         lb(const double& l, const double& b);
    void         lb_deg(const double& l, const double& b);
    void         rotate_deg(const double& phi, const double& theta);
    GSkyDir      dir(void) const;
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
    double       dist(GCTAInstDir& dir) const;
    double       dist_deg(GCTAInstDir& dir) const;
    double       posang(const GSkyDir& dir) const;
    double       posang_deg(const GSkyDir& dir) const;
    double       posang(const GCTAInstDir& dir) const;
    double       posang_deg(const GCTAInstDir& dir) const;
};


/***********************************************************************//**
 * @brief GCTAInstDir class extension
 ***************************************************************************/
%extend GCTAInstDir {
    GCTAInstDir copy() {
        return (*self);
    }
};
