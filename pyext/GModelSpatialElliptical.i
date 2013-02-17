/***************************************************************************
 *   GModelSpatialElliptical.i - Abstract elliptical spatial model class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialElliptical.i
 * @brief Abstract elliptical model base class Python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialElliptical.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialElliptical
 *
 * @brief Abstract elliptical spatial model base class
 ***************************************************************************/
class GModelSpatialElliptical : public GModelSpatial {
public:
    // Constructors and destructors
    GModelSpatialElliptical(void);
    GModelSpatialElliptical(const GModelSpatialElliptical& model);
    explicit GModelSpatialElliptical(const GXmlElement& xml);
    virtual ~GModelSpatialElliptical(void);

    // Pure virtual methods
    virtual void                     clear(void) = 0;
    virtual GModelSpatialElliptical* clone(void) const = 0;
    virtual std::string              type(void) const = 0;
    virtual double                   eval(const double& theta,
                                          const double& posangle) const = 0;
    virtual double                   eval_gradients(const double& theta,
                                                    const double& posangle) const = 0;
    virtual GSkyDir                  mc(GRan& ran) const = 0;
    virtual double                   theta_max(void) const = 0;

    // Implemented virtual methods
    virtual double eval(const GSkyDir& srcDir) const;
    virtual double eval_gradients(const GSkyDir& srcDir) const;
    virtual void   read(const GXmlElement& xml);
    virtual void   write(GXmlElement& xml) const;

    // Other methods
    double  ra(void) const;
    double  dec(void) const;
    double  posangle(void) const;
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    posangle(const double& posangle);
};


/***********************************************************************//**
 * @brief GModelSpatialElliptical class extension
 ***************************************************************************/
%extend GModelSpatialElliptical {
};
