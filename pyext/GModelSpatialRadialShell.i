/***************************************************************************
 *   GModelSpatialRadialShell.i - Radial spatial shell source model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2015 by Christoph Deil                              *
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
 * @file GModelSpatialRadialShell.i
 * @brief Radial spatial shell model class Python interface definition
 * @author Christoph Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialShell.hpp"
%}

/**************************************************************************
 * @class GModelSpatialRadialShell
 *
 * @brief Radial shell source model class
 ***************************************************************************/
class GModelSpatialRadialShell : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialShell(void);
    GModelSpatialRadialShell(const GSkyDir& dir,
                             const double&  radius,
                             const double&  width);
    explicit GModelSpatialRadialShell(const GXmlElement& xml);
    GModelSpatialRadialShell(const GModelSpatialRadialShell& model);
    virtual ~GModelSpatialRadialShell(void);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialRadialShell* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const double&  theta,
                                           const GEnergy& energy,
                                           const GTime& time) const;
    virtual double                    eval_gradients(const double& theta,
                                                     const GEnergy& energy,
                                                     const GTime& time) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual double                    theta_max(void) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    double  radius(void) const;
    double  width(void) const;
    void    radius(const double& radius);
    void    width(const double& width);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialShell class extension
 *
 * The eval(GSkyDir&) and eval_gradients(GSkyDir&) need to be defined in the
 * extension to force swig to build also the interface for these methods that
 * are implemented in the base class only. It's not clear to me why these
 * methods are not inherited automatically. Maybe this could also be handled
 * by a %typemap(typecheck) construct.
 ***************************************************************************/
%extend GModelSpatialRadialShell {
    GModelSpatialRadialShell copy() {
        return (*self);
    }
    double eval(const GPhoton& photon) const {
        return self->GModelSpatialRadial::eval(photon);
    }
    double eval_gradients(const GPhoton& photon) const {
        return self->GModelSpatialRadial::eval_gradients(photon);
    }
};
