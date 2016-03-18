/***************************************************************************
 *     GModelSpatialRadialProfile.i - Radial profile source model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialRadialProfile.i
 * @brief Radial profile model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialProfile.hpp"
%}


/**************************************************************************
 * @class GModelSpatialRadialProfile
 *
 * @brief Radial profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for an arbitrary radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfile : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialProfile(void);
    explicit GModelSpatialRadialProfile(const GXmlElement& xml);
    GModelSpatialRadialProfile(const GModelSpatialRadialProfile& model);
    virtual ~GModelSpatialRadialProfile(void);

    // Pure virtual methods
    virtual void                        clear(void) = 0;
    virtual GModelSpatialRadialProfile* clone(void) const = 0;
    virtual std::string                 classname(void) const = 0;
    virtual std::string                 type(void) const = 0;
    virtual double                      theta_min(void) const = 0;
    virtual double                      theta_max(void) const = 0;

    // Implemented pure virtual base class methods
    virtual double  eval(const double&  theta,
                         const GEnergy& energy,
                         const GTime&   time) const;
    virtual double  eval_gradients(const double&  theta,
                                   const GEnergy& energy,
                                   const GTime&   time) const;
    virtual GSkyDir mc(const GEnergy& energy,
                       const GTime&   time,
                       GRan&          ran) const;
    virtual bool    contains(const GSkyDir& dir,
                             const double&  margin = 0.0) const;

    // Implement other methods
    int  num_nodes(void) const;
    void num_nodes(const int& number);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialProfile class extension
 *
 * The eval(GSkyDir&) and eval_gradients(GSkyDir&) need to be defined in the
 * extension to force swig to build also the interface for these methods that
 * are implemented in the base class only. It's not clear to me why these
 * methods are not inherited automatically. Maybe this could also be handled
 * by a %typemap(typecheck) construct.
 ***************************************************************************/
%extend GModelSpatialRadialProfile {
    double eval(const GPhoton& photon) const {
        return self->GModelSpatialRadial::eval(photon);
    }
    double eval_gradients(const GPhoton& photon) const {
        return self->GModelSpatialRadial::eval_gradients(photon);
    }
};
