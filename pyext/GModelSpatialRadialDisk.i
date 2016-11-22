/***************************************************************************
 *        GModelSpatialRadialDisk.i - Radial disk source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Christoph Deil                              *
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
 * @file GModelSpatialRadialDisk.i
 * @brief Radial disk model class Python interface definition
 * @author Christoph Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialDisk.hpp"
%}

/**************************************************************************
 * @class GModelSpatialRadialDisk
 *
 * @brief Disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a disk source, i.e. constant surface brightness within some
 * radius and no emission outside.
 ***************************************************************************/
class GModelSpatialRadialDisk : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialDisk(void);
    GModelSpatialRadialDisk(const GSkyDir& dir, const double& radius);
    explicit GModelSpatialRadialDisk(const GXmlElement& xml);
    GModelSpatialRadialDisk(const GModelSpatialRadialDisk& model);
    virtual ~GModelSpatialRadialDisk(void);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpatialRadialDisk* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual double                   eval(const double&  theta,
                                          const GEnergy& energy,
                                          const GTime&   time,
                                          const bool&    gradients = false) const;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime&   time,
                                        GRan&          ran) const;
    virtual bool                     contains(const GSkyDir& dir,
                                              const double&  margin = 0.0) const;
    virtual double                   theta_max(void) const;
    virtual GSkyRegion*              region(void) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;

    // Other methods
    double radius(void) const;
    void   radius(const double& radius);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialDisk class extension
 *
 * The eval(GPhoton&, bool&) method need to be defined in the extension to
 * force swig to build also the interface for these methods that are
 * implemented in the base class only. It's not clear to me why these
 * methods are not inherited automatically. Maybe this could also be handled
 * by a %typemap(typecheck) construct.
 ***************************************************************************/
%extend GModelSpatialRadialDisk {
    GModelSpatialRadialDisk copy() {
        return (*self);
    }
    double eval(const GPhoton& photon, const bool& gradients) const {
        return self->GModelSpatialRadial::eval(photon, gradients);
    }
};
