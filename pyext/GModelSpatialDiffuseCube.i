/***************************************************************************
 *        GModelSpatialDiffuseCube.i - Spatial map cube model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseCube.i
 * @brief Spatial map cube model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialDiffuseCube.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialDiffuseCube
 *
 * @brief Spatial map cube model
 ***************************************************************************/
class GModelSpatialDiffuseCube  : public GModelSpatialDiffuse {
public:
    // Constructors and destructors
    GModelSpatialDiffuseCube(void);
    explicit GModelSpatialDiffuseCube(const GXmlElement& xml);
    GModelSpatialDiffuseCube(const std::string& filename,
                             const double&      value = 1.0);
    GModelSpatialDiffuseCube(const GSkymap&   map,
                             const GEnergies& energies,
                             const double&    value = 1.0);
    GModelSpatialDiffuseCube(const GModelSpatialDiffuseCube& model);
    virtual ~GModelSpatialDiffuseCube(void);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialDiffuseCube* clone(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GPhoton& photon) const;
    virtual double                    eval_gradients(const GPhoton& photon) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    int                        maps(void) const;
    int                        pixels(void) const;
    void                       load(const std::string& filename);
    double                     value(void) const;
    void                       value(const double& value);
    const std::string&         filename(void) const;
    void                       filename(const std::string& filename);
    const GSkymap&             cube(void) const;
    void                       cube(const GSkymap& cube);
    GEnergies                  energies(void);
    void                       energies(const GEnergies& energies);
    const GModelSpectralNodes& spectrum(void) const;
    void                       set_mc_cone(const GSkyDir& centre,
                                           const double&  radius);
};


/***********************************************************************//**
 * @brief GModelSpatialDiffuseCube class extension
 ***************************************************************************/
%extend GModelSpatialDiffuseCube {
    GModelSpatialDiffuseCube copy() {
        return (*self);
    }
};
