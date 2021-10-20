/***************************************************************************
 *            GModelSpatial.i - Spatial model abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2021 by Juergen Knoedlseder                         *
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
 * @file GModelSpatial.i
 * @brief Abstract spatial model base class interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatial.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatial
 *
 * @brief Abstract spatial model base class
 ***************************************************************************/
class GModelSpatial : public GBase {

public:
    // Constructors and destructors
    GModelSpatial(void);
    GModelSpatial(const GModelSpatial& model);
    virtual ~GModelSpatial(void);

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GModelSpatial* clone(void) const = 0;
    virtual std::string    classname(void) const = 0;
    virtual double         eval(const GPhoton& photon,
                                const bool& gradients = false) const = 0;
    virtual GSkyDir        mc(const GEnergy& energy, const GTime& time,
                              GRan& ran) const = 0;
    virtual double         mc_norm(const GSkyDir& dir,
                                   const double&  radius) const = 0;
    virtual bool           contains(const GSkyDir& dir,
                                    const double&  margin = 0.0) const = 0;
    virtual void           read(const GXmlElement& xml) = 0;
    virtual void           write(GXmlElement& xml) const = 0;

    // Virtual methods
    virtual double flux(const GSkyRegion&    region,
                        const GEnergy&       srcEng  = GEnergy(),
                        const GTime&         srcTime = GTime(),
                        const GPolarization& srcPol = GPolarization()) const;

    // Methods
    std::string       type(void) const;
    void              type(const std::string& type);
    GModelPar&        at(const int& index);
    bool              has_par(const std::string& name) const;
    bool              has_free_pars(void) const;
    int               size(void) const;
    void              autoscale(void);
    const GSkyRegion* region(void) const;
};


/***********************************************************************//**
 * @brief GModelSpatial class extension
 ***************************************************************************/
%extend GModelSpatial {
    GModelPar& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Model index",
                                           index, self->size());
        }
    }
    GModelPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GModelPar& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Model index",
                                           index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GModelPar& val) {
        (*self)[name] = val;
        return;
    }
};
