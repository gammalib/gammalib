/***************************************************************************
 *              GModelSky.hpp - Abstract sky model base class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @file GModelSky.hpp
 * @brief Abstract sky model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSKY_HPP
#define GMODELSKY_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModel.hpp"
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPhoton.hpp"
#include "GPhotons.hpp"
#include "GRan.hpp"
#include "GVector.hpp"
#include "GXmlElement.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModelSky
 *
 * @brief Abstract sky model class
 *
 * This class implements a sky model that is factorised in a spatial,
 * a spectral and a temporal component.
 * The class has two methods for model evaluation. The eval() method
 * evaluates the model for a given observed photon direction, photon energy
 * and photon arrival time, given a reponse function and a pointing. The
 * eval_gradients() also evaluates the model, but also sets the gradients
 * of the model.
 * Protected methods are implemented to handle source parameter integrations
 * depending on the requirements. Integration of the model (method fct) is
 * first done over all sky directions (method spatial), then over all
 * energies (method spectral) and then over all times (method temporal).
 * The eval() and eval_gradients() methods call temporal() to perform the
 * nested integrations.
 ***************************************************************************/
class GModelSky : public GModel {

public:
    // Constructors and destructors
    GModelSky(void);
    explicit GModelSky(const GXmlElement& xml);
    explicit GModelSky(const GXmlElement& spatial, const GXmlElement& spectral);
    GModelSky(const GModelSky& model);
    virtual ~GModelSky(void);

    // Operators
    virtual GModelSky&  operator=(const GModelSky& model);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModelSky*  clone(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented pure virtual methods
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const;
    virtual void        read(const GXmlElement& xml);
    virtual void        write(GXmlElement& xml) const;

    // Other methods
    GModelSpatial*      spatial(void) const { return m_spatial; }
    GModelSpectral*     spectral(void) const { return m_spectral; }
    GModelTemporal*     temporal(void) const { return m_temporal; }
    double              value(const GSkyDir& srcDir, const GEnergy& srcEng,
                              const GTime& srcTime);
    GVector             gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                                  const GTime& srcTime);
    GPhotons            mc(const double& area, const GSkyDir& dir, const double& radius,
                           const GEnergy& emin, const GEnergy& emax,
                           const GTime& tmin, const GTime& tmax,
                           GRan& ran) const;

protected:
    // Protected methods
    void            init_members(void);
    void            copy_members(const GModelSky& model);
    void            free_members(void);
    void            set_pointers(void);
    GModelSpatial*  xml_spatial(const GXmlElement& spatial) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal* xml_temporal(const GXmlElement& temporal) const;
    double          spatial(const GEvent& event, const GEnergy& srcEng,
                            const GTime& srcTime, const GObservation& obs,
                            bool grad) const;
    double          spectral(const GEvent& event, const GTime& srcTime,
                             const GObservation& obs, bool grad) const;
    double          temporal(const GEvent& event, const GObservation& obs,
                             bool grad) const;
    bool            valid_model(void) const;
    std::string     print_model(void) const;

    // Proteced data members
    GModelSpatial*  m_spatial;       //!< Spatial model
    GModelSpectral* m_spectral;      //!< Spectral model
    GModelTemporal* m_temporal;      //!< Temporal model
};

#endif /* GMODELSKY_HPP */
