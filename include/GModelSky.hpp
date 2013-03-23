/***************************************************************************
 *                     GModelSky.hpp - Sky model class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @brief Sky model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSKY_HPP
#define GMODELSKY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModel.hpp"
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPhotons.hpp"
#include "GRan.hpp"
#include "GVector.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModelSky
 *
 * @brief Sky model class
 *
 * This class implements a sky model that is factorised into a spatial,
 * a spectral and a temporal component.
 *
 * The class has two methods for model evaluation. The eval() method
 * evaluates the model for a given observed photon direction, photon energy
 * and photon arrival time, given a reponse function and a pointing. The
 * eval_gradients() also evaluates the model, but also sets the gradients
 * of the model.
 *
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
    explicit GModelSky(const std::string& type);
    explicit GModelSky(const GXmlElement& xml);
    explicit GModelSky(const GXmlElement& spatial,
                       const GXmlElement& spectral);
    explicit GModelSky(const GXmlElement& spatial,
                       const GXmlElement& spectral,
                       const GXmlElement& temporal);
    explicit GModelSky(const GModelSpatial& spatial,
                       const GModelSpectral& spectral);
    explicit GModelSky(const GModelSpatial& spatial,
                       const GModelSpectral& spectral,
                       const GModelTemporal& temporal);
    GModelSky(const GModelSky& model);
    virtual ~GModelSky(void);

    // Operators
    virtual GModelSky& operator=(const GModelSky& model);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GModelSky*  clone(void) const;
    virtual std::string type(void) const;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const;
    virtual double      npred(const GEnergy& obsEng,
                              const GTime& obsTime,
                              const GObservation& obs) const;
    virtual void        read(const GXmlElement& xml);
    virtual void        write(GXmlElement& xml) const;
    virtual std::string print(void) const;

    // Other methods
    GModelSpatial*      spatial(void) const;
    GModelSpectral*     spectral(void) const;
    GModelTemporal*     temporal(void) const;
    double              value(const GSkyDir& srcDir,
                              const GEnergy& srcEng,
                              const GTime& srcTime);
    GVector             gradients(const GSkyDir& srcDir,
                                  const GEnergy& srcEng,
                                  const GTime& srcTime);
    GPhotons            mc(const double& area,
                           const GSkyDir& dir, const double& radius,
                           const GEnergy& emin, const GEnergy& emax,
                           const GTime& tmin, const GTime& tmax,
                           GRan& ran) const;

protected:
    // Protected methods
    void            init_members(void);
    void            copy_members(const GModelSky& model);
    void            free_members(void);
    void            set_pointers(void);
    void            set_type(void);
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
    std::string     m_type;       //!< Model type
    GModelSpatial*  m_spatial;    //!< Spatial model
    GModelSpectral* m_spectral;   //!< Spectral model
    GModelTemporal* m_temporal;   //!< Temporal model
};


/***********************************************************************//**
 * @brief Return sky model type
 *
 * @return Sky model type.
 *
 * Returns the type of the sky model. The type is an arbitrary string that
 * is used in the XML declaration of the model to describe the model type.
 ***************************************************************************/
inline
std::string GModelSky::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return spatial model component
 *
 * @return Pointer to spatial model component.
 *
 * Returns a pointer to the spatial model component of the model. The pointer
 * is of type GModelSpatial. Note that a NULL pointer may be returned if the
 * sky model has no spatial model component.
 ***************************************************************************/
inline
GModelSpatial* GModelSky::spatial(void) const
{
    return (m_spatial);
}


/***********************************************************************//**
 * @brief Return spectral model component
 *
 * @return Pointer to spectral model component.
 *
 * Returns a pointer to the spectral model component of the model. The
 * pointer is of type GModelSpectral. Note that a NULL pointer may be
 * returned if the sky model has no spectral model component.
 ***************************************************************************/
inline
GModelSpectral* GModelSky::spectral(void) const
{
    return (m_spectral);
}


/***********************************************************************//**
 * @brief Return temporal model component
 *
 * @return Pointer to temporal model component.
 *
 * Returns a pointer to the temporal model component of the model. The
 * pointer is of type GModelTemporal. Note that a NULL pointer may be
 * returned if the sky model has no temporal model component.
 ***************************************************************************/
inline
GModelTemporal* GModelSky::temporal(void) const
{
    return (m_temporal);
}

#endif /* GMODELSKY_HPP */
