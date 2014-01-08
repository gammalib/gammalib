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
#include "GPhoton.hpp"
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
 * This class implements a sky model that is factorized into a spatial,
 * a spectral and a temporal component. The factorization is given by
 *
 * \f[
 *    S(\vec{p}, E, t) = S_{\rm p}(\vec{p} | E, t) \,
 *                       S_{\rm E}(E | t) \,
 *                       S_{\rm t}(t)
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\vec{p} | E, t)\f$ is the spatial component,
 * - \f$S_{\rm E}(E | t)\f$ is the spectral component, and
 * - \f$S_{\rm t}(t)\f$ is the temporal component.
 *
 * The spatial component describes the energy and time dependent morphology
 * of the source. It satisfies
 * \f[
 *    \int_{\Omega} S_{\rm p}(\vec{p} | E, t) d\Omega = 1
 * \f]
 * for all \f$E\f$ and \f$t\f$, hence the spatial component does not
 * impact the spatially integrated spectral and temporal properties of the
 * source.
 *
 * The spectral component describes the spatially integrated time dependent
 * spectral distribution of the source. It satisfies
 * \f[
 *    \int_{E} S_{\rm E}(E | t) dE = \Phi
 * \f]
 * for all \f$t\f$, where \f$\Phi\f$ is the spatially and spectrally
 * integrated total source flux. The spectral component does not impact
 * the temporal properties of the integrated flux \f$\Phi\f$.
 *
 * The temporal component describes the temporal variation of the total
 * integrated flux \f$\Phi\f$ of the source.
 *
 * The class has two methods for model evaluation that evaluate the model
 * for a specific event, given an observation. The eval() method returns
 * the model value, the eval_gradients() returns the model value and sets
 * the analytical gradients for all model parameters.
 * Note that the eval() and eval_gradients() methods call protected
 * methods that handle time dispersion, energy dispersion and the point
 * spread function (spatial dispersion). Dispersion is handled by
 * integrating over the relevant intervals of the source properties.
 * Integration of the model is done by three nested integrals.
 * The outmost integral integrates over time (method integrate_time()),
 * the next integral integrates over energy (method integrate_energy()),
 * and the innermost integral integrates over the point spread function
 * (method integrate_dir()).
 *
 * The npred() method returns the integral over the model for a given
 * observed energy and time.
 *
 * The read() and write() methods allow reading of model information from
 * and writing to an XML element. The type() method returns the model type
 * that has been found in an XML element.
 *
 * The model factorization is implemented by the abstract model component
 * classes GModelSpatial, GModelSpectral and GModelTemporal. The GModelSky
 * class holds pointers to derived instances of these classes, which can
 * be accessed using the spatial(), spectral() and temporal() methods. Note
 * that these pointers can be NULL (for example if no model has been yet
 * defined), so the validity of the pointers needs to be checked before
 * using them.
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
    virtual bool        is_constant(void) const;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const;
    virtual double      npred(const GEnergy& obsEng,
                              const GTime& obsTime,
                              const GObservation& obs) const;
    virtual void        read(const GXmlElement& xml);
    virtual void        write(GXmlElement& xml) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GModelSpatial*      spatial(void) const;
    GModelSpectral*     spectral(void) const;
    GModelTemporal*     temporal(void) const;
    double              value(const GPhoton& photon);
    GVector             gradients(const GPhoton& photon);
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
    double          integrate_time(const GEvent& event,
                                   const GObservation& obs,
                                   bool grad) const;
    double          integrate_energy(const GEvent& event,
                                     const GTime& srcTime,
                                     const GObservation& obs,
                                     bool grad) const;
    double          integrate_dir(const GEvent& event,
                                  const GEnergy& srcEng,
                                  const GTime& srcTime,
                                  const GObservation& obs,
                                  bool grad) const;
    bool            valid_model(void) const;

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
 * @brief Signals if sky model is temporally constant
 *
 * @return True if sky model is temporally constant, false otherwise.
 *
 * Signals if the sky model is temporally constant. A temporally constant
 * model is a model that has a temporal component of type "Constant".
 ***************************************************************************/
inline
bool GModelSky::is_constant(void) const
{
    return (m_temporal != NULL && m_temporal->type() == "Constant");
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
