/***************************************************************************
 *                         GModel.hpp  -  Model class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModel.hpp
 * @brief GModel class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODEL_HPP
#define GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPhoton.hpp"
#include "GRan.hpp"
#include "GVector.hpp"
#include "GXmlElement.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModel
 *
 * @brief GModel class interface defintion.
 *
 * This class implements a source model that is factorised in a spatial,
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
 *
 * @todo Do we need a separate eval_gradients() method or can we just use
 * eval() since we anyway need gradients for function optimization?
 ***************************************************************************/
class GModel {

    // Friend classes
    friend class GModels;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModel& model);
    friend GLog&         operator<< (GLog& log, const GModel& model);

public:
    // Constructors and destructors
    GModel(void);
    GModel(const GModel& model);
    explicit GModel(const GModelSpatial& spatial, const GModelSpectral& spectral);
    explicit GModel(const GXmlElement& spatial, const GXmlElement& spectral);
    virtual ~GModel(void);

    // Operators
    GModelPar&       operator() (int index);
    const GModelPar& operator() (int index) const;
    GModel&          operator= (const GModel& model);

    // Methods
    void            clear(void);
    GModel*         clone(void) const;
    int             size(void) const { return m_npars; }
    std::string     name(void) const { return m_name; }
    void            name(const std::string& name) { m_name=name; return; }
    void            instruments(const std::string& instruments);
    GModelSpatial*  spatial(void) const { return m_spatial; }
    GModelSpectral* spectral(void) const { return m_spectral; }
    GModelTemporal* temporal(void) const { return m_temporal; }
    double          value(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime);
    GVector         gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                              const GTime& srcTime);
    double          eval(const GEvent& event, const GObservation& obs);
    double          eval_gradients(const GEvent& event, const GObservation& obs);
    GPhotons        mc(const double& area, const GSkyDir& dir, const double& radius,
                       const GEnergy& emin, const GEnergy& emax,
                       const GTime& tmin, const GTime& tmax,
                       GRan& ran);
    bool            isvalid(const std::string& name) const;
    std::string     print(void) const;

protected:
    // Protected methods
    void            init_members(void);
    void            copy_members(const GModel& model);
    void            free_members(void);
    void            set_pointers(void);
    GModelSpatial*  xml_spatial(const GXmlElement& spatial) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    double          spatial(const GEvent& event, const GEnergy& srcEng, const GTime& srcTime,
                            const GObservation& obs, bool grad = false);
    double          spectral(const GEvent& event, const GTime& srcTime,
                             const GObservation& obs, bool grad = false);
    double          temporal(const GEvent& event,
                             const GObservation& obs, bool grad = false);
    bool            valid_model(void) const;

    // Proteced data members
    std::string              m_name;          //!< Model name
    std::vector<std::string> m_instruments;   //!< Instruments to which model applies
    int                      m_npars;         //!< Total number of model parameters
    GModelPar**              m_par;           //!< Pointers to all model parameters
    GModelSpatial*           m_spatial;       //!< Spatial model
    GModelSpectral*          m_spectral;      //!< Spectral model
    GModelTemporal*          m_temporal;      //!< Temporal model
};

#endif /* GMODEL_HPP */
