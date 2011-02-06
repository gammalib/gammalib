/***************************************************************************
 *           GModelData.hpp  -  Abstract virtual data model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelData.hpp
 * @brief GModelData class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELDATA_HPP
#define GMODELDATA_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GModel.hpp"
#include "GModelPar.hpp"
#include "GXmlElement.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModelData
 *
 * @brief Abstract virtual data model class interface defintion
 *
 * This abstract virtual base class implements methods to access model
 * parameters.
 *
 * @todo These methods could also be implemented on the GModel level.
 *       Before we can do this we also have to use std::vector<GModelPar*>
 *       in GModelSky. Once unified we can move this part to the base
 *       class.
 * @todo eval() and eval_gradients() should also be const. This has to be
 *       changed at the GModel level and has also to be propagated to the
 *       variuos derived classes.
 ***************************************************************************/
class GModelData : public GModel {

public:
    // Constructors and destructors
    GModelData(void);
    explicit GModelData(const GXmlElement& xml);
    GModelData(const GModelData& model);
    virtual ~GModelData(void);

    // Operators
    GModelPar&       operator() (int index);
    const GModelPar& operator() (int index) const;
    GModelData&      operator= (const GModelData& model);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModelData* clone(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual double      eval(const GEvent& event, const GObservation& obs) = 0;
    virtual double      eval_gradients(const GEvent& event, const GObservation& obs) = 0;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented pure virtual methods
    int size(void) const { return m_pars.size(); }

    // Other methods
    
protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelData& model);
    void         free_members(void);
    virtual void set_pointers(void) = 0;

    // Proteced data members
    std::vector<GModelPar*> m_pars;  //!< Pointers to all model parameters
};

#endif /* GMODELDATA_HPP */
