/***************************************************************************
 *                    GModels.hpp  -  Model container class                *
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
 * @file GModels.hpp
 * @brief Model container class definition
 * @author J. Knodlseder
 */

#ifndef GMODELS_HPP
#define GMODELS_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GOptimizerPars.hpp"
#include "GModel.hpp"
#include "GXml.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModels
 *
 * @brief Model container class interface defintion.
 *
 * This container class collects models of gamma-ray data that are used
 * for maximum likelihood fitting. It derives from the optimizer parameter
 * class GOptimizerPars.
 ***************************************************************************/
class GModels : public GOptimizerPars {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModels& models);
    friend GLog&         operator<<(GLog& log,        const GModels& models);

public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    explicit GModels(const std::string& filename);
    virtual ~GModels(void);

    // Operators
    GModel*       operator() (int index);
    const GModel* operator() (int index) const;
    GModels&      operator= (const GModels& models);

    // Methods
    void          clear(void);
    GModels*      clone(void) const;
    int           size(void) const { return m_models.size(); }
    void          append(const GModel& model);
    void          load(const std::string& filename);
    void          save(const std::string& filename) const;
    void          read(const GXml& xml);
    void          write(GXml& xml) const;
    double        eval(const GEvent& event, const GObservation& obs);
    double        eval_gradients(const GEvent& event, const GObservation& obs);
    std::string   print(void) const;

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModels& models);
    void          free_members(void);
    void          set_pointers(void);

    // Proteced members
    std::vector<GModel*> m_models;  //!< List of models
};

#endif /* GMODELS_HPP */
