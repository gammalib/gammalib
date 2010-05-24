/***************************************************************************
 *                    GModels.hpp  -  Model container class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
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
 * @brief GModels container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELS_HPP
#define GMODELS_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GOptimizerPars.hpp"
#include "GModel.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GModels
 *
 * @brief GModels class interface defintion.
 ***************************************************************************/
class GModels : public GOptimizerPars {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModels& models);

public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    ~GModels(void);
 
    // Operators
    GModels& operator= (const GModels& models);

    // Methods
    void   add(const GModel& model);
    double eval(GSkyDir* dir, GEnergy* energy, GTime* time);
    double eval_gradients(GSkyDir* dir, GEnergy* energy, GTime* time);
    int    size(void) const { return m_elements; }
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModels& models);
    void free_members(void);
    void set_pointers(void);

    // Proteced data members
    int     m_elements;          //!< Total number of models
    GModel* m_model;             //!< Array of models
};

#endif /* GMODELS_HPP */
