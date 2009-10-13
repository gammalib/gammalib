/***************************************************************************
 *                    GModels.hpp  -  Model container class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
#include "GException.hpp"
#include "GModel.hpp"
#include "GModelPar.hpp"
#include <iostream>


/***********************************************************************//**
 * @class GModels
 *
 * @brief GModels class interface defintion.
 ***************************************************************************/
class GModels {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModels& models);

public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    ~GModels();
 
    // Operators
    GModels& operator= (const GModels& models);

    // Methods
    void       add(const GModel& model);
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModels& models);
    void free_members(void);
    void set_pointers(void);

    // Proteced data members
    int         m_elements;      //!< Total number of models
    int         m_npars;         //!< Total number of model parameters
    GModel*     m_model;         //!< Array of models
    GModelPar** m_par;           //!< Pointers to all model parameters
};

#endif /* GMODELS_HPP */
