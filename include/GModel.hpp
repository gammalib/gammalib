/***************************************************************************
 *                         GModel.hpp  -  Model class                      *
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
 * @file GModel.hpp
 * @brief GModel class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODEL_HPP
#define GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include <iostream>


/***********************************************************************//**
 * @class GModel
 *
 * @brief GModel class interface defintion.
 ***************************************************************************/
class GModel {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModel& model);

public:
    // Constructors and destructors
    GModel(void);
    GModel(const GModelSpatial& spatial, const GModelSpectral& spectral);
    GModel(const GModel& model);
    ~GModel();
 
    // Operators
    GModel& operator= (const GModel& model);

    // Methods
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModel& model);
    void free_members(void);
    void set_pointers(void);

    // Proteced data members
    int             m_npars;         //!< Total number of model parameters
    GModelPar**     m_par;           //!< Pointers to all model parameters
    GModelSpatial*  m_spatial;       //!< Spatial model
    GModelSpectral* m_spectral;      //!< Spectral model
    GModelTemporal* m_temporal;      //!< Temporal model
};

#endif /* GMODEL_HPP */
