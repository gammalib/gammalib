/***************************************************************************
 *                         GModel.hpp  -  Model class                      *
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
 * @file GModel.hpp
 * @brief GModel class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODEL_HPP
#define GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GModel
 *
 * @brief GModel class interface defintion.
 *
 * This class implements a source model.
 ***************************************************************************/
class GModel {

    // Friend classes
    friend class GModels;
    
    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModel& model);

public:
    // Constructors and destructors
    GModel(void);
    GModel(const GModelSpatial& spatial, const GModelSpectral& spectral);
    GModel(const GModel& model);
    ~GModel(void);
 
    // Operators
    GModel& operator= (const GModel& model);

    // Methods
    std::string name(void) const { return m_name; }
    void        name(const std::string& name) { m_name=name; return; }
    int         npars(void) const { return m_npars; }
    GModelPar*  par(int index) const;
    double      eval(const GInstDir& dir, const GEnergy& energy, const GTime& time);
    double      eval_gradients(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                               const GResponse& rsp, const GPointing& pnt);
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModel& model);
    void free_members(void);
    void set_pointers(void);

    // Proteced data members
    std::string     m_name;          //!< Model name
    int             m_npars;         //!< Total number of model parameters
    GModelPar**     m_par;           //!< Pointers to all model parameters
    GModelSpatial*  m_spatial;       //!< Spatial model
    GModelSpectral* m_spectral;      //!< Spectral model
    GModelTemporal* m_temporal;      //!< Temporal model
};

#endif /* GMODEL_HPP */
