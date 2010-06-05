/***************************************************************************
 *        GModelTemporalConst.hpp  -  Temporal constant model class        *
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
 * @file GModelTemporalConst.hpp
 * @brief GModelTemporalConst class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELTEMPORALCONST_HPP
#define GMODELTEMPORALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelTemporal.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief Interface definition for the constant model class.
 *
 * This class implements the temporal component of the factorised model for
 * a model that is constant in time.
 ***************************************************************************/
class GModelTemporalConst  : public GModelTemporal {

    // Friend classes
    friend class GModel;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModelTemporalConst& model);

public:
    // Constructors and destructors
    GModelTemporalConst(void);
    GModelTemporalConst(const GModelTemporalConst& model);
    virtual ~GModelTemporalConst(void);
 
    // Operators
    GModelTemporalConst& operator= (const GModelTemporalConst& model);

    // Methods
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
    double     eval(const GTime& obsTime, const GSkyDir& srcDir,
                    const GEnergy& srcEng, const GTime& srcTime,
                    const GResponse& rsp, const GPointing& pnt);
    double     eval_gradients(const GTime& obsTime, const GSkyDir& srcDir,
                              const GEnergy& srcEng, const GTime& srcTime,
                              const GResponse& rsp, const GPointing& pnt);
  
protected:
    // Protected methods
    void                 init_members(void);
    void                 copy_members(const GModelTemporalConst& model);
    void                 free_members(void);
    GModelTemporalConst* clone(void) const;
    
    // Data area
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[1];          //!< Pointers to parameters
    GModelPar  m_norm;            //!< Constant
};

#endif /* GMODELTEMPORALCONST_HPP */
