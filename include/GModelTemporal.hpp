/***************************************************************************
 *        GModelTemporal.hpp  -  Abstract temporal model base class        *
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
 * @file GModelTemporal.hpp
 * @brief GModelTemporal abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELTEMPORAL_HPP
#define GMODELTEMPORAL_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GModelTemporal
 *
 * @brief Abstract interface definition for the temporal model class.
 *
 * This class implements the temporal component of the factorised model.
 ***************************************************************************/
class GModelTemporal {

    // Friend classes
    friend class GModel;

public:
    // Constructors and destructors
    GModelTemporal(void);
    GModelTemporal(const GModelTemporal& model);
    virtual ~GModelTemporal(void);
 
    // Operators
    virtual GModelTemporal& operator= (const GModelTemporal& model);

    // Virtual methods
    virtual int        npars(void) const = 0;
    virtual GModelPar* par(int index) const = 0;
    virtual double     eval(const GTime& obsTime, const GSkyDir& srcDir,
                            const GEnergy& srcEng, const GTime& srcTime,
                            const GResponse& rsp, const GPointing& pnt) = 0;
    virtual double     eval_gradients(const GTime& obsTime, const GSkyDir& srcDir,
                                      const GEnergy& srcEng, const GTime& srcTime,
                                      const GResponse& rsp, const GPointing& pnt) = 0;
  
protected:
    // Protected methods
    void                    init_members(void);
    void                    copy_members(const GModelTemporal& model);
    void                    free_members(void);
    virtual GModelTemporal* clone(void) const = 0;
};

#endif /* GMODELTEMPORAL_HPP */
