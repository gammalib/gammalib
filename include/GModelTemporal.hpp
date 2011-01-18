/***************************************************************************
 *        GModelTemporal.hpp  -  Abstract temporal model base class        *
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
 * @file GModelTemporal.hpp
 * @brief GModelTemporal abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELTEMPORAL_HPP
#define GMODELTEMPORAL_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GTime.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelTemporal
 *
 * @brief Abstract interface definition for the temporal model class
 *
 * This class implements the temporal component of the factorised model.
 * The temporal component of the factorised model is supposed to describe the
 * relative variation of the source flux with respect to the mean value
 * that is given by the spectral component. Normally, this model will have
 * a mean value of 1.
 ***************************************************************************/
class GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporal(void);
    GModelTemporal(const GModelTemporal& model);
    virtual ~GModelTemporal(void);

    // Operators
    virtual GModelPar&       operator() (int index);
    virtual const GModelPar& operator() (int index) const;
    virtual GModelTemporal&  operator= (const GModelTemporal& model);

    // Virtual methods
    virtual void               clear(void) = 0;
    virtual GModelTemporal*    clone(void) const = 0;
    virtual int                size(void) const = 0;
    virtual std::string        type(void) const = 0;
    virtual double             eval(const GTime& srcTime) = 0;
    virtual double             eval_gradients(const GTime& srcTime) = 0;
    virtual std::vector<GTime> mc(const double& rate,
                                  const GTime& tmin, const GTime& tmax) = 0;
    virtual void               read(const GXmlElement& xml) = 0;
    virtual void               write(GXmlElement& xml) const = 0;
    virtual std::string        print(void) const = 0;

    // Implemented methods
    virtual void mc_seed(unsigned long long int seed) { m_ran.seed(seed); }
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporal& model);
    void free_members(void);

    // Pure virtual methods
    virtual GModelPar** par(void) = 0;

    // Protected members
    GRan m_ran;     //!< Random number generator for time generation
};

#endif /* GMODELTEMPORAL_HPP */
