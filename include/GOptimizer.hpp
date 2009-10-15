/***************************************************************************
 *           GOptimizer.hpp  -  Abstract base class for optimizer          *
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
 * @file GOptimizer.hpp
 * @brief GOptimizer abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZER_HPP
#define GOPTIMIZER_HPP

/* __ Includes ___________________________________________________________ */
#include "GOptimizerPars.hpp"
#include "GOptimizerFunction.hpp"
//#include "GVector.hpp"


/***********************************************************************//**
 * @class GOptimizer
 *
 * @brief GOptimizer abstract base class interface defintion.
 ***************************************************************************/
class GOptimizer {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GOptimizer& opt);

public:

    // Constructors and destructors
    GOptimizer();
    GOptimizer(const GOptimizerFunction& fct, const GOptimizerPars &pars);
    GOptimizer(const GOptimizerFunction& fct, const GModels &models);
    GOptimizer(const GOptimizer& opt);
    virtual ~GOptimizer();

    // Operators
    virtual GOptimizer& operator= (const GOptimizer& opt);

    // Virtual methods
    virtual void optimize(void) const = 0;
 
    // Implement methods

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizer& opt);
    void free_members(void);
    
    // Protected data area
    const GOptimizerFunction* m_fct;     //!< Optimizer function
    GOptimizerPars            m_pars;    //!< Optimizer parameters

};

#endif /* GOPTIMIZER_HPP */
