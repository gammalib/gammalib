/***************************************************************************
 *              GOptimizerPars.hpp  -  Parameter container class           *
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
 * @file GOptimizerPars.hpp
 * @brief GOptimizerPars container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZERPARS_HPP
#define GOPTIMIZERPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GModelPar.hpp"


/***********************************************************************//**
 * @class GOptimizerPars
 *
 * @brief GOptimizerPars container class interface defintion.
 ***************************************************************************/
class GOptimizerPars {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GOptimizerPars& pars);

public:
    // Constructors and destructors
    GOptimizerPars(void);
    GOptimizerPars(const GOptimizerPars& pars);
    ~GOptimizerPars();
 
    // Operators
    GOptimizerPars& operator= (const GOptimizerPars& pars);

    // Methods
    int        npars(void) const { return m_npars; }
    int        nfree(void) const;
    GModelPar* par(int index) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerPars& pars);
    void free_members(void);

    // Proteced data members
    int         m_npars;         //!< Total number of model parameters
    GModelPar** m_par;           //!< Pointers to model parameters
};

#endif /* GOPTIMIZERPARS_HPP */
