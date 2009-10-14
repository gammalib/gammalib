/***************************************************************************
 *              test_GOptimizer.hpp  -  test GOptimizer class              *
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

#ifndef TEST_GOPTIMIZER_HPP
#define TEST_GOPTIMIZER_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception


/***********************************************************************//**
 * @class MyOptFct
 *
 * @brief MyOptFct class interface defintion.
 ***************************************************************************/
class MyOptFct : public GOptimizerFunction {

public:

    // Constructors and destructors
    MyOptFct() { init_members(); return; }
    MyOptFct(const MyOptFct& fct) { init_members(); copy_members(fct); return; }
    virtual ~MyOptFct() { free_members(); return; }

    // Operators
    MyOptFct& operator= (const MyOptFct& fct);

    // Methods
    void     first_item(void);     //!< Move to first item
    bool     next_item(void);      //!< Move to next element (true if end)
    double   get_data(void) { return m_event->counts(); }
    double   get_model(void) { return m_event->counts(); }
    GVector  get_gradient(void) ;  //!< Get gradient reference
    void     set_data(const GData& data);
 
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const MyOptFct& fct);
    void free_members(void) { return; }
    
    // Protected data area
    GData           m_data;
    GData::iterator m_event;   //!< Iterator on actual event
    GData::iterator m_end;     //!< Iterator on end (speeds up iterating)

};

#endif /* TEST_GOPTIMIZER_HPP */
