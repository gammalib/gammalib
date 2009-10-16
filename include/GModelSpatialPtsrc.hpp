/***************************************************************************
 *      GModelSpatialPtsrc.hpp  -  Spatial point source model class        *
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
 * @file GModelSpatialPtsrc.hpp
 * @brief GModelSpatialPtsrc class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALPTSRC_HPP
#define GMODELSPATIALPTSRC_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GSkyDir.hpp"
#include <iostream>


/***********************************************************************//**
 * @class GModelSpatialPtsrc
 *
 * @brief Powerlaw interface definition
 ***************************************************************************/
class GModelSpatialPtsrc  : public GModelSpatial {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModelSpatialPtsrc& model);

public:
    // Constructors and destructors
    GModelSpatialPtsrc(void);
    GModelSpatialPtsrc(GSkyDir dir);
    GModelSpatialPtsrc(const GModelSpatialPtsrc& model);
    virtual ~GModelSpatialPtsrc();
 
    // Operators
    GModelSpatialPtsrc& operator= (const GModelSpatialPtsrc& model);

    // Methods
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
    void       eval_gradients(void);
    double     ra(void) const { return m_ra.real_value(); }
    double     dec(void) const { return m_dec.real_value(); }
  
protected:
    // Protected methods
    void                init_members(void);
    void                copy_members(const GModelSpatialPtsrc& model);
    void                free_members(void);
    GModelSpatialPtsrc* clone(void) const;
    
    // Data area
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[2];          //!< Pointers to parameters
    GModelPar  m_ra;              //!< Right Ascension (deg)
    GModelPar  m_dec;             //!< Declination (deg)
};

#endif /* GMODELSPATIALPTSRC_HPP */
