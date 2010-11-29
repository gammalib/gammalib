/***************************************************************************
 *           GRoi.hpp  -  Region of interest abstract base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GRoi.hpp
 * @brief GRoi abstract base class definition.
 * @author J. Knodlseder
 */

#ifndef GROI_HPP
#define GROI_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GRoi
 *
 * @brief Abstract interface for the region of interest classes.
 *
 * The region of interest class holds instrument specific information about
 * the spatial region in detector or telescopes coordinates that is used
 * for an analysis. In particular, the definition of a region of interest
 * is required for an unbinned analysis.
 ***************************************************************************/
class GRoi {

  // Friend classes
  friend class GObservation;

public:
    // Constructors and destructors
    GRoi(void);
    GRoi(const GRoi& roi);
    virtual ~GRoi(void);

    // Operators
    virtual GRoi& operator= (const GRoi& roi);

    // Pure virtual methods
    virtual GRoi* clone(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GRoi& roi);
    void free_members(void);
};

#endif /* GROI_HPP */
