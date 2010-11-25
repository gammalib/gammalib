/***************************************************************************
 *              GEventAtom.hpp  -  Event atom abstract base class          *
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
 * @file GEventAtom.hpp
 * @brief GEventAtom abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTATOM_HPP
#define GEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvent.hpp"
#include "GModels.hpp"
#include "GVector.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GEventAtom
 *
 * @brief Abstract interface for the event atom class.
 *
 * An event atom is a single event occuring in an instrument. It has two
 * generic attributes (m_time and m_energy) that can be access via the
 * energy() and time() methods (both methods return constant pointers).
 * Furthermore, the counts() method returns the number of atoms in the event,
 * which by definition is 1.
 ***************************************************************************/
class GEventAtom : public GEvent {

    // Friend classes
    friend class GEvents;

public:
    // Constructors and destructors
    GEventAtom(void);
    GEventAtom(const GEventAtom& atom);
    virtual ~GEventAtom(void);

    // Operators
    virtual GEventAtom& operator= (const GEventAtom& atom);

    // Pure virtual methods
    virtual double           model(GModels& models, GVector* gradient = NULL) const = 0;
    virtual const GInstDir*  dir(void) const = 0;
    virtual const GPointing* pnt(void) const = 0;
    virtual const GResponse* rsp(void) const = 0;

    // Implemented methods
    bool           isatom(void) const { return true; }
    bool           isbin(void) const { return false; }
    double         counts(void) const { return 1.0; }
    double         size(void) const { return 1.0; }
    const GEnergy* energy(void) const { return &m_energy; }
    const GTime*   time(void) const { return &m_time; }
    
protected:
    // Protected methods
    void                init_members(void);
    void                copy_members(const GEventAtom& atom);
    void                free_members(void);
    virtual GEventAtom* clone(void) const = 0;

    // Protected data area
    GTime      m_time;         //!< Event time
	GEnergy    m_energy;       //!< Event energy

private:
};

#endif /* GEVENTATOM_HPP */
