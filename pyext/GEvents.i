/***************************************************************************
 *              GEvents.i  -  Event container class python I/F             *
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
 * @file GEvents.i
 * @brief GEvents class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEvents.hpp"
%}


/***********************************************************************//**
 * @class GEvents
 *
 * @brief GEvents container class interface defintion.
 *
 * This class is an abstract container base class for events. Events are
 * generally associated to an observation, and the class keeps the pointer
 * to an existing observations as a member.
 ***************************************************************************/
class GEvents {
public:
    // Constructors and destructors
    GEvents();
    GEvents(const GEvents& events);
    virtual ~GEvents();

    // Pure virtual methods
    virtual void     clear(void) = 0;
    virtual GEvents* clone(void) const = 0;
    virtual void     load(const std::string& filename) = 0;
    virtual GEvent*  pointer(int index) = 0;
    virtual int      number(void) const = 0;
    virtual int      size(void) const = 0;
    virtual bool     islist(void) const = 0;
    virtual bool     iscube(void) const = 0;
};


/***********************************************************************//**
 * @brief GEvents class extension
 ***************************************************************************/
%extend GEvents {
    /*
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    */
};
