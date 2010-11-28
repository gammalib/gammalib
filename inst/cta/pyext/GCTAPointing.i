/***************************************************************************
 *          GCTAPointing.i  -  CTA pointing class python bindings          *
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
 * @file GCTAPointing.i
 * @brief GCTAPointing class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAPointing.hpp"
%}
//%include stl.i


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief Interface for the CTA pointing class
 ***************************************************************************/
class GCTAPointing : public GPointing {
public:
    // Constructors and destructors
    GCTAPointing(void);
    GCTAPointing(const GCTAPointing& pnt);
    ~GCTAPointing(void);

    // Methods
    void          clear(void);
    GCTAPointing* clone(void) const;
};
