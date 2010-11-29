/***************************************************************************
 *              GRoi.i  -  Region of interest class python I/F             *
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
 * @file GRoi.i
 * @brief GRoi class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GRoi.hpp"
%}


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
public:
    // Constructors and destructors
    GRoi(void);
    GRoi(const GRoi& roi);
    virtual ~GRoi(void);

    // Pure virtual methods
    virtual GRoi* clone(void) const = 0;
};


/***********************************************************************//**
 * @brief GRoi class extension
 ***************************************************************************/
%extend GRoi {
    /*
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    */
};
