/***************************************************************************
 * GModelSpatialPtsrc.i  -  Spatial point source model class SWIG interface*
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
 * @file GModelSpatialPtsrc.i
 * @brief GModelSpatialPtsrc class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialPtsrc.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialPtsrc
 *
 * @brief Point source model interface definition.
 ***************************************************************************/
class GModelSpatialPtsrc  : public GModelSpatial {
public:
    // Constructors and destructors
    explicit GModelSpatialPtsrc(void);
    explicit GModelSpatialPtsrc(const GSkyDir& dir);
    GModelSpatialPtsrc(const GModelSpatialPtsrc& model);
    ~GModelSpatialPtsrc(void);

    // Methods
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
    double     eval(const GSkyDir& srcDir);
    double     eval_gradients(const GSkyDir& srcDir);
    bool       isptsource(void) const { return true; }
    double     ra(void) const { return m_ra.real_value(); }
    double     dec(void) const { return m_dec.real_value(); }
};


/***********************************************************************//**
 * @brief GModelSpatialPtsrc class extension
 ***************************************************************************/
%extend GModelSpatialPtsrc {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    GModelSpatialPtsrc copy() {
        return (*self);
    }
};
