/***************************************************************************
 *                 GLATLtCube.i  -  Fermi LAT lifetime cube                *
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
 * @file GLATLtCube.i
 * @brief GLATLtCube class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATLtCube.hpp"
%}


/***********************************************************************//**
 * @class GLATLtCube
 *
 * @brief Interface for the Fermi LAT lifetime cube.
 ***************************************************************************/
class GLATLtCube {

public:
    // Constructors and destructors
    GLATLtCube(void);
    GLATLtCube(const std::string& filename);
    GLATLtCube(const GLATLtCube& cube);
    virtual ~GLATLtCube(void);

    // Operators
    double operator() (const GSkyDir& dir, const GEnergy& energy, _ltcube_ctheta fct);
    double operator() (const GSkyDir& dir, const GEnergy& energy, _ltcube_ctheta_phi fct);

    // Methods
    void        clear(void);
    GLATLtCube* clone(void) const;
    void        load(const std::string& filename);
    void        save(const std::string& filename, bool clobber=false) const;
};


/***********************************************************************//**
 * @brief GLATLtCube class extension
 ***************************************************************************/
%extend GLATLtCube {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GLATLtCube copy() {
        return (*self);
    }
};
