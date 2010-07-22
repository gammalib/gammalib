/***************************************************************************
 *                 GPars.i - Application parameters SWIG file              *
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
 * @file GPars.i
 * @brief Application parameters SWIG file.
 * @author Jurgen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPars.hpp"
%}
//%include stl.i


/***********************************************************************//**
 * @class GPars
 *
 * @brief Application parameters interface defintion.
 ***************************************************************************/
class GPars {

public:
    // Constructors and destructors
    GPars(void);
    GPars(const std::string& filename);
    GPars(const std::string& filename, const std::vector<std::string>& args);
    GPars(const GPars& pars);
    ~GPars(void);
 
    // Methods
    void  load(const std::string& filename);
    void  load(const std::string& filename, const std::vector<std::string>& args);
    void  save(const std::string& filename);
    GPar* par(const std::string& name);
};


/***********************************************************************//**
 * @brief GPars class SWIG extension
 ***************************************************************************/
%extend GPars {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GPars copy() {
        return (*self);
    }
}
