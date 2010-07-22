/***************************************************************************
 *        GApplication.i - GammaLib application base class SWIG file       *
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
 * @file GApplication.i
 * @brief GammaLib application SWIG file.
 * @author Jurgen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GApplication.hpp"
%}
//%include stl.i


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application interface defintion.
 ***************************************************************************/
class GApplication {
public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name, const std::string& version);
    GApplication(const std::string& name, const std::string& version,
                 int argc, char* argv[]);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Methods
    std::string name(void) const;
    std::string version(void) const;
    double      telapse(void) const;
    GPar*       par(const std::string& name);
};


/***********************************************************************//**
 * @brief GApplication class SWIG extension
 ***************************************************************************/
%extend GApplication {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GApplication copy() {
        return (*self);
    }
}
