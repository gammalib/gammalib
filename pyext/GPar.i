/***************************************************************************
 *               GPar.i - Application parameter class SWIG file            *
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
 * @file GPar.i
 * @brief Application parameter class SWIG file.
 * @author Jurgen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPar.hpp"
%}
//%include stl.i


/***********************************************************************//**
 * @class GPar
 *
 * @brief Application parameter class interface defintion.
 ***************************************************************************/
class GPar {

    // Friend classes
    friend class GPars;

public:
    // Constructors and destructors
    GPar(void);
    GPar(const std::string& name, const std::string& type,
         const std::string& mode, const std::string& value,
         const std::string& min, const std::string& max, 
         const std::string& prompt);
    GPar(const GPar& par);
    ~GPar(void);
 
    // Methods
    void        type(const std::string& type);
    void        mode(const std::string& mode);
    void        value(const std::string& value);
    std::string name(void) const { return m_name; }
    std::string type(void) const { return m_type; }
    std::string mode(void) const { return m_mode; }
    std::string value(void);
    std::string min(void) const { return m_min; }
    std::string max(void) const { return m_max; }
    std::string prompt(void) const { return m_prompt; }
    bool        islearn(void) const;
    bool        isquery(void) const;
};


/***********************************************************************//**
 * @brief GPar class SWIG extension
 ***************************************************************************/
%extend GPar {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GPar copy() {
        return (*self);
    }
}
