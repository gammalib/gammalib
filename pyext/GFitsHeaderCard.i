/***************************************************************************
 *         GFitsHeaderCard.hpp  - FITS header card class SWIG file         *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHeaderCard.i
 * @brief GFitsHeaderCard class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHeaderCard.hpp"
%}

%include stl.i


/***********************************************************************//**
 * @class GFitsHeaderCard
 *
 * @brief Implements FITS header card class SWIG interface
 *
 * This class implements FITS header card. Each card consists of a
 * keyname (string), a value (string, floating pointer, integer or logical)
 * and a comment (string). COMMENT or HISTORY cards do not have any value.
 ***************************************************************************/
class GFitsHeaderCard {
public:
    // Constructors & Destructors
    GFitsHeaderCard(void);
    GFitsHeaderCard(const std::string& keyname, const std::string& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const double& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const int& value,
                    const std::string& comment);
    GFitsHeaderCard(const GFitsHeaderCard& card);
    virtual ~GFitsHeaderCard(void);

    // Methods to set card properties
    void         keyname(const std::string& keyname);
    void         value(const std::string& value);
    void         value(const bool& value);
    void         value(const float& value);
    void         value(const double& value);
    void         value(const unsigned short& value);
    void         value(const short& value);
    void         value(const unsigned int& value);
    void         value(const int& value);
    void         value(const long& value);
    void         value(const unsigned long& value);
    void         value(const long long& value);
    void         unit(const std::string& unit);
    void         comment(const std::string& comment);

    // Methods to get card properties
    std::string  keyname(void) const;
    std::string  value(void) const;
    int          decimals(void) const;
    std::string  unit(void) const;
    std::string  comment(void) const;
    std::string  string(void);
    double       real(void);
    int          integer(void);
};


/***********************************************************************//**
 * @brief GFitsHeaderCard class SWIG extension
 ***************************************************************************/
%extend GFitsHeaderCard {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GFitsHeaderCard copy() {
        return (*self);
    }
}
