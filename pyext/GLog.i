/***************************************************************************
 *                        GLog.i - Information logger                      *
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
 * @file GLog.i
 * @brief Information logger python bindings
 * @author Jurgen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLog.hpp"
%}


/***********************************************************************//**
 * @class GLog
 *
 * @brief Information logger interface defintion.
 ***************************************************************************/
class GLog {
public:
    // Constructors and destructors
    GLog(void);
    GLog(const std::string& filename, bool clobber = false);
    GLog(const GLog& log);
    ~GLog(void);

    // Operators
    void  operator()(const char *msgFormat, ...);
    GLog& operator<<(const GLog& log);
    GLog& operator<<(const std::string& str);
    GLog& operator<<(const char* str);
    GLog& operator<<(const char& value);
    GLog& operator<<(const unsigned char& value);
    GLog& operator<<(const bool& value);
    GLog& operator<<(const int& value);
    GLog& operator<<(const unsigned int& value);
    GLog& operator<<(const double& value);
    GLog& operator<<(std::ostream& (*fn)(std::ostream&));

    // Methods
    void        clear(void);
    int         size(void) const;
    void        open(const std::string& filename, bool clobber = false);
    void        close(void);
    void        date(bool flag);
    void        cout(bool flag);
    void        cerr(bool flag);
    void        name(const std::string& name);
    void        max_size(int size);
    void        indent(int indent);
    void        header0(const std::string& arg) { header(arg, 0); }
    void        header1(const std::string& arg) { header(arg, 1); }
    void        header2(const std::string& arg) { header(arg, 2); }
    bool        date(void) const { return m_use_date; }
    bool        cout(void) const { return m_stdout; }
    bool        cerr(void) const { return m_stderr; }
    std::string name(void) const { return m_name; }
    int         max_size(void) const { return m_max_length; }
    int         indent(void) const { return m_indent; }
};


/***********************************************************************//**
 * @brief GLog class SWIG extension
 ***************************************************************************/
%extend GLog {
    GLog copy() {
        return (*self);
    }
}
