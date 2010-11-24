/***************************************************************************
 *                    GPar.cpp - Application parameter                     *
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
 * @file GPar.cpp
 * @brief Application parameter class implementation
 * @author Jurgen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_LIBREADLINE
#include <stdio.h>             //!< Needed for declaration of FILE in readline
#include <readline/readline.h>
#endif
#include "GPar.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_BOOLEAN                                           "GPar::boolean()"
#define G_INTEGER                                           "GPar::integer()"
#define G_REAL                                                 "GPar::real()"
#define G_CHECK_TYPE                          "GPar::check_type(std::string)"
#define G_CHECK_MODE                          "GPar::check_mode(std::string)"
#define G_CHECK_VALUE_BOOL              "GPar::check_value_bool(std::string)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GPar::GPar(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 ***************************************************************************/
GPar::GPar(const std::string& name, const std::string& type,
           const std::string& mode, const std::string& value,
           const std::string& min, const std::string& max, 
           const std::string& prompt)
{
    // Initialise private members for clean destruction
    init_members();

    // Set parameter attributes
    m_name = name;
    this->mode(mode);
    this->type(type);
    m_min    = min;
    m_max    = max;
    this->value(value);
    m_prompt = prompt;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Object from which the instance should be built.
 ***************************************************************************/
GPar::GPar(const GPar& par)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPar::~GPar(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] par Object which should be assigned.
 ***************************************************************************/
GPar& GPar::operator= (const GPar& par)
{
    // Execute only if object is not identical
    if (this != &par) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(par);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set parameter type
 *
 * @param[in] type Parameter type.
 *
 * Valid parameter types are: b,i,r,s,f,fr,fw,fe,fn.
 ***************************************************************************/
void GPar::type(const std::string& type)
{
    // Verify that type is valid
    check_type(type);

    // Set mode
    m_type  = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter mode
 *
 * @param[in] mode Parameter mode.
 *
 * Valid parameter modes are: a,h,l,q,hl,ql,lh,lq.
 ***************************************************************************/
void GPar::mode(const std::string& mode)
{
    // Verify that mode is valid
    check_mode(mode);

    // Set mode
    m_mode  = mode;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 ***************************************************************************/
void GPar::value(const std::string& value)
{
    // Verify that valie is valid
    check_value(value);

    // Set mode
    m_value  = value;

    // Signal update
    m_update = true;

    // Don't query parameter again
    if (m_mode == "q")  m_mode = "h";
    if (m_mode == "ql") m_mode = "hl";
    if (m_mode == "lq") m_mode = "lh";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns parameter value (as string)
 ***************************************************************************/
std::string GPar::value(void)
{
    // Query parameter
    query();

    // Return value
    return m_value;
}


/***********************************************************************//**
 * @brief Returns boolean
 *
 * @exception GException::par_error
 *            Parameter is not of Boolean type.
 *
 * Returns the value of a Boolean parameter.
 ***************************************************************************/
bool GPar::boolean(void)
{
    // Query parameter
    query();

    // Check if parameter is boolean
    if (m_type != "b")
        throw GException::par_error(G_BOOLEAN, 
              "Boolean parameter type expected, \""+m_type+"\" found.");

    // Set result
    bool result = (toupper(m_value) == "YES" || toupper(m_value) == "TRUE");

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns integer
 *
 * @exception GException::par_error
 *            Parameter is not of integer type.
 *
 * Returns the value of an integer parameter.
 ***************************************************************************/
int GPar::integer(void)
{
    // Query parameter
    query();

    // Check if parameter is integer
    if (m_type != "i")
        throw GException::par_error(G_INTEGER, 
              "Integer parameter type expected, \""+m_type+"\" found.");

    // Set result
    int result = toint(m_value);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns real
 *
 * @exception GException::par_error
 *            Parameter is not of real type.
 *
 * Returns the value of a real parameter.
 ***************************************************************************/
double GPar::real(void)
{
    // Query parameter
    query();

    // Check if parameter is integer
    if (m_type != "r")
        throw GException::par_error(G_REAL, 
              "Real parameter type expected, \""+m_type+"\" found.");

    // Set result
    double result = todouble(m_value);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns if parameter mode is 'learn'
 ***************************************************************************/
bool GPar::islearn(void) const
{
    // Assign result
    bool result = (m_mode == "hl" || m_mode == "ql" || m_mode == "lh" || 
                   m_mode == "lq");

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns if parameter mode is 'query'
 ***************************************************************************/
bool GPar::isquery(void) const
{
    // Assign result
    bool result = (m_mode == "q" || m_mode == "ql" || m_mode == "lq");

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GPar::init_members(void)
{
    // Initialise members
    m_update = false;
    m_name.clear();
    m_type.clear();
    m_mode.clear();
    m_value.clear();
    m_min.clear();
    m_max.clear();
    m_prompt.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par Object from which members which should be copied.
 ***************************************************************************/
void GPar::copy_members(const GPar& par)
{
    // Copy attributes
    m_update = par.m_update;
    m_name   = par.m_name;
    m_type   = par.m_type;
    m_mode   = par.m_mode;
    m_value  = par.m_value;
    m_min    = par.m_min;
    m_max    = par.m_max;
    m_prompt = par.m_prompt;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPar::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of type string
 *
 * @param[in] type Type string.
 *
 * @exception GException::par_error
 *            Invalid parameter type specified.
 *
 * The parameter type has to be one of b,i,r,s,f,fr,fw,fe,fn. 
 * The fr,fw,fe,fn types test for read access, write access, file existence,
 * and file absence, respectively.
 ***************************************************************************/
void GPar::check_type(const std::string& type) const
{
    // Check if type is valid
    if (type != "b"  && type != "i"  && type != "r"  && type != "s" &&
        type != "f"  && type != "fr" && type != "fw" && type != "fe" &&
        type != "fn")
        throw GException::par_error(G_CHECK_TYPE, 
                                    "invalid parameter type '"+type+"'");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of mode string
 *
 * @param[in] mode Mode string.
 *
 * @exception GException::par_error
 *            Invalid parameter mode.
 *
 * The parameter mode has to be one of a,h,l,q,hl,ql,lh,lq.
 ***************************************************************************/
void GPar::check_mode(const std::string& mode) const
{
    // Check of mode is valid
    if (mode != "a"  && mode != "h"  && mode != "q" && mode != "hl" &&
        mode != "ql" && mode != "lh" && mode != "lq")
        throw GException::par_error(G_CHECK_MODE, 
                                    "invalid parameter mode '"+mode+"'");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of value string
 *
 * @param[in] value Value string.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 ***************************************************************************/
void GPar::check_value(const std::string& value) const
{
    // Make type dependend check
    if (m_type == "b")
        check_value_bool(value);
    else if (m_type == "i")
        check_value_int(value);
    else if (m_type == "r")
        check_value_real(value);
    else if (m_type == "s")
        check_value_string(value);
    else if (m_type == "f"  || m_type == "fr" || m_type == "fw" ||
             m_type == "fe" || m_type == "fn")
        check_value_filename(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of boolean value string
 *
 * @param[in] value Value string.
 *
 * @exception GException::par_error
 *            Boolean value string is not valid.
 *
 * The Boolean value string has to be one of y,n,yes,no (case insensitive).
 ***************************************************************************/
void GPar::check_value_bool(const std::string& value) const
{
    // Turn value to lower case for testing
    std::string lvalue = tolower(value);

    // Test for validity
    if (lvalue != "y" && lvalue != "yes" && lvalue != "n" && lvalue != "no")
        throw GException::par_error(G_CHECK_VALUE_BOOL,
                        "invalid Boolean value '"+value+"'; use y/n/yes/no");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of integer value string
 *
 * @param[in] value Value string.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 *
 * @todo Implement method.
 ***************************************************************************/
void GPar::check_value_int(const std::string& value) const
{
    // Convert value to integer
    int ivalue = toint(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of real value string
 *
 * @param[in] value Value string.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 *
 * @todo Implement method.
 ***************************************************************************/
void GPar::check_value_real(const std::string& value) const
{

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of string value string
 *
 * @param[in] value Value string.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 *
 * @todo Implement method.
 ***************************************************************************/
void GPar::check_value_string(const std::string& value) const
{

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of filename value string
 *
 * @param[in] value Value string.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 *
 * @todo Implement method.
 * @todo NONE is equivalent to an empty string.
 ***************************************************************************/
void GPar::check_value_filename(const std::string& value) const
{

    // Return
    return;
}


/***********************************************************************//**
 * @brief Query parameter if required
 *
 * This method queries the parameter from the stardard input if it is needed
 * to be input by the user.
 *
 * @todo Optional use of readline library
 ***************************************************************************/
void GPar::query(void)
{
    // Continue only if parameter has query mode
    if (isquery()) {

        // Dump prompt string
        std::string prompt = m_prompt;
        if (m_min.length() > 0 && m_max.length() > 0)
            prompt += " ("+m_min+"-"+m_max+")";
        else if (m_min.length() > 0)
            prompt += " ("+m_min+")";
        else if (m_max.length() > 0)
            prompt += " ("+m_max+")";
        prompt += " ["+m_value+"] ";

        // Get value
        #ifdef HAVE_LIBREADLINE
        std::string value;
        char* line = readline(prompt.c_str());
        if (line != NULL) {
            value = std::string(line);
            delete line;
        }
        #else
        std::cout << prompt;
        char line[1000];
        std::cin.getline(line, 1000);
        std::string value = std::string(line);
        #endif

        // Update value if value is not the default
        if (value.length() > 0) {
            m_value  = value;
            m_update = true;
        }

        // Don't query parameter again
        if (m_mode == "q")  m_mode = "h";
        if (m_mode == "ql") m_mode = "hl";
        if (m_mode == "lq") m_mode = "lh";

    } // endif: parameter had query mode

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] par Object to be dumped
 ***************************************************************************/
std::ostream& operator<<(std::ostream& os, const GPar& par)
{
    // Put object in stream
    os << par.name() << "=" << par.m_value;
    if (par.min().length() > 0 && par.max().length() > 0)
        os << " (" << par.min() << "-" << par.max() << ")";
    else if (par.min().length() > 0)
        os << " (>" << par.min() << ")";
    else if (par.max().length() > 0)
        os << " (<" << par.max() << ")";
    os << " t=" << par.type() << " m=" << par.mode();
    os << " " << par.prompt();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Write parameter into logger
 *
 * @param[in] log Logger.
 * @param[in] par Parameter to be written.
 ***************************************************************************/
GLog& operator<<(GLog& log, const GPar& par)
{
    // Write parameter in logger
    log << par.name() << "=" << par.m_value;
    if (par.min().length() > 0 && par.max().length() > 0)
        log << " (" << par.min() << "-" << par.max() << ")";
    else if (par.min().length() > 0)
        log << " (>" << par.min() << ")";
    else if (par.max().length() > 0)
        log << " (<" << par.max() << ")";
    log << " t=" << par.type() << " m=" << par.mode();
    log << " " << par.prompt();

    // Return logger
    return log;
}

