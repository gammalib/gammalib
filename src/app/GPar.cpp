/***************************************************************************
 *                    GPar.cpp - Application parameter                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GPar.cpp
 * @brief Application parameter class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_LIBREADLINE
#include <cstdio>             //!< Needed for declaration of FILE in readline
#include <readline/readline.h>
#endif
#include "GPar.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_STRING_SET                             "GPar::string(std::string&)"
#define G_FILENAME_SET                         "GPar::filename(std::string&)"
#define G_BOOLEAN_SET                                  "GPar::boolean(bool&)"
#define G_INTEGER_SET                                   "GPar::integer(int&)"
#define G_REAL_SET                                      "GPar::real(double&)"
#define G_STRING_GET                                         "GPar::string()"
#define G_FILENAME_GET                                     "GPar::filename()"
#define G_BOOLEAN_GET                                       "GPar::boolean()"
#define G_INTEGER_GET                                       "GPar::integer()"
#define G_REAL_GET                                             "GPar::real()"
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
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] name Parameter name.
 * @param[in] type Parameter type.
 * @param[in] mode Parameter mode.
 * @param[in] value Parameter value.
 * @param[in] min Parameter minimum.
 * @param[in] max Parameter maximum.
 * @param[in] prompt Parameter prompt string.
 *
 * Constructs a parameter from parameter attributes.
 ***************************************************************************/
GPar::GPar(const std::string& name, const std::string& type,
           const std::string& mode, const std::string& value,
           const std::string& min, const std::string& max, 
           const std::string& prompt)
{
    // Initialise members
    init_members();

    // Set parameter attributes
    m_name   = name;
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
 * @param[in] par Parameter.
 ***************************************************************************/
GPar::GPar(const GPar& par)
{
    // Initialise members
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
 * @param[in] par Parameter.
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
 * @brief Clone instance
 ***************************************************************************/
GPar* GPar::clone(void) const
{
    return (new GPar(*this));
}


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

    // Set type
    m_type = type;

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
    m_mode = mode;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 *
 * This method sets the parameter value dependent on the parameter type.
 * It is a parser method that verifies that the given parameter value is
 * compatible with the parameter type.
 ***************************************************************************/
void GPar::value(const std::string& value)
{
    // Verify that value is valid
    check_value(value);

    // Set value
    set_value(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set string parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::par_error
 *            Parameter is not of string type.
 *
 * This method sets a string parameter. The method only applies to parameters
 * of type "s". Other parameter types will produce an exception.
 ***************************************************************************/
void GPar::string(const std::string& value)
{
    // Check if parameter is a string parameter
    if (m_type != "s")
        throw GException::par_error(G_STRING_SET, name(),
              "attempted to set "+par_type_string(m_type)+
              " parameter with string value.");

    // Verify that value is a valid string
    check_value_string(value);

    // Set value
    set_value(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set filename parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::par_error
 *            Parameter is not of filename type.
 *
 * This method sets a filename parameter. The method only applies to filename
 * parameters. Other parameter types will produce an exception.
 ***************************************************************************/
void GPar::filename(const std::string& value)
{
    // Check if parameter is a filename parameter
    if (!isfilename()) {
        throw GException::par_error(G_FILENAME_SET, name(),
              "attempted to set "+par_type_string(m_type)+
              " parameter with filename.");
    }

    // Verify that value is a valid filename
    check_value_filename(value);

    // Set value
    set_value(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set bool parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::par_error
 *            Parameter is not of boolean type.
 *
 * This method sets a boolean parameter. The method only applies to parameters
 * of type "b". Other parameter types will produce an exception.
 ***************************************************************************/
void GPar::boolean(const bool& value)
{
    // Check if parameter is boolean
    if (m_type != "b") {
        throw GException::par_error(G_BOOLEAN_SET, name(),
              "attempted to set "+par_type_string(m_type)+
              " parameter with boolean value.");
    }

    // Set value string
    std::string value_string = (value) ? "yes" : "no";

    // Set value
    set_value(value_string);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set integer parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::par_error
 *            Parameter is not of integer type.
 *
 * This method sets an integer parameter. The method only applies to
 * parameters of type "i". Other parameter types will produce an exception.
 ***************************************************************************/
void GPar::integer(const int& value)
{
    // Check if parameter is integer
    if (m_type != "i") {
        throw GException::par_error(G_INTEGER_SET, name(),
              "attempted to set "+par_type_string(m_type)+
              " parameter with integer value.");
    }

    // Set value string
    std::string value_string = str(value);

    // Set value
    set_value(value_string);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set real parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::par_error
 *            Parameter is not of real type.
 *
 * This method sets a real parameter. The method only applies to parameters
 * of type "r". Other parameter types will produce an exception.
 ***************************************************************************/
void GPar::real(const double& value)
{
    // Check if parameter is boolean
    if (m_type != "r") {
        throw GException::par_error(G_REAL_SET, name(),
              "attempted to set "+par_type_string(m_type)+
              " parameter with real value.");
    }

    // Set value string
    std::string value_string = str(value);

    // Set value
    set_value(value_string);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns parameter value as string
 *
 * This method queries any parameter and returns it as a string. No parameter
 * type checking is performed.
 ***************************************************************************/
std::string GPar::value(void)
{
    // Query parameter
    query();

    // Return value
    return m_value;
}


/***********************************************************************//**
 * @brief Returns string parameter value
 *
 * @exception GException::par_error
 *            Parameter is not of string type.
 *
 * This method queries and returns a string parameter. The method only
 * applies to parameters of type "s". Other parameter types will produce an
 * exception.
 ***************************************************************************/
std::string GPar::string(void)
{
    // Check if parameter is a string parameter
    if (m_type != "s") {
        throw GException::par_error(G_STRING_GET, name(),
              "attempted reading "+par_type_string(m_type)+
              " parameter as a string value.");
    }

    // Query parameter
    query();

    // Return value
    return m_value;
}


/***********************************************************************//**
 * @brief Returns filename parameter value
 *
 * @exception GException::par_error
 *            Parameter is not of filename type.
 *
 * This method queries and returns a filename parameter. The method only
 * applies to filename parameters. Other parameter types will produce an
 * exception. Any environment variables that are encountered within the
 * filename are expanded automatically.
 ***************************************************************************/
std::string GPar::filename(void)
{
    // Check if parameter is a filename parameter
    if (!isfilename()) {
        throw GException::par_error(G_FILENAME_GET, name(),
              "attempted reading "+par_type_string(m_type)+
              " parameter as a filename.");
    }

    // Query parameter
    query();
    
    // Return value
    return (expand_env(m_value));
}


/***********************************************************************//**
 * @brief Returns boolean
 *
 * @exception GException::par_error
 *            Parameter is not of boolean type.
 *
 * This method queries and returns a boolean parameter. The method only
 * applies to parameters of type "b". Other parameter types will produce an
 * exception.
 ***************************************************************************/
bool GPar::boolean(void)
{
    // Check if parameter is boolean
    if (m_type != "b") {
        throw GException::par_error(G_BOOLEAN_GET, name(),
              "attempted reading "+par_type_string(m_type)+
              " parameter as a boolean value.");
    }

    // Query parameter
    query();

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
 * This method queries and returns an integer parameter. The method only
 * applies to parameters of type "i". Other parameter types will produce an
 * exception.
 ***************************************************************************/
int GPar::integer(void)
{
    // Check if parameter is integer
    if (m_type != "i") {
        throw GException::par_error(G_INTEGER_GET, name(),
              "attempted reading "+par_type_string(m_type)+
              " parameter as an integer value.");
    }

    // Query parameter
    query();

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
 * This method queries and returns a real parameter. The method only
 * applies to parameters of type "r". Other parameter types will produce an
 * exception.
 ***************************************************************************/
double GPar::real(void)
{
    // Check if parameter is integer
    if (m_type != "r") {
        throw GException::par_error(G_REAL_GET, name(),
              "attempted reading "+par_type_string(m_type)+
              " parameter as a real value.");
    }

    // Query parameter
    query();

    // Set result
    double result = todouble(m_value);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Signals if parameter mode is "learn"
 *
 * A parameter is in mode learn when it has one of the following modes:
 * hl,ql,lh, or lq.
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
 * @brief Signals if parameter mode is "query"
 *
 * A parameter will be queried when it has one of the following modes:
 * q, ql, or lq.
 ***************************************************************************/
bool GPar::isquery(void) const
{
    // Assign result
    bool result = (m_mode == "q" || m_mode == "ql" || m_mode == "lq");

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Signals if parameter mode is of type "filename"
 *
 * A parameter is of type "filename" if it has one of the following types:
 * f, fr, fw, fe, or fn.
 ***************************************************************************/
bool GPar::isfilename(void) const
{
    // Assign result
    bool result = (m_type == "f"  || m_type == "fr" || m_type == "fw" ||
                   m_type == "fe" || m_type == "fn");

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print parameter
 *
 * @todo Implement writing of parameter options. I guess that parameter
 *       options occur if there is only a minimum but not a maximum
 *       parameter. Check IRAF documentation.
 ***************************************************************************/
std::string GPar::print(void) const
{
    // Write parameter name
    std::string result = parformat(name());

    // Write value (use m_value here to avoid querying when printing)
    result.append(m_value);

    // Write limits
    if (min().length() > 0 && max().length() > 0) {
        result.append(" ("+min()+"-"+max()+")");
    }
    else if (min().length() > 0) {
        result.append(" (>"+min()+")");
    }
    else if (max().length() > 0) {
        result.append(" (<"+max()+")");
    }

    // Write type information
    result.append(" [t="+type()+", m="+mode()+"]");

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
 *            Invalid parameter type.
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
        type != "fn") {
        throw GException::par_error(G_CHECK_TYPE, name(),
              "invalid parameter type \""+type+"\"");
    }

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
        mode != "ql" && mode != "lh" && mode != "lq") {
        throw GException::par_error(G_CHECK_MODE, name(),
              "invalid parameter mode \""+mode+"\"");
    }

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
    // Make type dependent check
    if (m_type == "b") {
        check_value_bool(value);
    }
    else if (m_type == "i") {
        check_value_int(value);
    }
    else if (m_type == "r") {
        check_value_real(value);
    }
    else if (m_type == "s") {
        check_value_string(value);
    }
    else if (isfilename()) {
        check_value_filename(value);
    }

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
 * The Boolean value string has to be one of y,n,yes,no,true,false
 * (case insensitive).
 ***************************************************************************/
void GPar::check_value_bool(const std::string& value) const
{
    // Turn value to lower case for testing
    std::string lvalue = tolower(value);

    // Test for validity
    if (lvalue != "y" && lvalue != "yes" && lvalue != "true" &&
        lvalue != "n" && lvalue != "no"  && lvalue != "false") {
        throw GException::par_error(G_CHECK_VALUE_BOOL, name(),
              "invalid Boolean value \""+value+"\"; use y/n/yes/no/true/false");
    }

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
    //int ivalue = toint(value);

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
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 *
 * Set parameter value, signal it for update and disable parameter querying.
 *
 * @todo Implement parameter verification (min, max, options)
 ***************************************************************************/
void GPar::set_value(const std::string& value)
{
    // Set value
    m_value = value;

    // Signal update
    m_update = true;

    // Don't query parameter again
    stop_query();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Query parameter if required
 *
 * This method queries the parameter from the stardard input if it is needed
 * to be input by the user.
 *
 * @todo Add environment variable substitution.
 ***************************************************************************/
void GPar::query(void)
{
    // Continue only if parameter has query mode
    if (isquery()) {

        // Dump prompt string
        std::string prompt = m_prompt;
        if (m_min.length() > 0 && m_max.length() > 0) {
            prompt += " ("+m_min+"-"+m_max+")";
        }
        else if (m_min.length() > 0) {
            prompt += " ("+m_min+")";
        }
        else if (m_max.length() > 0) {
            prompt += " ("+m_max+")";
        }
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

        // Strip any whitespace from string
        value = strip_whitespace(value);

        // Update value if value is not the default
        if (value.length() > 0) {
            m_value  = value;
            m_update = true;
        }

        // Don't query parameter again
        stop_query();

    } // endif: parameter had query mode

    // Return
    return;
}


/***********************************************************************//**
 * @brief Don't query parameter
 ***************************************************************************/
void GPar::stop_query(void)
{
    // Don't query parameter again
    if (m_mode == "q")  m_mode = "h";
    if (m_mode == "ql") m_mode = "hl";
    if (m_mode == "lq") m_mode = "lh";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return type string
 *
 * @param[in] type Parameter type.
 *
 * Translates parameter type character(s) into a human readable type string.
 * Valid parameter types are: b,i,r,s,f,fr,fw,fe,fn. If the parameter type
 * is not valid, the method returns "unknown".
 ***************************************************************************/
std::string GPar::par_type_string(const std::string& type) const
{
    // Initialize empty type string
    std::string type_string = "";

    // Set type dependent string
    if (type == "b") {
        type_string.append("boolean");
    }
    else if (type == "i") {
        type_string.append("integer");
    }
    else if (type == "r") {
        type_string.append("real");
    }
    else if (type == "s") {
        type_string.append("string");
    }
    else if (type == "f"  || type == "fr" || type == "fw" ||
             type == "fe" || type == "fn") {
        type_string.append("filename");
    }
    else {
        type_string.append("unknown");
    }

    // Return type string
    return type_string;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream.
 * @param[in] par Parameter.
 ***************************************************************************/
std::ostream& operator<<(std::ostream& os, const GPar& par)
{
    // Write parameter in output stream
    os << par.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Write parameter into logger
 *
 * @param[in] log Logger.
 * @param[in] par Parameter.
 ***************************************************************************/
GLog& operator<<(GLog& log, const GPar& par)
{
    // Write parameter in logger stream
    log << par.print();

    // Return logger
    return log;
}
