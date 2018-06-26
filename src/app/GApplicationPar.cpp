/***************************************************************************
 *               GApplicationPar.cpp - Application parameter               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GApplicationPar.cpp
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
#include <climits>            //!< Needed for declaration of LONG_MAX
#include "GApplicationPar.hpp"
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GTime.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_STRING_SET                  "GApplicationPar::string(std::string&)"
#define G_FILENAME_SET                "GApplicationPar::filename(GFilename&)"
#define G_TIME_SET                            "GApplicationPar::time(GTime&)"
#define G_BOOLEAN_SET                       "GApplicationPar::boolean(bool&)"
#define G_INTEGER_SET                        "GApplicationPar::integer(int&)"
#define G_REAL_SET                           "GApplicationPar::real(double&)"
#define G_STRING_GET                              "GApplicationPar::string()"
#define G_FILENAME_GET                          "GApplicationPar::filename()"
#define G_TIME_GET                                  "GApplicationPar::time()"
#define G_BOOLEAN_GET                            "GApplicationPar::boolean()"
#define G_INTEGER_GET                            "GApplicationPar::integer()"
#define G_REAL_GET                                  "GApplicationPar::real()"
#define G_CHECK_TYPE               "GApplicationPar::check_type(std::string)"
#define G_CHECK_MODE               "GApplicationPar::check_mode(std::string)"
#define G_CHECK_VALUE_BOOL   "GApplicationPar::check_value_bool(std::string)"
#define G_CHECK_VALUE_INT    "GApplicationPar::check_value_int(std::string&)"
#define G_CHECK_VALUE_REAL  "GApplicationPar::check_value_real(std::string&)"
#define G_CHECK_OPTIONS        "GApplicationPar::check_options(std::string&)"

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
GApplicationPar::GApplicationPar(void)
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
GApplicationPar::GApplicationPar(const std::string& name,
                                 const std::string& type,
                                 const std::string& mode,
                                 const std::string& value,
                                 const std::string& min,
                                 const std::string& max, 
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

    // Since the call to the value method stops the query, we reset the
    // queried flag so that the parameter will get queried.
    m_queried = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Parameter.
 ***************************************************************************/
GApplicationPar::GApplicationPar(const GApplicationPar& par)
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
GApplicationPar::~GApplicationPar(void)
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
 * @return Parameter.
 ***************************************************************************/
GApplicationPar& GApplicationPar::operator=(const GApplicationPar& par)
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
 * @brief Clear parameter
 ***************************************************************************/
void GApplicationPar::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone parameter
 *
 * @return Pointer to deep copy of parameter.
 ***************************************************************************/
GApplicationPar* GApplicationPar::clone(void) const
{
    return (new GApplicationPar(*this));
}


/***********************************************************************//**
 * @brief Set parameter type
 *
 * @param[in] type Parameter type.
 *
 * Sets the parameter type. Valid parameter types are:
 * b,i,r,s,f,fr,fw,fe,fn,t.
 ***************************************************************************/
void GApplicationPar::type(const std::string& type)
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
 * Sets the parameter model. Valid parameter modes are: a,h,l,q,hl,ql,lh,lq.
 ***************************************************************************/
void GApplicationPar::mode(const std::string& mode)
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
 * Sets the parameter value dependent on the parameter type.
 *
 * It is a parser method that verifies that the given parameter value is
 * compatible with the parameter type.
 ***************************************************************************/
void GApplicationPar::value(const std::string& value)
{
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
 * @exception GException::invalid_value
 *            Parameter is not of string type.
 *
 * This method sets a string parameter. The method only applies to parameters
 * of type "s". Other parameter types will produce an exception.
 ***************************************************************************/
void GApplicationPar::string(const std::string& value)
{
    // Check if parameter is a string parameter
    if (m_type != "s") {
        std::string msg = "Attempt to set "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" with string value.";
        throw GException::invalid_value(G_STRING_SET, msg);
    }

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
 * @exception GException::invalid_value
 *            Parameter is not of filename type.
 *
 * This method sets a filename parameter. The method only applies to filename
 * parameters. Other parameter types will produce an exception.
 ***************************************************************************/
void GApplicationPar::filename(const GFilename& value)
{
    // Check if parameter is a filename parameter
    if (!is_filename()) {
        std::string msg = "Attempt to set "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" with filename value.";
        throw GException::invalid_value(G_FILENAME_SET, msg);
    }

    // Set value
    set_value(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::invalid_value
 *            Parameter is not of time type.
 *
 * This method sets a time parameter with a UTC string. The method only
 * applies to time parameters. Other parameter types will produce an
 * exception.
 ***************************************************************************/
void GApplicationPar::time(const GTime& value)
{
    // Check if parameter is a filename parameter
    if (m_type != "t") {
        std::string msg = "Attempt to set "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" with filename value.";
        throw GException::invalid_value(G_TIME_SET, msg);
    }

    // Set value
    set_value(value.utc());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set bool parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::invalid_value
 *            Parameter is not of boolean type.
 *
 * This method sets a boolean parameter. The method only applies to parameters
 * of type "b". Other parameter types will produce an exception.
 ***************************************************************************/
void GApplicationPar::boolean(const bool& value)
{
    // Check if parameter is boolean
    if (m_type != "b") {
        std::string msg = "Attempt to set "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" with boolean value.";
        throw GException::invalid_value(G_BOOLEAN_SET, msg);
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
 * @exception GException::invalid_value
 *            Parameter is not of integer type.
 *
 * This method sets an integer parameter. The method only applies to
 * parameters of type "i". Other parameter types will produce an exception.
 ***************************************************************************/
void GApplicationPar::integer(const int& value)
{
    // Check if parameter is integer
    if (m_type != "i") {
        std::string msg = "Attempt to set "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" with integer value.";
        throw GException::invalid_value(G_INTEGER_SET, msg);
    }

    // Set value string
    std::string value_string = gammalib::str(value);

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
 * @exception GException::invalid_value
 *            Parameter is not of real type.
 *
 * This method sets a real parameter. The method only applies to parameters
 * of type "r". Other parameter types will produce an exception.
 ***************************************************************************/
void GApplicationPar::real(const double& value)
{
    // Check if parameter is boolean
    if (m_type != "r") {
        std::string msg = "Attempt to set "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" with real value.";
        throw GException::invalid_value(G_REAL_SET, msg);
    }

    // Set value string
    std::string value_string = gammalib::str(value);

    // Set value
    set_value(value_string);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Query parameter if required
 *
 * This method queries the parameter from the stardard input if it is needed
 * to be input by the user.
 ***************************************************************************/
void GApplicationPar::query(void)
{
    // Continue only if parameter has query mode and it was not yet queried
    if (is_query() && !was_queried()) {

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
        value = gammalib::strip_whitespace(value);

        // Update value if value is not the default
        if (value.length() > 0) {
            set_value(value);
            m_update = true;
        }

        // Don't query parameter again
        stop_query();

    } // endif: parameter had query mode

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns parameter value as string
 *
 * This method queries any parameter and returns it as a string. No parameter
 * type checking is performed.
 ***************************************************************************/
std::string GApplicationPar::value(void)
{
    // Query parameter
    query();

    // Return value
    return m_value;
}


/***********************************************************************//**
 * @brief Returns string parameter value
 *
 * @exception GException::invalid_value
 *            Parameter is not of string type.
 *
 * This method queries and returns a string parameter. The method only
 * applies to parameters of type "s". Other parameter types will produce an
 * exception.
 ***************************************************************************/
std::string GApplicationPar::string(void)
{
    // Check if parameter is a string parameter
    if (m_type != "s") {
        std::string msg = "Attempt to read "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" as a string value.";
        throw GException::invalid_value(G_STRING_GET, msg);
    }

    // Query parameter
    query();

    // Check if parameter is valid
    if (m_status != ST_VALID) {
        std::string msg = "Parameter \""+m_name+"\" is "+
                          par_status_string()+". Please specify a valid"
                          " parameter value.";
        throw GException::invalid_value(G_STRING_GET, msg);
    }

    // Return value
    return m_value;
}


/***********************************************************************//**
 * @brief Returns filename parameter value
 *
 * @exception GException::invalid_value
 *            Parameter is not of filename type.
 *
 * This method queries and returns a filename parameter. The method only
 * applies to filename parameters. Other parameter types will produce an
 * exception. Any environment variables that are encountered within the
 * filename are expanded automatically.
 ***************************************************************************/
GFilename GApplicationPar::filename(void)
{
    // Check if parameter is a filename parameter
    if (!is_filename()) {
        std::string msg = "Attempt to read "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" as a filename value.";
        throw GException::invalid_value(G_FILENAME_GET, msg);
    }

    // Query parameter
    query();
    
    // Check if parameter is valid
    if (m_status != ST_VALID) {
        std::string msg = "Parameter \""+m_name+"\" is "+
                          par_status_string()+". Please specify a valid"
                          " parameter value.";
        throw GException::invalid_value(G_FILENAME_GET, msg);
    }

    // Return value
    return (GFilename(m_value));
}


/***********************************************************************//**
 * @brief Returns time parameter value
 *
 * @param[in] ref Time reference system.
 * @return Time.
 *
 * @exception GException::invalid_value
 *            Parameter is not of time type.
 *
 * This method queries and returns a time parameter. The method only applies
 * to time parameters. Other parameter types will produce an exception. If
 * the time is specified as "Mission Elapsed Time" (MET) the specified
 * time reference system will be used for conversion.
 ***************************************************************************/
GTime GApplicationPar::time(const GTimeReference& ref)
{
    // Check if parameter is a time parameter
    if (m_type != "t") {
        std::string msg = "Attempt to read "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" as a time value.";
        throw GException::invalid_value(G_TIME_GET, msg);
    }

    // Query parameter
    query();
    
    // Check if parameter is valid
    if (m_status != ST_VALID) {
        std::string msg = "Parameter \""+m_name+"\" is "+
                          par_status_string()+". Please specify a valid"
                          " parameter value.";
        throw GException::invalid_value(G_TIME_GET, msg);
    }

    // Return value
    return (GTime(m_value, ref));
}


/***********************************************************************//**
 * @brief Returns boolean
 *
 * @exception GException::invalid_value
 *            Parameter is not of boolean type.
 *
 * This method queries and returns a boolean parameter. The method only
 * applies to parameters of type "b". Other parameter types will produce an
 * exception.
 ***************************************************************************/
bool GApplicationPar::boolean(void)
{
    // Check if parameter is boolean
    if (m_type != "b") {
        std::string msg = "Attempt to read "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" as a boolean value.";
        throw GException::invalid_value(G_BOOLEAN_GET, msg);
    }

    // Query parameter
    query();

    // Check if parameter is valid
    if (m_status != ST_VALID) {
        std::string msg = "Parameter \""+m_name+"\" is "+
                          par_status_string()+". Please specify a valid"
                          " parameter value.";
        throw GException::invalid_value(G_BOOLEAN_GET, msg);
    }

    // Convert boolean value to upper case
    std::string uvalue = gammalib::toupper(m_value);

    // Set result
    bool result = (uvalue == "YES"  || uvalue == "Y" ||
                   uvalue == "TRUE" || uvalue == "T");

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns integer
 *
 * @exception GException::invalid_value
 *            Parameter is not of integer type.
 *
 * This method queries and returns an integer parameter. The method only
 * applies to parameters of type "i". Other parameter types will produce an
 * exception.
 ***************************************************************************/
int GApplicationPar::integer(void)
{
    // Check if parameter is integer
    if (m_type != "i") {
        std::string msg = "Attempt to read "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" as a integer value.";
        throw GException::invalid_value(G_INTEGER_GET, msg);
    }

    // Query parameter
    query();

    // Check if parameter is valid
    if (m_status != ST_VALID) {
        std::string msg = "Parameter \""+m_name+"\" is "+
                          par_status_string()+". Please specify a valid"
                          " parameter value.";
        throw GException::invalid_value(G_INTEGER_GET, msg);
    }

    // Set result
    int result = gammalib::toint(m_value);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns real
 *
 * @exception GException::invalid_value
 *            Parameter is not of real type.
 *
 * This method queries and returns a real parameter. The method only
 * applies to parameters of type "r". Other parameter types will produce an
 * exception.
 ***************************************************************************/
double GApplicationPar::real(void)
{
    // Check if parameter is integer
    if (m_type != "r") {
        std::string msg = "Attempt to read "+par_type_string(m_type)+
                          " parameter \""+m_name+"\" as a real value.";
        throw GException::invalid_value(G_REAL_GET, msg);
    }

    // Query parameter
    query();

    // Check if parameter is valid
    if (m_status != ST_VALID) {
        std::string msg = "Parameter \""+m_name+"\" is "+
                          par_status_string()+". Please specify a valid"
                          " parameter value.";
        throw GException::invalid_value(G_REAL_GET, msg);
    }

    // Set result
    double result = gammalib::todouble(m_value);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Signals if parameter mode is "learn"
 *
 * A parameter is in mode learn when it has one of the following modes:
 * hl,ql,lh, or lq.
 ***************************************************************************/
bool GApplicationPar::is_learn(void) const
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
bool GApplicationPar::is_query(void) const
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
bool GApplicationPar::is_filename(void) const
{
    // Assign result
    bool result = (m_type == "f"  || m_type == "fr" || m_type == "fw" ||
                   m_type == "fe" || m_type == "fn");

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Signals if parameter is valid
 ***************************************************************************/
bool GApplicationPar::is_valid(void)
{
    // Query parameter
    query();

    // Return validity
    return (m_status == ST_VALID);
}


/***********************************************************************//**
 * @brief Signals if parameter is undefined
 ***************************************************************************/
bool GApplicationPar::is_undefined(void)
{
    // Query parameter
    query();

    // Return validity
    return (m_status == ST_UNDEFINED);
}


/***********************************************************************//**
 * @brief Signals if parameter is not a number
 ***************************************************************************/
bool GApplicationPar::is_notanumber(void)
{
    // Query parameter
    query();

    // Return validity
    return (m_status == ST_NAN);
}


/***********************************************************************//**
 * @brief Set class from pickled string vector
 *
 * @param[in] string String vector containing class information.
 ***************************************************************************/
void GApplicationPar::pickle(const std::vector<std::string>& string)
{
    // Clear object
    clear();

    // Extract members
    m_update  = (bool)gammalib::toint(string[0]);
    m_queried = (bool)gammalib::toint(string[1]);
    m_name    = string[2];
    m_type    = string[3];
    m_mode    = string[4];
    m_value   = string[5];
    m_min     = string[6];
    m_max     = string[7];
    m_prompt  = string[8];
    m_status  = (Status)gammalib::toint(string[9]);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pickled string vector
 *
 * @return String vector containing class information.
 ***************************************************************************/
std::vector<std::string> GApplicationPar::pickle(void) const
{
    // Allocate vector of strings with 10 elements
    std::vector<std::string> string(10);

    // Set vector elements from class members
    string[0] = gammalib::str(m_update);
    string[1] = gammalib::str(m_queried);
    string[2] = m_name;
    string[3] = m_type;
    string[4] = m_mode;
    string[5] = m_value;
    string[6] = m_min;
    string[7] = m_max;
    string[8] = m_prompt;
    string[9] = gammalib::str((int)m_status);

    // Return string vector
    return (string);
}


/***********************************************************************//**
 * @brief Print parameter
 *
 * @param[in] chatter Chattiness.
 * @return String containing parameter information.
 ***************************************************************************/
std::string GApplicationPar::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Write parameter name
        result = gammalib::parformat(name());

        // Write value (use m_value here to avoid querying when printing)
        switch (m_status) {
        case ST_VALID:
            result.append(m_value);
            break;
        case ST_UNDEFINED:
            result.append("undefined");
            break;
        case ST_NAN:
            result.append("nan");
            break;
        case ST_UNDERFLOW:
            result.append(m_value+" (underflow)");
            break;
        case ST_OVERFLOW:
            result.append(m_value+" (overflow)");
            break;
        }

        // Write limits
        if (min().length() > 0 && max().length() > 0) {
            result.append(" ("+min()+"-"+max()+")");
        }
        else if (min().length() > 0) {
            result.append(" ("+min()+")");
        }

        // Write type information
        result.append(" [t="+type()+", m="+mode()+"]");

    } // endif: chatter was not silent

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
void GApplicationPar::init_members(void)
{
    // Initialise members
    m_update  = false;
    m_queried = false;
    m_name.clear();
    m_type.clear();
    m_mode.clear();
    m_value.clear();
    m_min.clear();
    m_max.clear();
    m_prompt.clear();
    m_status = ST_VALID;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par Parameter.
 ***************************************************************************/
void GApplicationPar::copy_members(const GApplicationPar& par)
{
    // Copy attributes
    m_update  = par.m_update;
    m_queried = par.m_queried;
    m_name    = par.m_name;
    m_type    = par.m_type;
    m_mode    = par.m_mode;
    m_value   = par.m_value;
    m_min     = par.m_min;
    m_max     = par.m_max;
    m_prompt  = par.m_prompt;
    m_status  = par.m_status;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplicationPar::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of type string
 *
 * @param[in] type Type string.
 *
 * @exception GException::invalid_value
 *            Invalid parameter type.
 *
 * The parameter type has to be one of b,i,r,s,f,fr,fw,fe,fn,t.
 * The fr,fw,fe,fn types test for read access, write access, file existence,
 * and file absence, respectively.
 ***************************************************************************/
void GApplicationPar::check_type(const std::string& type) const
{
    // Check if type is valid
    if (type != "b"  && type != "i"  && type != "r"  && type != "s" &&
        type != "f"  && type != "fr" && type != "fw" && type != "fe" &&
        type != "fn" && type != "t") {
        std::string msg = "Invalid parameter type \""+type+"\" encountered"
                          " for parameter \""+m_name+"\".";
        throw GException::invalid_value(G_CHECK_TYPE, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of mode string
 *
 * @param[in] mode Mode string.
 *
 * @exception GException::invalid_value
 *            Invalid parameter mode.
 *
 * The parameter mode has to be one of a,h,l,q,hl,ql,lh,lq.
 ***************************************************************************/
void GApplicationPar::check_mode(const std::string& mode) const
{
    // Check of mode is valid
    if (mode != "a"  && mode != "h"  && mode != "q" && mode != "hl" &&
        mode != "ql" && mode != "lh" && mode != "lq") {
        std::string msg = "Invalid parameter mode \""+mode+"\" encountered"
                          " for parameter \""+m_name+"\".";
        throw GException::invalid_value(G_CHECK_MODE, msg);
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
void GApplicationPar::check_value(const std::string& value) const
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
    else if (m_type == "t") {
        check_value_time(value);
    }
    else if (is_filename()) {
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
 * @exception GException::invalid_value
 *            Boolean value string is not valid.
 *
 * The Boolean value string has to be one of y,yes,true,t or n,no,false,f
 * (case insensitive).
 ***************************************************************************/
void GApplicationPar::check_value_bool(const std::string& value) const
{
    // Turn value to lower case for testing
    std::string lvalue = gammalib::tolower(value);

    // Test for validity
    if (lvalue != "y" && lvalue != "yes" && lvalue != "true"  && lvalue != "t" &&
        lvalue != "n" && lvalue != "no"  && lvalue != "false" && lvalue != "f") {
        std::string msg = "Invalid boolean value \""+value+"\" encountered"
                          " for parameter \""+m_name+"\". Use"
                          " y/n/yes/no/t/f/true/false";
        throw GException::invalid_value(G_CHECK_VALUE_BOOL, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of integer parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::invalid_value
 *            Integer parameter value outside validity range
 *
 * If either options or a validity range has been specified, check if the
 * integer parameter value satisfies the validity constraints.
 *
 * Requires that m_type, m_status, m_min and m_max are set. The method does
 * not verify if m_type is valid.
 ***************************************************************************/
void GApplicationPar::check_value_int(const std::string& value) const
{
    // Check only if value is not undefined
    if (m_status != ST_UNDEFINED) {

        // Check for options
        bool has_options = check_options(value);

        // If no options check has been done and if there is a m_min and
        // m_max value then perform an integer check
        if (!has_options && m_min.length() > 0 && m_max.length() > 0) {

            // Throw an exception if we have a NAN parameter
            if (m_status == ST_NAN) {
                std::string msg = "Parameter \""+m_name+"\" value \""+value+
                                  "\" outside validity range ["+m_min+","+
                                  m_max+"].";
                throw GException::invalid_value(G_CHECK_VALUE_INT, msg);
            }

            // Convert value and range to integers
            int ivalue = gammalib::toint(value);
            int imin   = gammalib::toint(m_min);
            int imax   = gammalib::toint(m_max);

            // Check if value is outside range
            if (imin > ivalue || imax < ivalue) {
                std::string msg = "Parameter \""+m_name+"\" value \""+value+
                                  "\" outside validity range ["+m_min+","+
                                  m_max+"].";
                throw GException::invalid_value(G_CHECK_VALUE_INT, msg);
            }

        } // endif: there was no options check and there was a range specified

    } // endif: value was not undefined

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of real parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::invalid_value
 *            Real parameter value outside validity range
 *
 * If either options or a validity range has been specified, check if the
 * real parameter value satisfies the validity constraints.
 *
 * Requires that m_type, m_status, m_min and m_max are set. The method does
 * not verify if m_type is valid.
 ***************************************************************************/
void GApplicationPar::check_value_real(const std::string& value) const
{
    // Check only if value is not undefined
    if (m_status != ST_UNDEFINED) {

        // Check for options
        bool has_options = check_options(value);

        // If no options check has been done and if there is a m_min and
        // m_max value then perform an integer check
        if (!has_options && m_min.length() > 0 && m_max.length() > 0) {

            // Throw an exception if we have a NAN parameter
            if (m_status == ST_NAN) {
                std::string msg = "Parameter \""+m_name+"\" value \""+value+
                                  "\" outside validity range ["+m_min+","+
                                  m_max+"].";
                throw GException::invalid_value(G_CHECK_VALUE_REAL, msg);
            }

            // Convert value and range to doubles
            double dvalue = gammalib::todouble(value);
            double dmin   = gammalib::todouble(m_min);
            double dmax   = gammalib::todouble(m_max);

            // Check if value is outside range
            if (dmin > dvalue || dmax < dvalue) {
                std::string msg = "Parameter \""+m_name+"\" value \""+value+
                                  "\" outside validity range ["+m_min+","+
                                  m_max+"].";
                throw GException::invalid_value(G_CHECK_VALUE_REAL, msg);
            }

        } // endif: there was no options check and there was a range specified

    } // endif: value was not undefined

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of string parameter value
 *
 * @param[in] value Parameter value.
 *
 * If options have been specified, check if the string parameter value
 * satisfies the options.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 ***************************************************************************/
void GApplicationPar::check_value_string(const std::string& value) const
{
    // Check for options
    check_options(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of filename parameter value
 *
 * @param[in] value Parameter value.
 *
 * If options have been specified, check if the filename parameter value
 * satisfies the options.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 *
 * @todo NONE is equivalent to an empty string.
 ***************************************************************************/
void GApplicationPar::check_value_filename(const std::string& value) const
{
    // Check for options
    check_options(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test validity of time parameter value
 *
 * @param[in] value Parameter value.
 *
 * If options have been specified, check if the time parameter value
 * satisfies the options.
 *
 * Requires that m_type, m_min and m_max are set. The method does not verify
 * if m_type is valid.
 ***************************************************************************/
void GApplicationPar::check_value_time(const std::string& value) const
{
    // Check for options
    check_options(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test if parameter value satisfies possible options
 *
 * @param[in] value Parameter value.
 * @return True if there were options, false otherwise.
 *
 * @exception GException::invalid_value
 *            Value does not satisfy one of the possible options.
 *
 * In case that the parameter has different options (m_max field not set and
 * at least one character chain in m_min field, multiple chains separated by
 * pipe symbol), check that the parameter value corresponds to one of the
 * values in the m_min field. For strings, and for strings only, the
 * parameter check is case insensitive.
 ***************************************************************************/
bool GApplicationPar::check_options(const std::string& value) const
{
    // Initialise options flags
    bool has_options = false;

    // Continue only if the m_max field is not set
    if (m_max.length() == 0) {

        // Get parameter options
        std::vector<std::string> options = gammalib::split(m_min, "|");

        // Determine number of options
        int num_options = options.size();

        // Continue only if there are options
        if (num_options > 0) {

            // Signal that there were options
            has_options = true;

            // Initalise found flag
            bool found = false;

            // Strip whitespace from all options, and for strings, convert
            // to upper case
            for (int i = 0; i < num_options; ++i) {
                options[i] = gammalib::strip_whitespace(options[i]);
                if (m_type == "s") {
                    options[i] = gammalib::toupper(options[i]);
                }
            }

            // Set parameter value
            std::string v = (m_type == "s") ? gammalib::toupper(value) : value;

            // Now check if we can find one of the options
            for (int i = 0; i < num_options; ++i) {
                if (options[i] == v) {
                    found = true;
                    break;
                }
            }

            // If no option was found then throw exception
            if (!found) {
                std::string msg = "Parameter \""+m_name+"\" value \""+value+
                                  "\" invalid. Must be one of \""+m_min+"\".";
                throw GException::invalid_value(G_CHECK_OPTIONS, msg);
            }

        } // endif: there were parameter options

    } // endif: the m_max field was not set

    // Return options flag
    return has_options;
}


/***********************************************************************//**
 * @brief Set parameter status
 *
 * @param[in] value Parameter value.
 * @return Updated parameter value.
 *
 * Set parameter status depending on the content of the value field.
 *
 * For an integer parameter, INDEF, NONE, UNDEF or UNDEFINED will result in
 * a status of "undefined". INF, INFINITY or NAN will be transformed in the
 * maximum long number.
 *
 * For a real parameter, INDEF, NONE, UNDEF or UNDEFINED will result in a
 * status of "undefined", while INF, INFINITY or NAN will result in a status
 * of "nan".
 *
 * For a time parameter, INDEF, NONE, UNDEF or UNDEFINED will result in a
 * status of "undefined".
 *
 * For a filename parameter, an empty string or INDEF, NONE, UNDEF or
 * UNDEFINED will result in a status of "undefined".
 ***************************************************************************/
std::string GApplicationPar::set_status(const std::string& value)
{
    // Initialise result value
    std::string result = value;

    // Set integer status. Catch the special values that signal that a
    // parameter is undefined. Any infinity or nan is set to the maximum
    // long value (APE standard)
    if (m_type == "i") {
        std::string lvalue = gammalib::tolower(value);
        if (lvalue == "indef" ||
            lvalue == "none"  ||
            lvalue == "undef" ||
            lvalue == "undefined") {
            m_status = ST_UNDEFINED;
        }
        else if (lvalue == "inf" ||
                 lvalue == "infinity" ||
                 lvalue == "nan") {
            result   = gammalib::str(LONG_MAX);
            m_status = ST_VALID;
        }
        else {
            m_status = ST_VALID;
        }
    }

    // Set real status. Catch the special values that signal that a
    // parameter is undefined, infinity or nan (APE standard).
    else if (m_type == "r") {
        std::string lvalue = gammalib::tolower(value);
        if (lvalue == "indef" ||
            lvalue == "none"  ||
            lvalue == "undef" ||
            lvalue == "undefined") {
            m_status = ST_UNDEFINED;
        }
        else if (lvalue == "inf" ||
                 lvalue == "infinity" ||
                 lvalue == "nan") {
            m_status = ST_NAN;
        }
        else {
            m_status = ST_VALID;
        }
    }

    // Set time status. Catch the special values that signal that a parameter
    // is undefined.
    else if (m_type == "t") {
        std::string lvalue = gammalib::tolower(value);
        if (lvalue == "indef" ||
            lvalue == "none"  ||
            lvalue == "undef" ||
            lvalue == "undefined") {
            m_status = ST_UNDEFINED;
        }
        else {
            m_status = ST_VALID;
        }
    }

    // Set filename status. Catch the special values that signal that a
    // parameter is undefined.
    else if (m_type == "f") {
        std::string lvalue = gammalib::tolower(value);
        if (lvalue == ""      ||
            lvalue == "indef" ||
            lvalue == "none"  ||
            lvalue == "undef" ||
            lvalue == "undefined") {
            m_status = ST_UNDEFINED;
        }
        else {
            m_status = ST_VALID;
        }
    }

    // Set other status
    else {
        m_status = ST_VALID;
    }

    // Return result value
    return result;
}


/***********************************************************************//**
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 *
 * Set parameter value, signal it for update and disable parameter querying.
 ***************************************************************************/
void GApplicationPar::set_value(const std::string& value)
{
    // Set parameter status
    std::string par_value = set_status(value);

    // Check parameter value
    check_value(par_value);

    // Set parameter value
    m_value = par_value;

    // Signal update
    m_update = true;

    // Don't query parameter again
    stop_query();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Don't query parameter again
 ***************************************************************************/
void GApplicationPar::stop_query(void)
{
    // Don't query parameter again
    m_queried = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return type string
 *
 * @param[in] type Parameter type.
 *
 * Translates parameter type character(s) into a human readable type string.
 * Valid parameter types are: b,i,r,s,f,fr,fw,fe,fn,t. If the parameter type
 * is not valid, the method returns "unknown".
 ***************************************************************************/
std::string GApplicationPar::par_type_string(const std::string& type) const
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
    else if (type == "t") {
        type_string.append("time");
    }
    else {
        type_string.append("unknown");
    }

    // Return type string
    return type_string;
}


/***********************************************************************//**
 * @brief Return status string
 *
 * @return Returns the parameter status in human readable form.
 ***************************************************************************/
std::string GApplicationPar::par_status_string(void) const
{
    // Allocate status string
    std::string status;

    // Set status
    switch (m_status) {
    case ST_VALID:
        status.append("valid");
        break;
    case ST_UNDEFINED:
        status.append("undefined");
        break;
    case ST_NAN:
        status.append("NaN");
        break;
    case ST_UNDERFLOW:
        status.append("underflow");
        break;
    case ST_OVERFLOW:
        status.append("overflow");
        break;
    }

    // Return status
    return status;
}
