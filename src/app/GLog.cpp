/***************************************************************************
 *                       GLog.cpp - Information logger                     *
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
 * @file GLog.hpp
 * @brief Information logger class implementation
 * @author Jurgen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdarg.h>     // for "va_list" type
#include "GLog.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

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
GLog::GLog(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Log file constructor
 *
 * @param[in] filename Name of log file.
 * @param[in] clobber If true overwrite any existing file, otherwise append.
 *
 * Construct a logger with an open log file.
 ***************************************************************************/
GLog::GLog(const std::string& filename, bool clobber)
{
    // Initialise private members for clean destruction
    init_members();

    // Open file
    open(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] log Object from which the instance should be built.
 ***************************************************************************/
GLog::GLog(const GLog& log)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(log);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLog::~GLog(void)
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
 * @param[in] log Object which should be assigned.
 ***************************************************************************/
GLog& GLog::operator= (const GLog& log)
{ 
    // Execute only if object is not identical
    if (this != &log) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(log);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/**************************************************************************//**
 * @brief Message logging operator
 *
 * @param[in] msgFormat Pointer to message format.
 * @param[in] ... Message to put into log.
 *
 * This operator provides message logging capabilities in C style. It expects
 * a C-style format string followed by an arbitrary number of arguments that
 * is used to fill the format.
 *
 * @todo A static buffer of length 8192 is allocated. Although this is likely
 * sufficient for any normal application there is nothing that prevents in
 * principle the overflow of the buffer. A more sophisticated scheme should be
 * implemented.
 ******************************************************************************/
void GLog::operator()(const char *msgFormat, ...)
{
    // Allocate message buffer
    char buffer[8192];

    // Put message into buffer
    va_list vl;
    va_start(vl, msgFormat);
    vsprintf(buffer, msgFormat, vl);
    va_end(vl);

    // Append message to string buffer
    m_buffer.append(buffer);
    m_buffer.append("\n");

    // Flush buffer
    flush(true);

    // Return
    return;

}


/***********************************************************************//**
 * @brief Insert logger into logger
 *
 * @param[in] log Logger to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const GLog& log)
{
    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert C++ string into logger
 *
 * @param[in] str C++ string to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const std::string& str)
{
    // Add string to buffer
    m_buffer.append(str);

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert C string into logger
 *
 * @param[in] str C string to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const char* str)
{
    // Add string to buffer
    m_buffer.append(str);

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert character value into logger
 *
 * @param[in] value Character value to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const char& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    m_buffer.append(oss.str());

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert unsigned character value into logger
 *
 * @param[in] value Unsigned character value to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const unsigned char& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    m_buffer.append(oss.str());

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert bool into logger
 *
 * @param[in] value Boolean to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const bool& value)
{
    // Add boolean to buffer
    if (value)
        m_buffer.append("true");
    else
        m_buffer.append("false");

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert integer into logger
 *
 * @param[in] value Integer to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const int& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    m_buffer.append(oss.str());

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert unsigned integer into logger
 *
 * @param[in] value Unsigned integer to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const unsigned int& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    m_buffer.append(oss.str());

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert double precision value into logger
 *
 * @param[in] value Double precision value to be inserted.
 ***************************************************************************/
GLog& GLog::operator<<(const double& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    m_buffer.append(oss.str());

    // Flush buffer
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Ostream function output operator
 *
 * @param[in] fn Ostream function.
 *
 * The function output operator is interpreted as a line feed operator.
 *
 * @todo This is a quick an dirty implementation of the std::endl operator.
 * A clean GLog::endl operator should be implemented instead, yet this needs
 * some deeper understanding of our the 
 ***************************************************************************/
GLog& GLog::operator<<(std::ostream& (*fn)(std::ostream&))
{
    // Append CR to buffer
    m_buffer.append("\n");

    // Flush buffer
    flush(true);

    // Return output stream
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 ***************************************************************************/
void GLog::clear(void)
{
    // Free members
    free_members();

    // Init members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return the number of characters that are actually in the buffer
 ***************************************************************************/
int GLog::size(void) const
{
    // Return
    return m_buffer.size();
}


/***********************************************************************//**
 * @brief Open log file
 *
 * @param[in] filename Name of log file.
 * @param[in] clobber If true overwrite any existing file, otherwise append.
 *
 * Opens a file for logging. If a log file was already open it is closed
 * before opening a new file.
 ***************************************************************************/
void GLog::open(const std::string& filename, bool clobber)
{
    // Close any existing file
    close();

    // Open file
    if (clobber)
        m_file = fopen(filename.c_str(), "w");
    else
        m_file = fopen(filename.c_str(), "a");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close log file
 *
 * Close log file. This method even works if no file has been opened. In any
 * case flushing of the buffer is enforced.
 ***************************************************************************/
void GLog::close(void)
{
    // Flush buffer
    flush(true);

    // Close any open file
    if (m_file != NULL) {
        fclose(m_file);
        m_file = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set data flag that controls date prefixing
 *
 * @param[in] flag Enable/disable date flag (true/false). 
 ***************************************************************************/
void GLog::date(bool flag)
{
    // Set flag
    m_use_date = flag;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Enables/disables logging into standard output stream
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * Enables or disables logging into the standard output stream.
 ***************************************************************************/
void GLog::cout(bool flag)
{
    // Set flag
    m_stdout = flag;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Enables/disables logging into standard error stream
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * Enables or disables logging into the standard error stream.
 ***************************************************************************/
void GLog::cerr(bool flag)
{
    // Set flag
    m_stderr = flag;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set name to put into prefix
 *
 * @param[in] name Name to put into prefix. 
 ***************************************************************************/
void GLog::name(const std::string& name)
{
    // Set name
    m_name = name;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set the maximum buffer size
 *
 * @param[in] size Maximum buffer size.
 *
 * Set the maximum number of characters allowed in the internal buffer before
 * the buffer is flushed into the specified streams or file.
 *
 * @todo No maximum buffer size checking is performed. This could be
 * implemented but is not really crucial since the physical buffer size can
 * indeed be larger than this maximum size.
 ***************************************************************************/
void GLog::max_size(int size)
{
    // Set size
    m_max_length = size;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set indentation
 *
 * @param[in] indent Indentation in number of characters.
 *
 * Set the indentation of the logger. This specifies the number of whitespace
 * characters that is added between the prefix and the log message of each
 * line. Indentation allows for easy formatting of text in the log file.
 ***************************************************************************/
void GLog::indent(int indent)
{
    // Set size
    m_indent = indent;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLog::init_members(void)
{
    // Initialise members
    m_max_length = 8192;
    m_indent     = 0;
    m_stdout     = false;
    m_stderr     = false;
    m_use_date   = false;
    m_newline    = true;
    m_file       = NULL;
    m_name.clear();
    m_buffer.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] log Object from which members which should be copied.
 ***************************************************************************/
void GLog::copy_members(const GLog& log)
{
    // Copy attributes
    m_max_length = log.m_max_length;
    m_indent     = log.m_indent;
    m_stdout     = log.m_stdout;
    m_stderr     = log.m_stderr;
    m_use_date   = log.m_use_date;
    m_newline    = log.m_newline;
    m_name       = log.m_name;
    m_buffer     = log.m_buffer;
    m_file       = log.m_file;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * Deletes class member after flushing the buffer and closing the log file.
 ***************************************************************************/
void GLog::free_members(void)
{
    // Flush buffer
    flush(true);

    // Close log file
    close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Flush string buffer into log file
 *
 * @param[in] force If true, force flushing.
 *
 * Flush string buffer if it is full or if flushing is enforced. This method
 * writes the output string in the relevant streams. It decomposes the
 * buffer in lines that a separated by a '\n' character and preprends, if
 * requested, the current date and name. The following streams are currently
 * implemented (and will be filled in parallel): stdout, stderr, and an ASCII
 * file.
 ***************************************************************************/
void GLog::flush(bool force)
{
    // Check if buffer should be flushed
    bool flush = (force || size() > m_max_length);

    // Flush buffer
    if (flush) {

        // Flush buffer until it is empty
        while (m_buffer.size() > 0) {

            // Initialise newline for next line
            bool newline = m_newline;

            // Find next CR
            std::string line;
            size_t pos = m_buffer.find_first_of("\n", 0);
            if (pos == std::string::npos) {
                line    = m_buffer;
                newline = false;
                m_buffer.clear();
            }
            else {
                line     = m_buffer.substr(0, pos) + "\n";
                m_buffer = m_buffer.substr(pos+1, m_buffer.size()-pos);
                newline  = true;
            }

            // Prepend prefix if we are at the beginning of a new line
            if (m_newline) {
                std::string prefix = "";
                if (m_use_date)
                    prefix.append(strdate());
                if (m_name.length() > 0)
                    prefix.append(" "+m_name);
                if (prefix.length() > 0)
                    prefix.append(": ");
                line = prefix + fill(" ", m_indent) + line;
            }

            // Put line into requested streams
            if (m_stdout)
                std::cout << line;
            if (m_stderr)
                std::cerr << line;
            if (m_file != NULL)
                fprintf(m_file, line.c_str());

            // Set newline flag for next line
            m_newline = newline;

        } // endwhile: flush until empty

    } // endif: flush was required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write string as header into logger
 *
 * @param[in] arg Header string.
 * @param[in] level Header level (0,1,...).
 *
 * Writes a string as header into the logger. Various header levels exist:
 * 0: 80 character wide centred header
 * 1: Header framed by '=' characters
 * 2: Header framed by '-' characters
 ***************************************************************************/
void GLog::header(const std::string& arg, int level)
{
    // Declare frame and text strings
    std::string frame;
    std::string text;

    // Create level dependent header strings
    switch (level) {
    case 0:
        break;
    case 1:
    case 2:
        text  = "| " + arg + " |";
        frame = (level == 1) ? "+" +  fill("=", text.length()-2) + "+"
                             : "+" +  fill("-", text.length()-2) + "+";
        break;
    default:
        break;
    }

    // Write header into logger
    if (text.length() > 0) {
        m_buffer.append(frame);
        m_buffer.append("\n");
        m_buffer.append(text);
        m_buffer.append("\n");
        m_buffer.append(frame);
        m_buffer.append("\n");
        flush(true);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return current date
 *
 * Returns the current date as string in the format yyyy-mm-ddThh:mm:ss.
 ***************************************************************************/
std::string GLog::strdate(void)
{
    // Allocate variables
    struct tm timeStruct;
    time_t    now;
    char      buffer[100];

    // Get time
    now = time(NULL);
    #ifdef HAVE_GMTIME_R   
    gmtime_r(&now, &timeStruct);
    #else
    memcpy(&timeStruct, gmtime(&now), sizeof(struct tm));
    #endif

    // Write message type, time and task name to buffer
    sprintf(buffer, "%04d-%02d-%02dT%02d:%02d:%02d",
                    timeStruct.tm_year + 1900,
                    timeStruct.tm_mon + 1,
                    timeStruct.tm_mday,
                    timeStruct.tm_hour,
                    timeStruct.tm_min,
                    timeStruct.tm_sec);

    // Build string from buffer
    std::string date = buffer;

    // Return date
    return date;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
