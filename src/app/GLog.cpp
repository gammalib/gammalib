/***************************************************************************
 *                       GLog.cpp - Information logger                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GLog.cpp
 * @brief Information logger class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <ctime>
#include <cstdarg>      // for std::va_list type
#include <cstdio>       // std::fopen, std::fgets, std::fclose, etc...
#include <cstring>      // std::memcpy
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
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
GLog::GLog(const std::string& filename, const bool& clobber)
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
    std::va_list vl;
    va_start(vl, msgFormat);
    std::vsprintf(buffer, msgFormat, vl);
    va_end(vl);

    // Append message to string buffer
    append(std::string(buffer)+"\n");

    // Return
    return;

}


/***********************************************************************//**
 * @brief Insert logger into logger
 *
 * @param[in] log Logger to be inserted.
 *
 * Insert a logger object into another logger object. This method handles
 * information that has either been written into a file, or that is still
 * present in the logger buffer.
 ***************************************************************************/
GLog& GLog::operator<<(GLog& log)
{
    // If log object has an open file then append the log buffer to this
    // buffer
    if (log.m_file == NULL) {
        m_buffer.append(log.m_buffer);
    }

    // ... otherwise read all information from file and append it to the
    // buffer
    else {

        // Force flush so that the log object writes the entire buffer
        // into the file
        log.flush(true);

        // Get the filename
        std::string filename = log.filename();

        // Close log's file
        log.close();

        // Open log's file in read mode
        std::ifstream file(filename.c_str(), std::ios::in);
        if (file) {

            // Append all lines from file to this buffer
            std::string line;
            while(std::getline(file,line)){
                m_buffer.append(line+"\n");
            }

            // Close log's file
            file.close();
        }

        // Open log's file in append mode
        log.open(filename,false);
    }

    // Flush
    flush();

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert C++ string into logger
 *
 * @param[in] str C++ string to be inserted.
 *
 * This method insert a C++ string into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const std::string& str)
{
    // Add string to buffer
    append(str);

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert C string into logger
 *
 * @param[in] str C string to be inserted.
 *
 * This method insert a C string into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const char* str)
{
    // Add string to buffer
    append(std::string(str));

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert character value into logger
 *
 * @param[in] value Character value to be inserted.
 *
 * This method inserts a character value into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const char& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    append(oss.str());

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert unsigned character value into logger
 *
 * @param[in] value Unsigned character value to be inserted.
 *
 * This method inserts an unsigned character value into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const unsigned char& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    append(oss.str());

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert bool into logger
 *
 * @param[in] value Boolean to be inserted.
 *
 * Depending on the argument, this method inserts either the character string
 * "true" or "false" into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const bool& value)
{
    // Add boolean to buffer
    if (value) {
        append("true");
    }
    else {
        append("false");
    }

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert integer into logger
 *
 * @param[in] value Integer to be inserted.
 *
 * This method inserts an integer into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const int& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    append(oss.str());

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert unsigned integer into logger
 *
 * @param[in] value Unsigned integer to be inserted.
 *
 * This method inserts an unsigned integer into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const unsigned int& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    append(oss.str());

    // Return logger
    return *this;
}


/***********************************************************************//**
 * @brief Insert double precision value into logger
 *
 * @param[in] value Double precision value to be inserted.
 *
 * This method inserts a double precision value into the logger.
 ***************************************************************************/
GLog& GLog::operator<<(const double& value)
{
    // Allocate output stream
    std::ostringstream oss;

    // Put value into stream
    oss << value;

    // Append stream to buffer
    append(oss.str());

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
    append("\n");

    // Return logger
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
void GLog::open(const std::string& filename, const bool& clobber)
{
    // Store the filename
    m_filename = filename;

    // Close any existing file
    close();

    // Open file
    if (clobber) {
        m_file = std::fopen(filename.c_str(), "w");
    }
    else {
        m_file = std::fopen(filename.c_str(), "a");
    }

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
        std::fclose(m_file);
        m_file = NULL;
    }

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
    m_file       = NULL;
    m_filename.clear();
    m_name.clear();
    m_buffer.clear();
    m_chatter    = NORMAL;

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
    m_file       = log.m_file;
    m_filename   = log.m_filename;
    m_name       = log.m_name;
    m_buffer     = log.m_buffer;
    m_chatter    = log.m_chatter;

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
 * buffer in lines that a separated by a '\n' character. The following
 * streams are currently implemented (and will be filled in parallel):
 * stdout, stderr, and an ASCII file.
 ***************************************************************************/
void GLog::flush(const bool& force)
{
    // Check if buffer should be flushed
    bool flush = (force || size() > m_max_length);

    // Flush buffer
    if (flush) {

        // Flush buffer until it is empty
        while (m_buffer.size() > 0) {

            // Find next CR
            std::string line;
            std::size_t pos = m_buffer.find_first_of("\n", 0);
            if (pos == std::string::npos) {
                line    = m_buffer;
                m_buffer.clear();
            }
            else {
                line     = m_buffer.substr(0, pos) + "\n";
                m_buffer = m_buffer.substr(pos+1, m_buffer.size()-pos);
            }

            // Put line into requested streams
            if (m_stdout) {
                std::cout << line;
            }
            if (m_stderr) {
                std::cerr << line;
            }
            if (m_file != NULL) {
                std::fprintf(m_file, "%s", line.c_str());
            }

        } // endwhile: flush until empty

    } // endif: flush was required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write string as header into logger
 *
 * @param[in] arg Header string.
 * @param[in] level Header level (0,1,2,3).
 *
 * Writes a string as header into the logger. Various header levels exist:
 * 0: 80 character wide centred header
 * 1: Header framed by '=' characters
 * 2: Header framed by '-' characters
 * 3: Header enclosed by '===' string
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
    case 3:
        text  = "=== " + arg + " ===";
        break;
    default:
        break;
    }

    // Write header into logger
    if (frame.length() > 0) {
        append(frame);
        append("\n");
    }
    if (text.length() > 0) {
        append(text);
        append("\n");
    }
    if (frame.length() > 0) {
        append(frame);
        append("\n");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return current date
 *
 * Returns the current date as string in the format yyyy-mm-ddThh:mm:ss.
 ***************************************************************************/
std::string GLog::strdate(void) const
{
    // Allocate variables
    struct std::tm timeStruct;
    std::time_t    now;
    char           buffer[100];

    // Get time
    now = std::time(NULL);
    #ifdef HAVE_GMTIME_R   
    std::gmtime_r(&now, &timeStruct);
    #else
    std::memcpy(&timeStruct, gmtime(&now), sizeof(struct tm));
    #endif

    // Write message type, time and task name to buffer
    std::sprintf(buffer, "%04d-%02d-%02dT%02d:%02d:%02d",
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


/***********************************************************************//**
 * @brief Return prefix
 *
 * Returns the prefix for each line.
 ***************************************************************************/
std::string GLog::prefix(void) const
{
    // Initialize prefix to empty string
    std::string prefix = "";

    // Add date if requested
    if (m_use_date){
        prefix.append(strdate());
    }

    // Add name if requested
    if (m_name.length() > 0){
        prefix.append(" "+m_name);
    }

    // If there is a prefix then add separator
    if (prefix.length() > 0){
        prefix.append(": ");
    }

    // Add any indent
    prefix.append(fill(" ", m_indent));

    // Return prefix
    return prefix;
}


/***********************************************************************//**
 * @brief Append string to the buffer
 *
 * @param[in] arg String to append
 *
 * This method appends a string to the buffer and prepends, if required, the
 * current date and name at the beginning of each line.
 ***************************************************************************/
void GLog::append(std::string arg)
{
    // If the buffer is empty or if the last charater is a \n, prepend a
    // prefix at the beginning of the string to be inserted.
    if (m_buffer.size() == 0 || m_buffer[m_buffer.size()-1] == '\n') {

        // Prepend prefix
        arg.insert(0, prefix());
    }

    // Search the first CR (\n)
    std::size_t pos = arg.find_first_of("\n",0);

    // Search all \n characters. Ignore the last CR.
    while (pos != std::string::npos && pos < arg.size()-1) {

        // Prepend prefix
        std::string pre = prefix();
        arg.insert(pos+1, pre);

        // Search next CR
        pos = arg.find_first_of("\n",pos+1+pre.size());

    } // endwhile

    // Add string to buffer
    m_buffer.append(arg);

    // Flush Buffer
    flush();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
