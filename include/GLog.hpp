/***************************************************************************
 *                       GLog.hpp - Information logger                     *
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
 * @file GLog.hpp
 * @brief Information logger class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLOG_HPP
#define GLOG_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GTypemaps.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GLog
 *
 * @brief Information logger interface defintion.
 *
 * This class implements an interface for logging of text messages to the
 * standard output and/or error streams and into an ASCII log file. It
 * provides a C++ style interface (allowing to stream messages using code
 * such as
 *
 *     log << "This is a message";
 *
 * and a C style interface that allows making use of C style text formating
 * using code such as
 *
 *     log("These are %d messages with %d items.", n, m);
 *
 * The logger interface implements an internal character buffer that is
 * flushed once a maximum size parameter is exceeded. The buffer limit can
 * be modified using the max_size() method; by default it is set to 8192
 * characters. The logger allows prepending of the current date and the task
 * name.
 ***************************************************************************/
class GLog {

public:
    // Constructors and destructors
    GLog(void);
    GLog(const std::string& filename, const bool& clobber = false);
    GLog(const GLog& log);
    virtual ~GLog(void);

    // Operators
    GLog& operator=(const GLog& log);
    void  operator()(const char *msgFormat, ...);
    GLog& operator<<(GLog& log);
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
    void               clear(void);
    int                size(void) const;
    void               open(const std::string& filename, const bool& clobber = false);
    void               close(void);
    void               flush(const bool& force = false);
    void               date(const bool& flag);
    void               cout(const bool& flag);
    void               cerr(const bool& flag);
    void               name(const std::string& name);
    void               max_size(const int& size);
    void               indent(const int& indent);
    void               chatter(const GChatter& chatter);
    void               header0(const std::string& arg);
    void               header1(const std::string& arg);
    void               header2(const std::string& arg);
    void               header3(const std::string& arg);
    const bool&        date(void) const;
    const bool&        cout(void) const;
    const bool&        cerr(void) const;
    const std::string& name(void) const;
    const int&         max_size(void) const;
    const int&         indent(void) const;
    const GChatter&    chatter(void) const;
    const std::string& filename(void) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GLog& log);
    void        free_members(void);
    void        header(const std::string& arg, int level);
    std::string strdate(void) const;
    std::string prefix(void) const;
    void        append(std::string arg);

    // Protected data members
    int         m_max_length; //!< Maximum buffer length
    int         m_indent;     //!< Indentation of text
    bool        m_stdout;     //!< Dump in standard output
    bool        m_stderr;     //!< Dump in standard error
    bool        m_use_date;   //!< Dump date in prefix
    FILE*       m_file;       //!< Log file pointer
    std::string m_filename;   //!< Log file name
    std::string m_name;       //!< Name for prefix
    std::string m_buffer;     //!< Output string buffer
    GChatter    m_chatter;    //!< Chattiness for print() method
};


/***********************************************************************//**
 * @brief Write string as centred header into logger
 *
 * @param[in] arg Header string.
 *
 * Writes a string as 80 character wide centred header.
 ***************************************************************************/
inline
void GLog::header0(const std::string& arg)
{
    header(arg, 0);
    return;
}


/***********************************************************************//**
 * @brief Write string as header framed by '=' characters into logger
 *
 * @param[in] arg Header string.
 *
 * Writes a string as header framed by '=' characters into logger.
 ***************************************************************************/
inline
void GLog::header1(const std::string& arg)
{
    header(arg, 1);
    return;
}


/***********************************************************************//**
 * @brief Write string as header framed by '-' characters into logger
 *
 * @param[in] arg Header string.
 *
 * Writes a string as header framed by '-' characters into logger.
 ***************************************************************************/
inline
void GLog::header2(const std::string& arg)
{
    header(arg, 2);
    return;
}


/***********************************************************************//**
 * @brief Write string as header enclosed by '===' into logger
 *
 * @param[in] arg Header string.
 *
 * Writes a string as header enclosed by '===' into logger.
 ***************************************************************************/
inline
void GLog::header3(const std::string& arg)
{
    header(arg, 3);
    return;
}


/***********************************************************************//**
 * @brief Set date flag that controls date prefixing
 *
 * @param[in] flag Enable/disable date flag (true/false). 
 *
 * The date flag specifies whether each new line will be prepended by the
 * actual date when the writing of the line was started.
 ***************************************************************************/
inline
void GLog::date(const bool& flag)
{
    m_use_date = flag;
    return;
}


/***********************************************************************//**
 * @brief Set standard output stream (cout) flag
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * The cout flag specifies whether the logging should go in addition to the
 * log file also into the standard output stream.
 ***************************************************************************/
inline
void GLog::cout(const bool& flag)
{
    m_stdout = flag;
    return;
}


/***********************************************************************//**
 * @brief Enables/disables logging into standard error stream
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * The cerr flag specifies whether the logging should go in addition to the
 * log file also into the standard error stream.
 ***************************************************************************/
inline
void GLog::cerr(const bool& flag)
{
    m_stderr = flag;
    return;
}


/***********************************************************************//**
 * @brief Set name to put into prefix
 *
 * @param[in] name Name to put into prefix. 
 *
 * Sets the @p name that will be prefixed to each new line. If the date flag
 * is true, the @p name will be inserted after the date at the beginning of
 * each line.
 *
 * Specifying an empty string will suppress the insertion of the name into
 * the logger.
 ***************************************************************************/
inline
void GLog::name(const std::string& name)
{
    m_name = name;
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
inline
void GLog::max_size(const int& size)
{
    m_max_length = size;
    return;
}


/***********************************************************************//**
 * @brief Set indentation
 *
 * @param[in] indent Number of indentation characters.
 *
 * Sets the number of whitespace characters that are added after the
 * prefix and before the log message content for each line.
 *
 * Indentation allows for easy formatting of text in the log file.
 ***************************************************************************/
inline
void GLog::indent(const int& indent)
{
    m_indent = indent;
    return;
}


/***********************************************************************//**
 * @brief Set chattiness
 *
 * @param[in] chatter Chattiness.
 *
 * Sets the chattiness that is applied to print() methods that are piped into
 * the logger. The following chattiness values are available:
 * - SILENT: no logging will be performed
 * - TERSE: minimum logging
 * - NORMAL: normal level of logging
 * - EXPLICIT: detailed information will be logged
 * - VERBOSE: maximum information will be logged
 ***************************************************************************/
inline
void GLog::chatter(const GChatter& chatter)
{
    m_chatter = chatter;
    return;
}


/***********************************************************************//**
 * @brief Return date flag
 *
 * @return True if date flag is set.
 *
 * The date flag specifies whether each new line will be prepended by the
 * actual date when the writing of the line was started.
 ***************************************************************************/
inline
const bool& GLog::date(void) const
{
    return m_use_date;
}


/***********************************************************************//**
 * @brief Return standard output stream (cout) flag
 *
 * @return True if cout flag is set.
 *
 * The cout flag specifies whether the logging should go in addition to the
 * log file also into the standard output stream.
 ***************************************************************************/
inline
const bool& GLog::cout(void) const
{
    return m_stdout;
}


/***********************************************************************//**
 * @brief Return standard error stream (cerr) flag
 *
 * @return True if cerr flag is set.
 *
 * The cerr flag specifies whether the logging should go in addition to the
 * log file also into the standard error stream.
 ***************************************************************************/
inline
const bool& GLog::cerr(void) const
{
    return m_stderr;
}


/***********************************************************************//**
 * @brief Return prefix name
 *
 * @return Name to put into prefix.
 *
 * Returns the name that will be prefixed to each new line. If the name is
 * empty, nothing will be prefixed.
 ***************************************************************************/
inline
const std::string& GLog::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Return the maximum buffer size
 *
 * @return Maximum buffer size.
 *
 * Returns the maximum number of characters allowed in the internal buffer
 * before the buffer is flushed into the specified streams or file.
 ***************************************************************************/
inline
const int& GLog::max_size(void) const
{
    return m_max_length;
}


/***********************************************************************//**
 * @brief Return indentation
 *
 * @return Number of indentation characters.
 *
 * Returns the number of whitespace characters that are added after the
 * prefix and before the log message content for each line.
 ***************************************************************************/
inline
const int& GLog::indent(void) const
{
    return m_indent;
}


/***********************************************************************//**
 * @brief Return chattiness
 *
 * @return Chattiness.
 *
 * Returns the chattiness that is applied to print() methods that are piped
 * into the logger. The following chattiness values are available:
 * - SILENT: no logging will be performed
 * - TERSE: minimum logging
 * - NORMAL: normal level of logging
 * - EXPLICIT: detailed information will be logged
 * - VERBOSE: maximum information will be logged
 ***************************************************************************/
inline
const GChatter& GLog::chatter(void) const
{
    return m_chatter;
}


/***********************************************************************//**
 * @brief Return log filename
 *
 * @return Log filename.
 *
 * Returns the filename of the log file.
 ***************************************************************************/
inline
const std::string& GLog::filename(void) const
{
    return m_filename;
}

#endif /* GLOG_HPP */
