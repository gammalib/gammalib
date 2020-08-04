/***************************************************************************
 *               GApplicationPars.cpp - Application parameters             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GApplicationPars.cpp
 * @brief Application parameter container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <pwd.h>           // user/passwd function
#include <fcntl.h>         // for file locking
#include <unistd.h>        // access() function
#include <sys/stat.h>      // mkdir() function
#include <cstdlib>         // std::getenv() function
#include <cstdio>          // std::fopen(), etc. functions
#include <sstream>
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GApplicationPars.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                 "GApplicationPars::operator[](std::string&)"
#define G_AT                    "GApplicationPar& GApplicationPars::at(int&)"
#define G_APPEND                 "GApplicationPars::append(GApplicationPar&)"
#define G_INSERT1          "GApplicationPar& GApplicationPars::insert(int&, "\
                                                          "GApplicationPar&)"
#define G_INSERT2  "GApplicationPar& GApplicationPars::insert(std::string&, "\
                                                          "GApplicationPar&)"
#define G_REMOVE1                            "GApplicationPars::remove(int&)"
#define G_REMOVE2                    "GApplicationPars::remove(std::string&)"
#define G_EXTEND                "GApplicationPars::extend(GApplicationPars&)"
#define G_LOAD1                          "GApplicationPars::load(GFilename&)"
#define G_LOAD2                         "GApplicationPars::load(GFilename&, "\
                                                 "std::vector<std::string>&)"
#define G_SAVE                           "GApplicationPars::save(GFilename&)"
#define G_OUTPATH                   "GApplicationPars::outpath(std::string&)"
#define G_READ                         "GApplicationPars::read(std::string&)"
#define G_WRITE                       "GApplicationPars::write(std::string&)"
#define G_PARSE                                   "GApplicationPars::parse()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_LOCK_PARFILE                //!< Enables parfile locking
//#define G_CHECK_LOCK_PARFILE          //!< Enables check of parfile locking

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GApplicationPars::GApplicationPars(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter file constructor
 *
 * @param[in] filename Parameter filename. 
 ***************************************************************************/
GApplicationPars::GApplicationPars(const GFilename& filename)
{
    // Initialise private members for clean destruction
    init_members();

    // Load parameters
    load(filename);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] filename Parameter filename. 
 * @param[in] args Command line arguments. 
 ***************************************************************************/
GApplicationPars::GApplicationPars(const GFilename&                filename,
                                   const std::vector<std::string>& args)
{
    // Initialise private members for clean destruction
    init_members();

    // Load parameters
    load(filename, args);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pars Parameter container.
 ***************************************************************************/
GApplicationPars::GApplicationPars(const GApplicationPars& pars)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(pars);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GApplicationPars::~GApplicationPars(void)
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
 * @param[in] pars Parameter container.
 * @return Parameter container.
 ***************************************************************************/
GApplicationPars& GApplicationPars::operator=(const GApplicationPars& pars)
{
    // Execute only if object is not identical
    if (this != &pars) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(pars);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Returns reference to parameter
 *
 * @param[in] name Parameter name.
 * @return Parameter.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found in container.
 ***************************************************************************/
GApplicationPar& GApplicationPars::operator[](const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name has not been found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" has not been found in "
                          "parameter file. Please specify a valid parameter "
                          "name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to parameter (const version)
 *
 * @param[in] name Parameter name.
 * @return Parameter.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found in container.
 ***************************************************************************/
const GApplicationPar& GApplicationPars::operator[](const std::string& name) const
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name has not been found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" has not been found in "
                          "parameter file. Please specify a valid parameter "
                          "name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return (m_pars[index]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear parameter container
 ***************************************************************************/
void GApplicationPars::clear(void)
{
    // Free members
    free_members();

    // Init members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone parameter container
 *
 * @return Pointer to deep copy of parameter container.
 ***************************************************************************/
GApplicationPars* GApplicationPars::clone(void) const
{
    return (new GApplicationPars(*this));
}


/***********************************************************************//**
 * @brief Returns reference to parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GApplicationPar& GApplicationPars::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }
    #endif

    // Return reference
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
const GApplicationPar& GApplicationPars::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }
    #endif

    // Return reference
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Append parameter to container
 *
 * @param[in] par Parameter.
 * @return Reference to deep copy of appended parameter.
 *
 * @exception GException::invalid_value
 *            Parameter with same name exists already in container.
 *
 * This method appends one parameter to the parameter container. The
 * parameter provided to the method can be deleted after calling this method.
 ***************************************************************************/
GApplicationPar& GApplicationPars::append(const GApplicationPar& par)
{
    // Check if a parameter with specified name does not yet exist
    int inx = get_index(par.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to append parameter with name \""+par.name()+"\" in"
            " parameter container, but a parameter with the same name exists"
            " already at index "+gammalib::str(inx)+" in the container.\n"
            "Every parameter in the parameter container needs a unique name.";
        throw GException::invalid_value(G_APPEND, msg);
    }
    
    // Append parameter
    m_pars.push_back(par);

    // Get reference of appended parameter
    GApplicationPar& parameter = m_pars[m_pars.size()-1];
    
    // Build parameter file line
    size_t      start = 0;
    size_t      stop  = 0;
    std::string line  = parline(parameter, &start, &stop);
    
    // Append parameter file line and attributes
    m_line.push_back(m_parfile.size());
    m_parfile.push_back(line);
    m_vstart.push_back(start);
    m_vstop.push_back(stop);

    // Return parameter reference
    return parameter;
}


/***********************************************************************//**
 * @brief Append standard parameters to container
 *
 * This method appends the standard parameters to the parameter container.
 * Standard parameters are: "chatter", "clobber", "debug" and "mode".
 ***************************************************************************/
void GApplicationPars::append_standard(void)
{
    // Append standard parameters
    append(GApplicationPar("chatter","i","h","2","0","4","Chattiness of output"));
	append(GApplicationPar("clobber","b","h","yes","","","Overwrite existing output files with new output files?"));
	append(GApplicationPar("debug","b","h","no","","","Debugging mode activated"));
	append(GApplicationPar("mode","s","h","ql","","","Mode of automatic parameters"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert parameter into container
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @param[in] par Parameter.
 * @return Reference to deep copy of inserted parameter.
 *
 * Inserts a parameter into the container before the parameter with the
 * specified @p index.
 ***************************************************************************/
GApplicationPar& GApplicationPars::insert(const int& index, const GApplicationPar& par)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT1, "Parameter index", index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT1, "Parameter index", index, size());
        }
    }
    #endif

    // Check if a parameter with specified name does not yet exist
    int inx = get_index(par.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert parameter with name \""+par.name()+"\" in"
            " parameter container before index "+gammalib::str(index)+
            ", but a parameter with the same name exists already at index "+
            gammalib::str(inx)+" in the container.\n"
            "Every parameter in the parameter container needs a unique name.";
        throw GException::invalid_value(G_INSERT1, msg);
    }

    // Inserts parameter
    m_pars.insert(m_pars.begin()+index, par);

    // Get reference of appended parameter
    GApplicationPar& parameter = m_pars[m_pars.size()-1];

    // Build parameter file line
    size_t      start = 0;
    size_t      stop  = 0;
    std::string line  = parline(parameter, &start, &stop);

    // Determine at which line number of the parameter file the parameter
    // should be inserted
    int line_number = m_line[index];

    // Insert parameter file line and parameter attributes
    m_parfile.insert(m_parfile.begin()+line_number, line);
    m_line.insert(m_line.begin()+index, line_number);
    m_vstart.insert(m_vstart.begin()+index, start);
    m_vstop.insert(m_vstop.begin()+index, stop);

    // Increment the line numbers for all parameters after the inserted one
    for (int i = index+1; i < size(); ++i) {
        m_line[i]++;
    }

    // Return parameter reference
    return parameter;
}


/***********************************************************************//**
 * @brief Insert parameter into container
 *
 * @param[in] name Parameter name.
 * @param[in] par Parameter.
 * @return Reference to deep copy of inserted parameter.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found in container.
 * @exception GException::invalid_value
 *            Name of parameter exists already in container.
 *
 * Inserts a parameter into the container before the parameter with the
 * specified @p name.
 ***************************************************************************/
GApplicationPar& GApplicationPars::insert(const std::string&     name,
                                          const GApplicationPar& par)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" has not been found in"
                          " parameter file.\n"
                          "Please specify a valid parameter name.";
        throw GException::invalid_argument(G_INSERT2, msg);
    }

    // Insert by index and return parameter reference
    return (insert(index, par));
}


/***********************************************************************//**
 * @brief Remove parameter from container
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Remove parameter with specified @p index from container.
 ***************************************************************************/
void GApplicationPars::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, "Parameter index", index, size());
    }
    #endif

    // Erase parameter from container
    m_pars.erase(m_pars.begin() + index);

    // Remove parameter file line and parameter attributes
    m_parfile.erase(m_parfile.begin() + m_line[index]);
    m_line.erase(m_line.begin() + index);
    m_vstart.erase(m_vstart.begin() + index);
    m_vstop.erase(m_vstop.begin() + index);

    // Decrement the line numbers for all parameters after the removed one
    for (int i = index; i < size(); ++i) {
        m_line[i]--;
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove parameter from container
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found in container.
 *
 * Remove parameter with specified @p name from container.
 ***************************************************************************/
void GApplicationPars::remove(const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" has not been found in"
                          " parameter file.\n"
                          "Please specify a valid parameter name.";
        throw GException::invalid_argument(G_REMOVE2, msg);
    }

    // Remove by index
    remove(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append parameter container
 *
 * @param[in] pars Parameter container.
 *
 * Append parameter container to the container.
 ***************************************************************************/
void GApplicationPars::extend(const GApplicationPars& pars)
{
    // Do nothing if parameter container is empty
    if (!pars.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = pars.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all parameters and append copies 
        for (int i = 0; i < num; ++i) {

            // Check if parameter name does not yet exist
            int inx = get_index(pars[i].name());
            if (inx != -1) {
                std::string msg =
                    "Attempt to append parameter with name \""+pars[i].name()+
                    "\" to parameter container, but a parameter with the same name"
                    " exists already at index "+gammalib::str(inx)+" in the"
                    " container.\n"
                    "Every parameter in the parameter container needs a unique"
                    " name.";
                throw GException::invalid_value(G_EXTEND, msg);
            }

            // Append parameter to container
            append(pars[i]);

        } // endfor: looped over all parameters

    } // endif: parameter container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Check parameter exists
 *
 * @param[in] name Parameter name.
 * @return True if parameter with specified name exists, false otherwise.
 *
 * Determines whether a parameter with the specified name exists already in
 * the parameter container.
 ***************************************************************************/
bool GApplicationPars::contains(const std::string& name) const
{
    // Get parameter index
    int inx = get_index(name);
    
    // Return test result
    return (inx != -1);
}


/***********************************************************************//**
 * @brief Load parameters
 *
 * @param[in] filename Parameter filename.
 *
 * @exception GException::par_file_not_found
 *            Parameter file not found.
 *
 * Loads all parameters from parameter file.
 *
 * If the syspfiles path is set then load the parameter file from that
 * location and update the parameters using a paramater file that is found
 * in the users "pfiles" directory.
 *
 * Otherwise, search the parameter file in the usual location (see the
 * GApplicationPars::inpath method for more information).
 ***************************************************************************/
void GApplicationPars::load(const GFilename& filename)
{
    // Reset parameters
    m_parfile.clear();

    // If a syspfiles path was set then get parfile from this path and
    // update the parameters using a copy found in the users "pfiles"
    // directory
    std::string path = syspfiles_path(filename.url());
    if (!path.empty()) {

        // Read parfile
        read(path);

        // Parse parfile
        parse();

        // Synchronize parfile using a copy found in the users "pfiles"
        // directory
        path = pfiles_path(filename.url());
        if (!path.empty()) {
            synchronise(path);
        }

    }

    // ... otherwise get path to parameter file for input
    else {

        // Get parfile file path
        path = inpath(filename.url());

        // If file path is empty then throw an exception
        if (path.empty()) {
            throw GException::par_file_not_found(G_LOAD1, filename.url());
        }

        // Read parfile
        read(path);

        // Parse parfile
        parse();

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load parameters
 *
 * @param[in] filename Parameter filename.
 * @param[in] args Command line arguments. 
 *
 * @exception GException::bad_cmdline_argument
 *            Invalid command line argument encountered.
 *
 * Loads all parameters from parameter file. Parameters are overwritten by
 * the values specified in the command line arguments.
 ***************************************************************************/
void GApplicationPars::load(const GFilename&                filename,
                            const std::vector<std::string>& args)
{
    // Load parameters
    load(filename);

    // Overwrite parameter values that are specified in the command line
    for (int i = 1; i < args.size(); ++i) {

        // Extract parameter name and value (empty values are permitted)
        size_t pos = args[i].find("=");
        if (pos == std::string::npos) {
            throw GException::bad_cmdline_argument(G_LOAD2, args[i],
                                                   "no \"=\" specified");
        }
        std::string name  = args[i].substr(0, pos);
        std::string value = args[i].substr(pos+1);
        if (name.length() < 1) {
            throw GException::bad_cmdline_argument(G_LOAD2, args[i],
                                       "no parameter name before \"=\"");
        }
        
        // Check if parameter exists
        if (!contains(name)) {
            throw GException::bad_cmdline_argument(G_LOAD2, args[i],
                                  "invalid parameter name \""+name+"\"");
        }

        // Assign value
        try {
            (*this)[name].value(value);
        }
        catch (GException::par_error &e) {
            throw GException::bad_cmdline_argument(G_LOAD2, args[i]);
        }

        // Set mode to hidden to prevent querying the parameter
        if ((*this)[name].mode() == "q") {
            (*this)[name].mode("h");
        }
        else if ((*this)[name].mode() == "ql") {
            (*this)[name].mode("hl");
        }
        else if ((*this)[name].mode() == "lq") {
            (*this)[name].mode("lh");
        }

    } // endfor: looped over all parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save parameters
 *
 * @param[in] filename Parameter filename.
 *
 * @exception GException::par_file_not_found
 *            No valid directory to write the parameter file has been found.
 ***************************************************************************/
void GApplicationPars::save(const GFilename& filename)
{
    // Get path to parameter file for output
    std::string path = outpath(filename.url());
    if (path.size() == 0) {
        throw GException::par_file_not_found(G_SAVE, filename.url());
    }

    // Update parameter file
    update();

    // Write parfile
    write(path);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set class from pickled string vector
 *
 * @param[in] string String vector containing class information.
 ***************************************************************************/
void GApplicationPars::pickle(const std::vector<std::string>& string)
{
    // Clear object
    clear();

    // Extract vector sizes
    int n_parfile = gammalib::toint(string[0]);
    int n_pars    = gammalib::toint(string[1]);
    int n_line    = gammalib::toint(string[2]);
    int n_vstart  = gammalib::toint(string[3]);
    int n_vstop   = gammalib::toint(string[4]);

    // Initialise string pointer
    int istring = 5;

    // Extract parameter file lines
    for (int i = 0; i < n_parfile; ++i, ++istring) {
        m_parfile.push_back(string[istring]);
    }

    // Extract application parameters
    for (int i = 0; i < n_pars; ++i) {
        int n_par   = gammalib::toint(string[istring]);
        int i_start = istring + 1;
        int i_end   = i_start + n_par;
        //std::vector<std::string> sub(string.begin()+i_start,
        //                             string.begin()+i_end);
        GApplicationPar par;
        par.pickle(std::vector<std::string>(string.begin()+i_start,
                                            string.begin()+i_end));
        m_pars.push_back(par);
        istring = i_end;
    }

    // Extract line numbers
    for (int i = 0; i < n_line; ++i, ++istring) {
        m_line.push_back(gammalib::toint(string[istring]));
    }

    // Extract columns of start values
    for (int i = 0; i < n_vstart; ++i, ++istring) {
        m_vstart.push_back((size_t)gammalib::toint(string[istring]));
    }

    // Extract columns of stop values
    for (int i = 0; i < n_vstop; ++i, ++istring) {
        m_vstop.push_back((size_t)gammalib::toint(string[istring]));
    }

    // Extract effective mode and location of syspfiles
    m_mode      = string[istring++];
    m_syspfiles = string[istring++];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pickled string vector
 *
 * @return String vector containing class information.
 ***************************************************************************/
std::vector<std::string> GApplicationPars::pickle(void) const
{
    // Allocate vector of strings
    std::vector<std::string> string;

    // Store vector lengths
    string.push_back(gammalib::str(m_parfile.size()));
    string.push_back(gammalib::str(m_pars.size()));
    string.push_back(gammalib::str(m_line.size()));
    string.push_back(gammalib::str(m_vstart.size()));
    string.push_back(gammalib::str(m_vstop.size()));

    // Store parfile lines
    for (int i = 0; i < m_parfile.size(); ++i) {
        string.push_back(m_parfile[i]);
    }

    // Store application parameters
    for (int i = 0; i < m_pars.size(); ++i) {
        std::vector<std::string> par = m_pars[i].pickle();
        string.push_back(gammalib::str(par.size()));
        for (int k = 0; k < par.size(); ++k) {
            string.push_back(par[k]);
        }
    }

    // Store line numbers
    for (int i = 0; i < m_line.size(); ++i) {
        string.push_back(gammalib::str(m_line[i]));
    }

    // Store columns of start values
    for (int i = 0; i < m_vstart.size(); ++i) {
        string.push_back(gammalib::str(m_vstart[i]));
    }

    // Store columns of stop values
    for (int i = 0; i < m_vstop.size(); ++i) {
        string.push_back(gammalib::str(m_vstop[i]));
    }

    // Store effective mode and location of syspfiles
    string.push_back(m_mode);
    string.push_back(m_syspfiles);

    // Return string vector
    return string;
}


/***********************************************************************//**
 * @brief Print parameters
 *
 * @param[in] chatter Chattiness.
 * @return String containing parameter information.
 ***************************************************************************/
std::string GApplicationPars::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GApplicationPars ===");

        // Append parameters
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i].print(chatter));
        }

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
void GApplicationPars::init_members(void)
{
    // Initialise members
    m_parfile.clear();
    m_pars.clear();
    m_line.clear();
    m_vstart.clear();
    m_vstop.clear();
    m_mode = "h";
    m_syspfiles.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pars Object from which members which should be copied.
 ***************************************************************************/
void GApplicationPars::copy_members(const GApplicationPars& pars)
{
    // Copy attributes
    m_parfile   = pars.m_parfile;
    m_pars      = pars.m_pars;
    m_line      = pars.m_line;
    m_vstart    = pars.m_vstart;
    m_vstop     = pars.m_vstop;
    m_mode      = pars.m_mode;
    m_syspfiles = pars.m_syspfiles;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GApplicationPars::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Determine filepath for parameter file input
 *
 * @param[in] filename Parameter filename to search for.
 *
 * Locates parameter file for input. 
 * The parameter file is first searched in the directories that are listed
 * in the PFILES environment variable. Directories may be separated by : or
 * by ; in PFILES.
 *
 * If the PFILES environment variable is not set the parameter file is
 * searched in the users pfiles directory.
 *
 * If still no parameter file is found, the parameter file is searched (in
 * the given order) in ${GAMMALIB}/syspfiles and in ${prefix}/syspfiles,
 * where ${prefix} is the path to the GammaLib installation.
 ***************************************************************************/
std::string GApplicationPars::inpath(const std::string& filename) const
{
    // Allocate result path
    std::string path;

    // Search for parameter file in PFILES directories if the PFILES
    // environment variable has been set
    char* ptr = std::getenv("PFILES");
    if (ptr != NULL) {

        // Extract directories from PFILES environment variable
        std::string              pfiles = ptr;
        std::vector<std::string> dirs   = gammalib::split(pfiles, ":;");

        // Search for first occurence of parameter file
        for (int i = 0; i < dirs.size(); ++i) {

            // Build filename
            std::string fname = dirs[i] + "/" + filename;

            // If file is accessible for reading then exit loop
            if (access(fname.c_str(), R_OK) == 0) {
                path = fname;
                break;
            }

        } // endfor: searched all directories given in PFILES

    } // endif: PFILES environment variable has been set

    // ... otherwise, if no PFILES environment variable has been set
    // then search in users pfiles directory
    else {
        uid_t uid         = geteuid();
        struct passwd* pw = getpwuid(uid);
        if (pw != NULL) {
            std::string fname = std::string(pw->pw_dir) + "/pfiles/" + filename;
            if (access(fname.c_str(), R_OK) == 0) {
                path = fname;
            }
        }    
    } // endif: searched in users pfiles directory

    // If we have no valid path so far then search file within GAMMALIB
    // repository (${GAMMALIB}/syspfiles)
    if (path.size() == 0) {
        ptr = std::getenv("GAMMALIB");
        if (ptr != NULL) {
            std::string fname = std::string(ptr) + "/syspfiles/" + filename;
            if (access(fname.c_str(), R_OK) == 0) {
                path = fname;
            }
        }
    }

    // If we have no valid path so far then search file within GammaLib
    // package (${prefix}/syspfiles)
    #ifdef PACKAGE_PREFIX
    if (path.size() == 0) {
        std::string fname = std::string(PACKAGE_PREFIX) + "/syspfiles/" +
                            filename;
        if (access(fname.c_str(), R_OK) == 0) {
            path = fname;
        }
    }
    #endif

    // Return path
    return path;
}


/***********************************************************************//**
 * @brief Return path to parfile in $PFILES or $HOME/pfiles folder
 *
 * @param[in] filename Parameter filename to search for.
 * @return Full parameter filename, including absolute access path.
 *
 * Locates parameter file in the directories that are listed in the PFILES
 * environment variable. Directories may be separated by : or by ; in the
 * PFILES environment variable.
 *
 * If the PFILES environment variable is not set the parameter file is
 * searched in the users pfiles directory.
 *
 * If parameter file was not found, or the parameter file is not
 * read-accessible, an empty string is returned.
 ***************************************************************************/
std::string GApplicationPars::pfiles_path(const std::string& filename) const
{
    // Allocate result path
    std::string path;

    // Search for parameter file in PFILES directories if the PFILES
    // environment variable has been set
    char* ptr = std::getenv("PFILES");
    if (ptr != NULL) {

        // Extract directories from PFILES environment variable
        std::string              pfiles = ptr;
        std::vector<std::string> dirs   = gammalib::split(pfiles, ":;");

        // Search for first occurence of parameter file
        for (int i = 0; i < dirs.size(); ++i) {

            // Build filename
            std::string fname = dirs[i] + "/" + filename;

            // If file is accessible for reading then exit loop
            if (access(fname.c_str(), R_OK) == 0) {
                path = fname;
                break;
            }

        } // endfor: searched all directories given in PFILES

    } // endif: PFILES environment variable has been set

    // ... otherwise, if no PFILES environment variable has been set
    // then search in users pfiles directory
    else {
        uid_t uid         = geteuid();
        struct passwd* pw = getpwuid(uid);
        if (pw != NULL) {
            std::string fname = std::string(pw->pw_dir) + "/pfiles/" + filename;
            if (access(fname.c_str(), R_OK) == 0) {
                path = fname;
            }
        }    
    } // endif: searched in users pfiles directory

    // Return path
    return path;
}


/***********************************************************************//**
 * @brief Return path to parfile in m_syspfiles folder
 *
 * @param[in] filename Parameter filename to search for.
 * @return Full parameter filename, including absolute access path.
 *
 * Locates parameter file in the directory that is specified by the
 * m_syspfiles string. If the string is empty, the specified file is not
 * found or read accessible, an empty string is returned.
 ***************************************************************************/
std::string GApplicationPars::syspfiles_path(const std::string& filename) const
{
    // Allocate result path
    std::string path;

    // Continue only if m_syspfiles is not empty
    if (!m_syspfiles.empty()) {
        std::string fname = m_syspfiles + "/" + filename;
        if (access(fname.c_str(), R_OK) == 0) {
            path = fname;
        }
    }

    // Return path
    return path;
}


/***********************************************************************//**
 * @brief Determine filepath for parameter file output
 *
 * @param[in] filename Parameter filename.
 *
 * @exception GException::home_not_found
 *            Unable to determine users home directory.
 * @exception GException::could_not_create_pfiles
 *            Unable to create pfiles directory.
 * @exception GException::pfiles_not_accessible
 *            Unable to make pfiles directory accessible to user.
 *
 * Searchs for first writable directory listed in PFILES environment
 * variable. If PFILES is not set then use pfiles directory in users
 * home directory. If pfiles directory does not exist then create it.
 * If directory exists but is not writable then make it writable.
 ***************************************************************************/
std::string GApplicationPars::outpath(const std::string& filename) const
{
   // Allocate result path
    std::string path;

    // Search for writeable PFILES directories
    char* ptr = std::getenv("PFILES");
    if (ptr != NULL) {

        // Extract directories from PFILES environment variable
        std::string              pfiles = ptr;
        std::vector<std::string> dirs   = gammalib::split(pfiles, ":;");

        // Search for first writeable
        for (int i = 0; i < dirs.size(); ++i) {

            // If directory is accessible for writing then exit loop
            if (access(dirs[i].c_str(), W_OK) == 0) {
                path = dirs[i] + "/" + filename;
                break;
            }

        }
    } // endif: PFILES environment variable exists

    // If no valid directory is found in PFILES environment variable then 
    // use pfiles directory in users home directory.
    if (path.size() == 0) {

        // Get users home directory
        uid_t uid         = geteuid();
        gid_t gid         = getegid();
        struct passwd* pw = getpwuid(uid);
        if (pw == NULL) {
            throw GException::home_not_found(G_OUTPATH);
        }

        // Set path
        path = std::string(pw->pw_dir) + "/pfiles";

        // If directory does not exist then create it
        if (access(path.c_str(), F_OK) != 0) {
            if (mkdir(path.c_str(), 
                S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
                throw GException::could_not_create_pfiles(G_OUTPATH, path);
            }
        }

        // If directory exists but is not writable then make it writable
        else if (access(path.c_str(), W_OK) != 0) {
            if (chown(path.c_str(), uid, gid) != 0 ||
                chmod(path.c_str(), 
                      S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
                throw GException::pfiles_not_accessible(G_OUTPATH, path);
            }
        }

        // Append filename
        path = path + "/" + filename;

    } // endif: no valid directory found in PFILES

    // Return path
    return path;
}


/***********************************************************************//**
 * @brief Read parameter file
 *
 * @param[in] filename Parameter filename (absolut path).
 *
 * @exception GException::par_file_open_error
 *            Unable to open parameter file (read access requested).
 *
 * Read all lines of the parameter file. Each line is terminated by a newline
 * character. The file will be locked to avoid simultaneous access by another
 * process.
 ***************************************************************************/
void GApplicationPars::read(const std::string& filename)
{
    // Put in OpenMP critical zone. Note that it is important that the zone
    // has an identical name to the corresponding zone in the write()
    // method so that multiple threads cannot at the same time read and write
    // the parameter file (see #3287).
    #pragma omp critical(GApplicationsPars_io)
    {
        // Allocate line buffer
        const int n = 1000; 
        char  line[n];

        // Trying to get file lock
        #if defined(G_LOCK_PARFILE)
        struct flock lock;
        lock.l_type   = F_RDLCK;  // Want a read lock
        lock.l_whence = SEEK_SET; // Want beginning of file
        lock.l_start  = 0;        // No offset, lock entire file ...
        lock.l_len    = 0;        // ... to the end
        lock.l_pid    = getpid(); // Current process ID
        int fd;
        if ((fd = open(filename.c_str(), O_RDONLY)) == -1) {
            throw GException::par_file_open_error(G_READ, filename,
                  "Could not open file for locking.");
        }
        #if defined(G_CHECK_LOCK_PARFILE)
        if (fcntl(fd, F_SETLKW, &lock) == -1) { // F_SETLKW: wait until unlocked
            throw GException::par_file_open_error(G_READ, filename,
                  "Could not get a lock on the file.");
        }
        #else
        fcntl(fd, F_SETLKW, &lock);
        #endif
        #endif

        // Open parameter file
        FILE* fptr = fopen(filename.c_str(), "r");
        if (fptr == NULL) {
            throw GException::par_file_open_error(G_READ, filename);
        }

        // Read lines
        while (fgets(line, n, fptr) != NULL) {
            m_parfile.push_back(std::string(line));
        }

        // Close file
        fclose(fptr);

        // Unlock file
        #if defined(G_LOCK_PARFILE)
        if (fd != -1) {
            lock.l_type = F_UNLCK;
            #if defined(G_CHECK_LOCK_PARFILE)
            if (fcntl(fd, F_SETLK, &lock) == -1) {
                throw GException::par_file_open_error(G_READ, filename,
                      "Could not unlock the file.");
            }
            #else
            fcntl(fd, F_SETLK, &lock);
            #endif
            close(fd);
        }
        #endif
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write parameter file
 *
 * @param[in] filename Parameter filename (absolut path).
 *
 * @exception GException::par_file_open_error
 *            Unable to open parameter file (write access requested).
 *
 * Writes all lines of the parameter file. The file will be locked to avoid
 * simultaneous access by another process.
 ***************************************************************************/
void GApplicationPars::write(const std::string& filename) const
{
    // Put in OpenMP critical zone Note that it is important that the zone
    // has an identical name to the corresponding zone in the read()
    // method so that multiple threads cannot at the same time read and write
    // the parameter file (see #3287).
    #pragma omp critical(GApplicationsPars_io)
    {
        // Trying to get file lock. We have to do this after opening the file
        // using the fopen function, as the file may not exist, hence it needs
        // to be created first.
        #if defined(G_LOCK_PARFILE)
        struct flock lock;
        lock.l_type   = F_WRLCK;  // Want a write lock
        lock.l_whence = SEEK_SET; // Want beginning of file
        lock.l_start  = 0;        // No offset, lock entire file ...
        lock.l_len    = 0;        // ... to the end
        lock.l_pid    = getpid(); // Current process ID
        int  fd;
        if ((fd = open(filename.c_str(), O_WRONLY)) != -1) {
            #if defined(G_CHECK_LOCK_PARFILE)
            if (fcntl(fd, F_SETLKW, &lock) == -1) { // F_SETLKW: wait until unlocked
                throw GException::par_file_open_error(G_WRITE, filename,
                      "Could not get a lock on the file.");
            }
            #else
            fcntl(fd, F_SETLKW, &lock);
            #endif
        }
        #endif

        // Open parameter file.
        FILE* fptr = fopen(filename.c_str(), "w");
        if (fptr == NULL) {
            throw GException::par_file_open_error(G_WRITE, filename);
        }

        // If file is not locked then lock it now.
        #if defined(G_LOCK_PARFILE)
        if (fd == -1) {
            if ((fd = open(filename.c_str(), O_WRONLY)) != -1) {
                #if defined(G_CHECK_LOCK_PARFILE)
                if (fcntl(fd, F_SETLKW, &lock) == -1) { // F_SETLKW: wait until unlocked
                    fclose(fptr);
                    throw GException::par_file_open_error(G_WRITE, filename,
                          "Could not get a lock on the file.");
                }
                #else
                fcntl(fd, F_SETLKW, &lock);
                #endif
            }
        }
        #endif

        // Write lines
        for (int i = 0; i < m_parfile.size(); ++i) {
            fprintf(fptr, "%s", m_parfile[i].c_str());
        }

        // Close file
        fclose(fptr);

        // Unlock file
        #if defined(G_LOCK_PARFILE)
        if (fd != -1) {
            lock.l_type = F_UNLCK;
            #if defined(G_CHECK_LOCK_PARFILE)
            if (fcntl(fd, F_SETLK, &lock) == -1) {
                throw GException::par_file_open_error(G_WRITE, filename,
                      "Could not unlock the file.");
            }
            #else
            fcntl(fd, F_SETLK, &lock);
            #endif
            close(fd);
        }
        #endif
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse parameter file
 *
 * @exception GException::par_file_syntax_error
 *            Syntax error encountered in parameter file.
 *
 * The parameter type has to be one b,i,r,s,f,fr,fw,fe,fn. The fr,fw,fe,fn
 * types test for read access, write access, file existence, and file
 * absence, respectively.
 * The parameter mode has to be one of a,h,l,q,hl,ql,lh,lq. For mode 'a' the
 * effective mode equals to the value given by the mode parameter, if it
 * exists. Without the presence of the mode parameter the effective mode
 * will be 'h'.
 ***************************************************************************/
void GApplicationPars::parse(void)
{
    // Preset effective mode to 'hidden'
    m_mode = "h";

    // Parse all lines
    for (int i = 0; i < m_parfile.size(); ++i) {

        // Get line without any leading and trailing whitespace
        std::string line = gammalib::strip_whitespace(m_parfile[i]);

        // If line contains a single linefeed then skip the line
        if (line.length() == 1 && line[0] == '\n')
            continue;

        // If line is empty or if line starts with # then skip the line
        if (line.length() == 0 || line[0] == '#') {
            continue;
        }

        // Get the 7 text fields of valid a parameter line
        std::string fields[7];
        int         quotes = 0;
        size_t      start  = 0;
        size_t      end    = line.length() - 1;
        int         index  = 0;
        size_t      vstart = 0;
        size_t      vstop  = 0;
        for (size_t pos = 0; pos < line.length(); ++pos) {

            // Toggle quotes
            if (line[pos] == '"') {
                quotes = 1-quotes;
            }

            // Search for comma only if we are outside quotes. If comma is
            // found or end of line is reached then extract a field and start
            // searching again from the position following the comma. Strip
            // leading and trailing quotes from the fields.
            if (quotes == 0) {
                if (line[pos] == ',' || pos == end) {
                    if (index < 7) {
                        fields[index] = 
                          gammalib::strip_chars(gammalib::strip_whitespace(line.substr(start, 
                                                                 pos-start)), "\"");
                        if (index == 3) {
                            vstart = start;
                            vstop  = pos;
                        }
                        start = pos + 1;
                    }
                    index++;
                }
            }

        } // endfor: looped over line

        // Throw an error if quotes are not balanced
        if (quotes != 0) {
            throw GException::par_file_syntax_error(G_PARSE, 
                                                    gammalib::strip_chars(line,"\n"),
                                                    "quotes are not balanced");
        }

        // Throw an error if line has not 7 fields
        if (index != 7) {
            throw GException::par_file_syntax_error(G_PARSE, 
                                                    gammalib::strip_chars(line,"\n"),
                                 "found "+gammalib::str(index)+" fields, require 7");
        }

        // Verify if parameter name does not yet exist
        if (contains(fields[0])) {
            throw GException::par_file_syntax_error(G_PARSE, 
                                                    gammalib::strip_chars(line,"\n"),
                          "redefinition of parameter name \""+fields[0]+"\"");
        }

        // Add parameter
        try {
            m_pars.push_back(GApplicationPar(fields[0], fields[1], fields[2],
                                             fields[3], fields[4], fields[5],
                                             fields[6]));
            m_line.push_back(i);
            m_vstart.push_back(vstart);
            m_vstop.push_back(vstop);
        }
        catch (GException::par_error &e) {
            throw GException::par_file_syntax_error(G_PARSE, 
                                                    gammalib::strip_chars(line,"\n"),
                                                                 e.what());
        }

        // If parameter name is mode then store the effective mode
        if (fields[0] == "mode") {
            if (fields[3] != "h"  && fields[3] != "q" &&
                fields[3] != "hl" && fields[3] != "ql" &&
                fields[3] != "lh" && fields[3] != "lq") {
                throw GException::par_file_syntax_error(G_PARSE, 
                                                    gammalib::strip_chars(line,"\n"),
                       "mode parameter has invalid value \""+fields[3]+"\"");
            }
            m_mode = fields[3];
        }

    } // endfor: looped over lines

    // Set effective mode for all parameters that have mode 'auto'
    for (int i = 0; i < m_pars.size(); ++i) {
        if (m_pars[i].mode() == "a") {
            m_pars[i].mode(m_mode);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update parameter file
 *
 * Update lines of parameter file according to the parameter values. This
 * method handles correctly formatted parameter files by replacing the value
 * at its original location within the line (preserving additional
 * whitespace).
 * Updating is only done of the parameter mode is 'learn'.
 ***************************************************************************/
void GApplicationPars::update(void)
{
    // Loop over all parameters
    for (int i = 0; i < m_pars.size(); ++i) {

        // Update only if requested and allowed
        if (m_pars[i].m_update && m_pars[i].is_learn()) {
            m_parfile[m_line[i]] = m_parfile[m_line[i]].substr(0, m_vstart[i]) +
                                   m_pars[i].m_value +
                                   m_parfile[m_line[i]].substr(m_vstop[i]);
            m_vstop[i] = m_vstart[i] + m_pars[i].m_value.length();
        }

    } // endfor: looped over all parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Synchronise parameter file with the parameter values in a another
 *        parameter file
 *
 * @param[in] filename File name of other parameter file.
 *
 * Copies over all values from another parameter file that correspond to
 * parameters with the same name and type, and that are signaled as "learn"
 * in the other parameter file.
 ***************************************************************************/
void GApplicationPars::synchronise(const std::string& filename)
{
    // Allocate other application parameters
    GApplicationPars pars;

    // Read other application parameters
    pars.read(filename);

    // Parse other parameters
    pars.parse();

    // Loop over all parameters in parameter file
    for (int i = 0; i < m_pars.size(); ++i) {

        // If a parameter with the same name and type exists in the other
        // parameter file
        int inx = pars.get_index(m_pars[i].name());
        if (inx != -1) {

            // If the parameter types match, the mode in the other parameter
            // file is "learn", and the parameter values differ, then copy
            // over the value and signal value updating
            if ((m_pars[i].type() == pars[inx].type()) &&
                 pars[inx].is_learn() &&
                 m_pars[i].m_value != pars[inx].m_value) {
                m_pars[i].m_value  = pars[inx].m_value;
                m_pars[i].m_status = pars[inx].m_status;
                m_pars[i].m_update = true;
            }

        } // endif: parameter with same name existed

    } // endfor: looped over all parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return parameter index by name
 *
 * @param[in] name Parameter name.
 * @return Parameter index (-1 if parameter name has not been found)
 *
 * Returns parameter index based on the specified @p name. If no parameter
 * with the specified @p name is found the method returns -1.
 ***************************************************************************/
int GApplicationPars::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search parameter with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i].name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Return parameter file line for a specific parameter
 *
 * @param[in] par Parameter.
 * @param[out] start Column of value start.
 * @param[out] stop Column of value stop.
 * @return Parameter file line (termined by \n).
 *
 * Constructs the parameter file line for a specific parameter and returns
 * the line as a string. The line is terminated by a \n character.
 ***************************************************************************/
std::string GApplicationPars::parline(GApplicationPar& par,
                                      size_t*          start,
                                      size_t*          stop) const
{
    // Declate line
    std::string line;

    // Build parameter file line
    line.append(par.name()+", ");
    line.append(par.type()+ ", ");
    line.append(par.mode()+ ",");
    *start = line.length();
    line.append(par.value()+ ",");
    *stop = line.length();
    line.append(par.min()+ ",");
    line.append(par.max()+ ",");
    line.append("\""+par.prompt()+"\"\n");

    // Return line
    return line;
}
