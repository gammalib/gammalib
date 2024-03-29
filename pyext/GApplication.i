/***************************************************************************
 *              GApplication.i - GammaLib application base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2023 by Juergen Knoedlseder                         *
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
 * @file GApplication.i
 * @brief GammaLib application base class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#define SWIG
#include "GApplication.hpp"
%}

// Include (int ARGC, char** ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application Python interface definition.
 *
 * The GApplication class is the base class for all tools and scripts that
 * will be based on GammaLib.
 ***************************************************************************/
class GApplication : public GBase {
public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name,
                 const std::string& version);
    GApplication(const std::string&      name,
                 const std::string&      version,
                 const GApplicationPars& pars);
    GApplication(const std::string& name,
                 const std::string& version,
                 int                ARGC,
                 char               **ARGV);
    GApplication(const GApplication& app);
    ~GApplication(void);

    // Public methods
    double telapse(void) const;
    double celapse(void) const;
    double gCO2e(const std::string& country) const;
    void   logFileOpen(const bool& clobber = true);
    void   logFileClose(void);

    // Ignore base class methods and make methods private in Python by
    // prepending an underscore
    %ignore                 clear;
    %ignore                 clone;
    %ignore                 classname;
    %rename(_name)          name;
    %rename(_version)       version;
    %rename(_add_celapse)   add_celapse;
    %rename(_logTerse)      logTerse;
    %rename(_logNormal)     logNormal;
    %rename(_logExplicit)   logExplicit;
    %rename(_logVerbose)    logVerbose;
    %rename(_logDebug)      logDebug;
    %rename(_clobber)       clobber;
    %rename(_par_filename)  par_filename;
    %rename(_log_filename)  log_filename;
    %rename(_log_header)    log_header;
    %rename(_log_trailer)   log_trailer;
    %rename(_need_help)     need_help;
    %rename(_statistics)    statistics;
    %rename(_args)          args;
    %rename(_stamp)         stamp;
    %rename(_log)           log;

    // Public methods
    void                            clear(void);
    GApplication*                   clone(void) const;
    std::string                     classname(void) const;
    const std::string&              name(void) const;
    const std::string&              version(void) const;
    void                            add_celapse(const double& celapse);
    bool                            logTerse(void) const;
    bool                            logNormal(void) const;
    bool                            logExplicit(void) const;
    bool                            logVerbose(void) const;
    bool                            logDebug(void) const;
    bool                            clobber(void) const;
    bool                            has_par(const std::string& name) const;
    const std::string&              par_filename(void) const;
    const std::string&              log_filename(void) const;
    void                            log_header(void);
    void                            log_trailer(void);
    const bool&                     need_help(void) const;
    void                            statistics(const bool& statistics);
    const bool&                     statistics(void) const;
    const GApplicationPars&         pars(void) const;
    void                            pars(const GApplicationPars& pars);
    void                            stamp(GFitsHDU& hdu) const;
    void                            stamp(GFits& fits) const;
    void                            stamp(const GFilename& filename) const;
    const std::vector<std::string>& args(void) const;

    // Public members
    GLog log;   //!< Application logger

    // Make methods private in Python by prepending an underscore
    %rename(_running) running;

    // Protected methods
    static int& running();
};


/***********************************************************************//**
 * @brief GApplication Python class extension
 ***************************************************************************/
%pythoncode %{
# Log the value of a parameter
def _log_value(self, chatter, name, value, *unit):
    if len(unit) > 0:
        string = gammalib.parformat(str(name))+str(value)+' '+str(unit[0])
    else:
        string = gammalib.parformat(str(name))+str(value)
    self._log_string(chatter, string)
GApplication._log_value = _log_value
%}


/***********************************************************************//**
 * @brief GApplication C++ class extension
 ***************************************************************************/
%extend GApplication {
    void _log_string(const int& chatter, const std::string& string,
                     const bool& linefeed = true) {
        self->log_string(GChatter(chatter), string, linefeed);
    }
    void _log_header1(const int& chatter, const std::string& header) {
        self->log_header1(GChatter(chatter), header);
    }
    void _log_header2(const int& chatter, const std::string& header) {
        self->log_header2(GChatter(chatter), header);
    }
    void _log_header3(const int& chatter, const std::string& header) {
        self->log_header3(GChatter(chatter), header);
    }
    void _log_parameters(const int& chatter) {
        self->log_parameters(GChatter(chatter));
    }
    GApplicationPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const std::string& name, const GApplicationPar& val) {
        (*self)[name] = val;
        return;
    }
    void __setitem__(const std::string& name, const bool& val) {
        GApplicationPar& par = (*self)[name];
        if (par.type() == "b") {
            par.boolean(val);
        }
        else {
            std::string msg = "Attempt to set \""+par.type()+
                              "\" parameter \""+name+"\" with boolean "
                              "value \""+gammalib::str(val)+"\".";
            throw GException::invalid_argument("__setitem__(std::string, int)",
                                               msg);
        }
        return;
    }
    void __setitem__(const std::string& name, const int& val) {
        GApplicationPar& par = (*self)[name];
        if (par.type() == "r") {
            par.real(val);
        } 
        else if (par.type() == "i") {
            par.integer(val);
        }
        else if (par.type() == "b") {
            par.boolean(val);
        }
        else if (par.type() == "t") {
            par.value(gammalib::str(val));
        }
        else {
            std::string msg = "Attempt to set \""+par.type()+
                              "\" parameter \""+name+"\" with integer "
                              "value \""+gammalib::str(val)+"\".";
            throw GException::invalid_argument("__setitem__(std::string, int)",
                                               msg);
        }
        return;
    }
    void __setitem__(const std::string& name, const double& val) {
        GApplicationPar& par = (*self)[name];
        if (par.type() == "t") {
            par.value(gammalib::str(val));
        }
        else {
            par.real(val);
        }
        return;
    }
    void __setitem__(const std::string& name, const std::string& val) {
        GApplicationPar& par = (*self)[name];
        if (par.type() == "s") {
            par.string(val);
        } 
        else if (par.type().substr(0,1) == "f") {
            par.filename(val);
        }
        else if ((par.type() == "r") || (par.type() == "i")) {
            std::string lval = gammalib::tolower(val);
            if (lval == "indef"     ||
                lval == "none"      ||
                lval == "undef"     ||
                lval == "undefined" ||
                lval == "inf"       ||
                lval == "infinity"  ||
                lval == "nan") {
                par.value(val);
            }
            else {
                std::string msg = "Attempt to set \""+par.type()+
                                  "\" parameter \""+name+"\" with string "
                                  "value \""+val+"\".";
                throw GException::invalid_argument("__setitem__(std::string, std::string)",
                                                   msg);
            }
        } 
        else if (par.type() == "t") {
            par.value(val);
        } 
        else {
            std::string msg = "Attempt to set \""+par.type()+
                              "\" parameter \""+name+"\" with string "
                              "value \""+val+"\".";
            throw GException::invalid_argument("__setitem__(std::string, std::string)",
                                               msg);
        } 
        return;
    }
    void __setitem__(const std::string& name, const GTime& val) {
        GApplicationPar& par = (*self)[name];
        if (par.type() == "t") {
            par.time(val);
        }
        else {
            std::string msg = "Attempt to set \""+par.type()+
                              "\" parameter \""+name+"\" with time "
                              "value \""+val.print()+"\".";
            throw GException::invalid_argument("__setitem__(std::string, GTime)",
                                               msg);
        }
        return;
    }
    GApplicationPar& __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->pars().size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->pars().size()+index >= 0) {
            return (*self)[self->pars().size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Application parameter index",
                                           index, self->pars().size());
        }
    }
    void __setitem__(const int& index, const GApplicationPar& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->pars().size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->pars().size()+index >= 0) {
            (*self)[self->pars().size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Application parameter index",
                                           index, self->pars().size());
        }
    }
    GApplication copy() {
        return (*self);
    }
    void _inc_running(void) {
        self->running()++;
        return;
    }
    void _dec_running(void) {
        self->running()--;
        return;
    }
%pythoncode {
    def __getstate__(self):
        pars = self.pars()
        if pars.contains('logfile'):
            pars['logfile'].filename('')
        state = (self._name(), self._version(), pars)
        return state
    def __setstate__(self, state):
        if state[0] and state[1]:
            self.__init__(state[0], state[1], state[2])
        else:
            self.__init__()
}
};
