/***************************************************************************
 *                 GEbounds.i - Energy boundaries class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2022 by Juergen Knoedlseder                         *
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
 * @file GEbounds.i
 * @brief Energy boundaries class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEbounds.hpp"
#include "GTools.hpp"
%}

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_ebounds = "EBOUNDS";
}


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Energy boundaries container class
 ***************************************************************************/
class GEbounds : public GContainer {
public:
    // Constructors and destructors
    GEbounds(void);
    explicit GEbounds(const GFilename& filename);
    explicit GEbounds(const GXmlElement& xml);
    explicit GEbounds(const GEnergies& energies);
    GEbounds(const GEbounds& ebds);
    GEbounds(const GEnergy& emin, const GEnergy& emax);
    GEbounds(const int&         num,
             const GEnergy&     emin,
             const GEnergy&     emax,
             const std::string& method = "LOG",
             const double&      gamma  = 1.0);
    virtual ~GEbounds(void);

    // Methods
    void           clear(void);
    GEbounds*      clone(void) const;
    std::string    classname(void) const;
    int            size(void) const;
    bool           is_empty(void) const;
    void           append(const GEnergy& emin, const GEnergy& emax);
    void           insert(const GEnergy& emin, const GEnergy& emax);
    void           merge(void);
    void           merge(const GEnergy& emin, const GEnergy& emax);
    void           remove(const int& index);
    void           remove(const GEnergy& emin, const GEnergy& emax);
    void           reserve(const int& num);
    void           extend(const GEbounds& ebds);
    void           set(const GEnergies& energies);
    void           set(const int&         num,
                       const GEnergy&     emin,
                       const GEnergy&     emax,
                       const std::string& method = "LOG",
                       const double&      gamma  = 1.0);
    void           load(const GFilename& filename);
    void           save(const GFilename& filename,
                        const bool& clobber = false,
                        const std::string& unit = "keV") const;
    void           read(const GFitsTable& table);
    void           write(GFits& file,
                         const std::string& extname = gammalib::extname_ebounds,
                         const std::string& unit = "keV") const;
    void           read(const GXmlElement& xml);
    void           write(GXmlElement& xml) const;
    int            index(const GEnergy& eng) const;
    const GEnergy& emin(void) const;
    const GEnergy& emax(void) const;
    void           emin(const int& index, const GEnergy& energy);
    void           emax(const int& index, const GEnergy& energy);
    GEnergy        emin(const int& index) const;
    GEnergy        emax(const int& index) const;
    GEnergy        emean(const int& index) const;
    GEnergy        elogmean(const int& index) const;
    GEnergy        ewidth(const int& index) const;
    bool           contains(const GEnergy& eng) const;
    bool           contains(const GEnergy& emin, const GEnergy& emax) const;
};


/***********************************************************************//**
 * @brief GEbounds class extension
 ***************************************************************************/
%extend GEbounds {
    GEbounds* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GEbounds* ebounds = new GEbounds;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        ebounds->append(self->emin(i), self->emax(i));
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        ebounds->append(self->emin(i), self->emax(i));
                    }
                }
                return ebounds;
            }
            else {
                throw GException::invalid_argument("__getitem__(PyObject)",
                                                   "Invalid slice indices");
            }
        }
        else {
            throw GException::invalid_argument("__getitem__(PyObject)","");
        }
    }
    bool __eq__(const GEbounds& bounds) const {
        return ((*self) == bounds);
    }
    bool __ne__(const GEbounds& bounds) const {
        return ((*self) != bounds);
    }
    GEbounds copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = tuple([self.emin(i) for i in range(self.size())]), \
                tuple([self.emax(i) for i in range(self.size())])
        return state
    def __setstate__(self, state):
        self.__init__()
        size = len(state[0])
        self.reserve(size)
        for i in range(size):
            self.append(state[0][i], state[1][i])
}
};
