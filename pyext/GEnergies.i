/***************************************************************************
 *                    GEnergies.i - Energy container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Juergen Knoedlseder                         *
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
 * @file GEnergies.i
 * @brief Energy container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEnergies.hpp"
#include "GException.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GEnergies
 *
 * @brief Container class for energies.
 ***************************************************************************/
class GEnergies : public GContainer {

public:
    // Constructors and destructors
    GEnergies(void);
    explicit GEnergies(const GFilename& filename);
    explicit  GEnergies(const GEbounds& ebounds);
    GEnergies(const GEnergies& energies);
    GEnergies(const int& num, const GEnergy& emin, const GEnergy& emax,
              const bool& log = true);
    virtual ~GEnergies(void);
 
    // Methods
    void        clear(void);
    GEnergies*  clone(void) const;
    std::string classname(void) const;
    GEnergy&    at(const int& index);
    int         size(void) const;
    bool        is_empty(void) const;
    GEnergy&    append(const GEnergy& energy);
    GEnergy&    insert(const int& index, const GEnergy& energy);
    void        remove(const int& index);
    void        reserve(const int& num);
    void        extend(const GEnergies& energies);
    void        set(const GEbounds& ebounds);
    void        set_lin(const int& num, const GEnergy& emin,
                                        const GEnergy& emax);
    void        set_log(const int& num, const GEnergy& emin,
                                        const GEnergy& emax);
    void        load(const GFilename& filename);
    void        save(const GFilename& filename,
                     const bool&      clobber = false) const;
    void        read(const GFitsTable& table);
    void        write(GFits& file,
                      const std::string& extname = gammalib::extname_energies) const;
};


/***********************************************************************//**
 * @brief GEnergies class extension
 ***************************************************************************/
%extend GEnergies {
    GEnergy& __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    GEnergies* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GEnergies* energies = new GEnergies;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        energies->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        energies->append((*self)[i]);
                    }
                }
                return energies;
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
    void __setitem__(const int& index, const GEnergy& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GEnergies copy() {
        return (*self);
    }
};
