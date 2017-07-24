/***************************************************************************
 *                GCTAEventList.i - CTA event list class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2017 by Juergen Knoedlseder                         *
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
 * @file GCTAEventList.i
 * @brief CTA event list class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventList.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventList
 *
 * @brief CTA event list class Python interface
 ***************************************************************************/
class GCTAEventList : public GEventList {

public:
    // Constructors and destructors
    GCTAEventList(void);
    explicit GCTAEventList(const GFilename& filename);
    GCTAEventList(const GCTAEventList& list);
    virtual ~GCTAEventList(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCTAEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& fits);
    virtual void           write(GFits& fits) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GCTARoi& roi(void) const;

    // Implement other methods
    void               append(const GCTAEventAtom& event);
    void               append_column(const GFitsTableCol& column);
    void               reserve(const int& number);
    void               remove(const int& index, const int& number = 1);
    void               write(GFits& fits,
                             const std::string& evtname,
                             const std::string& gtiname) const;
    void               fetch(void) const;
    void               dispose(void) const;
    double             irf_cache(const std::string& name,
                                 const int& index) const;
    void               irf_cache(const std::string& name,
                                 const int& index,
                                 const double& irf) const;
    const GPhases&     phases(void) const;
    void               phases(const GPhases& phases);
    const std::string& gtiname(void) const;
    void               has_phase(const bool& has_phase);
    void               has_detxy(const bool& has_detxy);
    void               has_mc_id(const bool& has_mc_id);
    const bool&        has_phase() const;
    const bool&        has_detxy() const;
    const bool&        has_mc_id() const;
    void               set_mc_id_names(const std::vector<int>&         ids,
                                       const std::vector<std::string>& names);
};


/***********************************************************************//**
 * @brief GCTAEventList class extension
 ***************************************************************************/
%extend GCTAEventList {
    GCTAEventAtom* __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Event index",
                                           index, self->size());
        }
    }
    GCTAEventList* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GCTAEventList* list = new GCTAEventList;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        list->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        list->append(*(*self)[i]);
                    }
                }
                return list;
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
    void __setitem__(const int& index, const GCTAEventAtom& event) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            *(*self)[index] = event;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            *(*self)[self->size()+index] = event;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "Event index",
                                           index, self->size());
        }
    }
    GCTAEventList copy() {
        return (*self);
    }
};
