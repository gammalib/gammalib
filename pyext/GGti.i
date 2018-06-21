/***************************************************************************
 *                    GGti.i - Good time interval class                    *
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
 * @file GGti.i
 * @brief Good time interval class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GGti.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GGti
 *
 * @brief Good Time Interval class
 ***************************************************************************/
class GGti : public GContainer {

public:
    // Constructors and destructors
    GGti(void);
    explicit GGti(const GFilename& filename);
    GGti(const GGti& gti);
    GGti(const GTime& tstart, const GTime& tstop);
    explicit GGti(const GXmlElement& xml);
    explicit GGti(const GTimeReference& ref);
    virtual ~GGti(void);

    // Methods
    void                  clear(void);
    GGti*                 clone(void) const;
    std::string           classname(void) const;
    int                   size(void) const;
    bool                  is_empty(void) const;
    void                  append(const GTime& tstart, const GTime& tstop);
    void                  insert(const GTime& tstart, const GTime& tstop);
    void                  merge(void);
    void                  merge(const GTime& tstart, const GTime& tstop);
    void                  reduce(const GTime& tstart, const GTime& tstop);
    void                  remove(const int& index);
    void                  reserve(const int& num);
    void                  extend(const GGti& gti);
    void                  load(const GFilename& filename);
    void                  save(const GFilename& filename,
                               const bool& clobber = false) const;
    void                  read(const GFitsTable& table);
    void                  write(GFits& fits,
                                const std::string& extname = gammalib::extname_gti) const;
    void                  read(const GXmlElement& xml);
    void                  write(GXmlElement& xml) const;
    const GTime&          tstart(void) const;
    const GTime&          tstop(void) const;
    const GTime&          tstart(const int& index) const;
    const GTime&          tstop(const int& index) const;
    const double&         telapse(void) const;
    const double&         ontime(void) const;
    void                  reference(const GTimeReference& ref);
    const GTimeReference& reference(void) const;
    bool                  contains(const GTime& time) const;
};


/***********************************************************************//**
 * @brief GGti class extension
 ***************************************************************************/
%extend GGti {
    GGti* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GGti* gti = new GGti;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        gti->append(self->tstart(i), self->tstop(i));
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        gti->append(self->tstart(i), self->tstop(i));
                    }
                }
                return gti;
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
    GGti copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        args = tuple([self.tstart(i) for i in range(self.size())]), \
               tuple([self.tstop(i)  for i in range(self.size())]), \
               self.reference()
        return args
    def __setstate__(self, state):
        self.__init__()
        size = len(state[0])
        self.reserve(size)
        self.reference(state[2])
        for i in range(size):
            self.append(state[0][i], state[1][i])
}
};
