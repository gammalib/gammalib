/***************************************************************************
 *                GPha.i - XSPEC Pulse Height Analyzer class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Juergen Knoedlseder                         *
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
 * @file GPha.i
 * @brief XSPEC Pulse Height Analyzer class
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPha.hpp"
%}

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_pha = "SPECTRUM";
}


/***********************************************************************//**
 * @class GPha
 *
 * @brief Pulse Height Analyzer class
 ***************************************************************************/
class GPha : public GBase {
public:
    // Constructors and destructors
    GPha(void);
    explicit GPha(const GFilename& filename);
    explicit GPha(const GEbounds& ebds);
    explicit GPha(const int& bins);
    GPha(const GPha& pha);
    virtual ~GPha(void);

    // Operators
    GPha&   operator+=(const GPha& pha);
    GPha&   operator-=(const GPha& pha);
    GPha&   operator*=(const double& scale);
    double& operator()(const int& index, const int& col);

    // Methods
    void               clear(void);
    GPha*              clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    int                columns(void) const;
    void               append(const std::string&         name,
                              const std::vector<double>& column);
    const GEbounds&    ebounds(void) const;
    double             counts(void) const;
    GNdarray           counts_spectrum(void) const;
    void               areascal(const int& index, const double& areascal);
    const double&      areascal(const int& index) const;
    void               backscal(const int& index, const double& backscal);
    const double&      backscal(const int& index) const;
    GNdarray           backscal_spectrum(void) const;
    const double&      underflow(void) const;
    const double&      overflow(void) const;
    const double&      outflow(void) const;
    void               exposure(const double& exposure);
    const double&      exposure(void) const;
    void               emin_obs(const GEnergy& emin_obs);
    const GEnergy&     emin_obs(void) const;
    void               emax_obs(const GEnergy& emax_obs);
    const GEnergy&     emax_obs(void) const;
    void               backfile(const std::string& backfile);
    const std::string& backfile(void) const;
    void               corrfile(const std::string& corrfile);
    const std::string& corrfile(void) const;
    void               respfile(const std::string& respfile);
    const std::string& respfile(void) const;
    void               ancrfile(const std::string& ancrfile);
    const std::string& ancrfile(void) const;
    void               fill(const GEnergy& energy, const double& value = 1.0);
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    void               read(const GFits& fits);
    void               read(const GFitsTable& table);
    void               write(GFits& fits) const;
    const GFilename&   filename(void) const;
    const GFitsHeader& header(void) const;
    void               header(const GFitsHeader& header);
};


/***********************************************************************//**
 * @brief GPha class extension
 ***************************************************************************/
%extend GPha {
    double __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    std::vector<double> __getitem__(const std::string& colname) {
        return (*self)[colname];
    }
    void __setitem__(const int& index, const double& value) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = value;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    int __len__() {
        return (self->size());
    }
    GPha __add__(const GPha& pha) const {
        return ((*self) + pha);
    }
    GPha __sub__(const GPha& pha) const {
        return ((*self) - pha);
    }
    GPha __mul__(const double& scale) const {
        return ((*self) * scale);
    }
    // Python 2.x
    GPha __div__(const double& scale) const {
        return ((*self) / scale);
    }
    GPha __idiv__(const double& scale) {
        self->operator/=(scale);
        return (*self);
    }
    // Python 3.x
    GPha __truediv__(const double& scale) const {
        return ((*self) / scale);
    }
    GPha __itruediv__(const double& scale) {
        self->operator/=(scale);
        return (*self);
    }
    GPha copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        fits = gammalib.GFits()
        self.write(fits)
        state = (fits,)
        return state
    def __setstate__(self, state):
        self.__init__()
        if not state[0].is_empty():
            self.read(state[0])
}
};
