/***************************************************************************
 *          GModelSpectralFunc.i - Spectral function model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralFunc.i
 * @brief Spectral function model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralFunc.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralFunc
 *
 * @brief Spectral function model class
 ***************************************************************************/
class GModelSpectralFunc : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralFunc(void);
    GModelSpectralFunc(const GFilename& filename, const double& norm);
    GModelSpectralFunc(const GModelSpectral& model, const GEnergies& energies);
    explicit GModelSpectralFunc(const GXmlElement& xml);
    GModelSpectralFunc(const GModelSpectralFunc& model);
    virtual ~GModelSpectralFunc(void);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralFunc* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime&   srcTime = GTime(),
                                     const bool&    gradients = false) const;
    virtual double              flux(const GEnergy& emin,
                                     const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin,
                                   const GEnergy& emax,
                                   const GTime&   time,
                                   GRan&          ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    int              nodes(void) const;
    bool             is_empty(void) const;
    void             append(const GEnergy& energy, const double& intensity);
    void             insert(const GEnergy& energy, const double& intensity);
    void             remove(const int& index);
    void             reserve(const int& num);
    void             extend(const GModelSpectralFunc& filefct);
    GEnergy          energy(const int& index) const;
    void             energy(const int& index, const GEnergy& energy);
    double           intensity(const int& index) const;
    void             intensity(const int& index, const double& intensity);
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
    double           norm(void) const;
    void             norm(const double& norm);
    void             save(const GFilename& filename,
                          const bool&      clobber = false) const;
};


/***********************************************************************//**
 * @brief GModelSpectralFunc class extension
 ***************************************************************************/
%extend GModelSpectralFunc {
    GModelSpectralFunc copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state       = {'nodes':       self.nodes(),
                       'energies':    [self.energy(i)    for i in range(self.nodes())],
                       'intensities': [self.intensity(i) for i in range(self.nodes())],
                       'filename':    self.filename(),
                       'norm':        self.norm()}
        return state
    def __setstate__(self, state):
        self.__init__()
        if state['filename'].is_empty():
            for i in range(state['nodes']):
                self.append(state['energies'][i], state['intensities'][i])
        else:
            self.filename(state['filename'])
        self.norm(state['norm'])
}
};
