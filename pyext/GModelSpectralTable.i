/***************************************************************************
 *           GModelSpectralTable.i - Spectral table model class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTable.i
 * @brief Spectral table model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralTable.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralTable
 *
 * @brief Spectral table model class
 ***************************************************************************/
class GModelSpectralTable : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralTable(void);
    GModelSpectralTable(const GFilename& filename, const double& norm);
    GModelSpectralTable(const GEbounds&                ebounds,
                        const GModelSpectralTablePars& pars,
                        const GNdarray&                spectra);
    explicit GModelSpectralTable(const GXmlElement& xml);
    GModelSpectralTable(const GModelSpectralTable& model);
    virtual ~GModelSpectralTable(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralTable* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng,
                                      const GTime&   srcTime = GTime(),
                                      const bool&    gradients = false) const;
    virtual double               flux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin,
                                       const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin,
                                    const GEnergy& emax,
                                    const GTime&   time,
                                    GRan&          ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;

    // Other methods
    int                          nspectra(void) const;
    GModelSpectralTablePar&      table_par(const int& index);
    GModelSpectralTablePar&      table_par(const std::string& name);
    double                       norm(void) const;
    void                         norm(const double& norm);
    const GEbounds&              ebounds(void) const;
    void                         load(const GFilename& filename);
    void                         save(const GFilename& filename,
                                      const bool&      clobber = false) const;
    const GFilename&             filename(void) const;
    void                         energy_scale(const std::string& name);
};


/***********************************************************************//**
 * @brief GModelSpectralTable class extension
 ***************************************************************************/
%extend GModelSpectralTable {
    GModelSpectralTable copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = xml,
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
