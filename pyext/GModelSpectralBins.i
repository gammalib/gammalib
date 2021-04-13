/***************************************************************************
 *            GModelSpectralBins.i - Spectral bins model class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GModelSpectralBins.i
 * @brief Spectral bins model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralBins.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralBins
 *
 * @brief Spectral bins model class
 ***************************************************************************/
class GModelSpectralBins : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralBins(void);
    GModelSpectralBins(const GModelSpectral& model,
                       const GEbounds&       ebounds,
                       const double&         index = -2.0);
    explicit GModelSpectralBins(const GXmlElement& xml);
    GModelSpectralBins(const GModelSpectralBins& model);
    virtual ~GModelSpectralBins(void);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralBins* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime& srcTime = GTime(),
                                     const bool& gradients = false) const;
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
    int     bins(void) const;
    void    append(const GEnergy& emin,
                   const GEnergy& emax,
                   const double&  intensity);
    void    insert(const int&     index,
                   const GEnergy& emin,
                   const GEnergy& emax,
                   const double&  intensity);
    void    remove(const int& index);
    void    reserve(const int& num);
    void    extend(const GModelSpectralBins& bins);
    GEnergy emin(const int& index) const;
    GEnergy emax(const int& index) const;
    void    emin(const int& index, const GEnergy& emin);
    void    emax(const int& index, const GEnergy& emax);
    double  intensity(const int& index) const;
    void    intensity(const int& index, const double& intensity);
    double  error(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelSpectralBins class extension
 ***************************************************************************/
%extend GModelSpectralBins {
    GModelSpectralBins copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = xml,
        return state
    def __setstate__(self, state):
        if state[0].elements('bin') == 0:
            self.__init__()
        else:
            self.__init__(state[0])
}
};
