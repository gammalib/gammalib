/***************************************************************************
 *    GModelSpectralExpInvPlaw.i - Exponential cut off power law model     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2018 by Alexander Ziegler                           *
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
 * @file GModelSpectralExpInvPlaw.i
 * @brief Exponential cut off power law spectral class interface definition
 * @author Alexander Ziegler
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralExpInvPlaw.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralExpInvPlaw
 *
 * @brief Exponential cut off power law spectral class
 ***************************************************************************/
class GModelSpectralExpInvPlaw : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralExpInvPlaw(void);
    GModelSpectralExpInvPlaw(const std::string& type,
                             const std::string& prefactor,
                             const std::string& index,
                             const std::string& pivot,
                             const std::string& inverse_cutoff);
    GModelSpectralExpInvPlaw(const double&  prefactor,
                             const double&  index,
                             const GEnergy& pivot,
                             const double&  inverse_cutoff);
    GModelSpectralExpInvPlaw(const double&  prefactor,
                             const double&  index,
                             const GEnergy& pivot,
                             const GEnergy& cutoff);
    explicit GModelSpectralExpInvPlaw(const GXmlElement& xml);
    GModelSpectralExpInvPlaw(const GModelSpectralExpInvPlaw& model);
    virtual ~GModelSpectralExpInvPlaw(void);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpectralExpInvPlaw* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GEnergy& srcEng,
                                           const GTime&   srcTime = GTime(),
                                           const bool&    gradients = false) const;
    virtual double                    flux(const GEnergy& emin,
                                           const GEnergy& emax) const;
    virtual double                    eflux(const GEnergy& emin,
                                            const GEnergy& emax) const;
    virtual GEnergy                   mc(const GEnergy& emin,
                                         const GEnergy& emax,
                                         const GTime&   time,
                                         GRan&          ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    void    type(const std::string& type);
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index(void) const;
    void    index(const double& index);
    double  inverse_cutoff(void) const;
    void    inverse_cutoff(const double& alpha);
    GEnergy cutoff(void) const;
    void    cutoff(const GEnergy& cutoff);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);
};


/***********************************************************************//**
 * @brief GModelSpectralExpInvPlaw class extension
 ***************************************************************************/
%extend GModelSpectralExpInvPlaw {
    GModelSpectralExpInvPlaw copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.type(), self[0], self[1], self[2], self[3])
        return state
    def __setstate__(self, state):
        self.__init__()
        self.type(state[0])
        self[0] = state[1]
        self[1] = state[2]
        self[2] = state[3]
        self[3] = state[4]
}
};
