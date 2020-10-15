/***************************************************************************
 *              GModel.i - Abstract virtual model base class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2020 by Juergen Knoedlseder                         *
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
 * @file GModel.i
 * @brief Abstract model base class python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModel.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModel
 *
 * @brief Abstract model base class python interface
 ***************************************************************************/
class GModel : public GBase {

public:
    // Constructors and destructors
    GModel(void);
    explicit GModel(const GXmlElement& xml);
    GModel(const GModel& model);
    virtual ~GModel(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModel*     clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual bool        is_constant(void) const = 0;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs,
                             const bool& gradients = false) const = 0;
    virtual GVector     eval(const GObservation& obs,
                             GMatrixSparse* gradients = NULL) const = 0;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;

    // Implemented methods
    int                     size(void) const;
    int                     scales(void) const;
    bool                    has_par(const std::string& name) const;
    bool                    has_scales(void) const;
    const std::string&      name(void) const;
    void                    name(const std::string& name);
    const double&           ts(void) const;
    void                    ts(const double& ts);
    const bool&             tscalc(void) const;
    void                    tscalc(const bool& tscalc);
    const bool&             has_ts(void) const;
    std::string             instruments(void) const;
    void                    instruments(const std::string& instruments);
    GModelPar&              scale(const int& index);
    GModelPar               scale(const std::string& instrument) const;
    void                    scale(const GModelPar& par);
    std::string             ids(void) const;
    void                    ids(const std::string& ids);
    bool                    is_valid(const std::string& instruments,
                                     const std::string& ids) const;
    const bool&             has_eval_indices(void) const;
    const std::vector<int>& eval_indices(void) const;
};


/***********************************************************************//**
 * @brief GModel class extension
 ***************************************************************************/
%extend GModel {
    GModelPar& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    GModelPar& __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GModelPar& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    void __setitem__(const std::string& name, const GModelPar& val) {
        (*self)[name] = val;
        return;
    }
%pythoncode {
    def __getstate__(self):
        state = (self.name(), self.instruments(), self.ids(), self.tscalc(),
                 self.has_ts(), self.ts())
        return state
    def __setstate__(self, state):
        self.__init__()
        self.name(state[0])
        self.instruments(state[1])
        self.ids(state[2])
        self.tscalc(state[3])
        if state[4]:
            self.ts(state[5])
}
};
