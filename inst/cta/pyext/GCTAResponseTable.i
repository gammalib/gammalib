/***************************************************************************
 *              GCTAResponseTable.i - CTA response table class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseTable.i
 * @brief CTA response table class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTools.hpp"
#include "GCTAResponseTable.hpp"
%}

/* Define std::vector<double> as valid return type (otherwise a memory leak
   occurs. */
%include "std_vector.i"
%template(VecDouble) std::vector<double>;


/***********************************************************************//**
 * @brief Tuple to index conversion to provide pixel access.
 *
 * The following function provides conversion between a Python tuple and
 * an integer array. This allows CTAResponseTable pixel access via tuples, such as in
 * a[3,5,10] = 10.0 or c = a[2,9].
 ***************************************************************************/
%{
static int rsp_table_tuple(PyObject *input, int *ptr) {
    if (PySequence_Check(input)) {
        int size = PyObject_Length(input);
        if (size > 2) {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return 0;
        }
        ptr[0] = size;
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem(input,i);
            if (!PyInt_Check(o)) {
                Py_XDECREF(o);
                PyErr_SetString(PyExc_ValueError,"Expecting a tuple of integers");
                return 0;
            }
            ptr[i+1] = (int)PyInt_AsLong(o);
            Py_DECREF(o);
        }
        return 1;
    }
    else {
        ptr[0] = 1;
        if (!PyInt_Check(input)) {
            PyErr_SetString(PyExc_ValueError,"Expecting an integer");
            return 0;
        }
        ptr[1] = (int)PyInt_AsLong(input);
        return 1;       
    }
}
%}

// This is the typemap that makes use of the function defined above
%typemap(in) int GCTAResponseTableInx[ANY](int temp[3]) {
   if (!rsp_table_tuple($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}

// This typecheck verifies that all arguments are integers. The typecheck
// is needed for using "int GCTAResponseTableInx" in overloaded methods.
%typemap(typecheck) int GCTAResponseTableInx[ANY] {
    $1 = 1;
    if (PySequence_Check($input)) {
        int size = PyObject_Length($input);
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem($input,i);
            if (!PyInt_Check(o)) {
                $1 = 0;
                break;
            }
        }
    }
    else {
        if (!PyInt_Check($input)) {
            $1 = 0;
        }
    }
}


/***********************************************************************//**
 * @class GCTAResponseTable
 *
 * @brief Interface for the CTA response table class
 ***************************************************************************/
class GCTAResponseTable : public GBase {

public:
    // Constructors and destructors
    GCTAResponseTable(void);
    GCTAResponseTable(const GCTAResponseTable& table);
    explicit GCTAResponseTable(const GFitsTable& hdu);
    virtual ~GCTAResponseTable(void);

    // Methods
    void               clear(void);
    GCTAResponseTable* clone(void) const;
    int                size(void) const;
    const int&         elements(void) const;
    const int&         axes(void) const;
    int                axis(const int& index) const;
    double             axis_lo(const int& index, const int& bin) const;
    double             axis_hi(const int& index, const int& bin) const;
    void               axis_linear(const int& index);
    void               axis_log10(const int& index);
    void               axis_radians(const int& index);
    std::string        axis_lo_name(const int& index) const;
    std::string        axis_hi_name(const int& index) const;
    std::string        axis_lo_unit(const int& index) const;
    std::string        axis_hi_unit(const int& index) const;
    std::string        unit(const int& index) const;
    const GNodeArray&  nodes(const int& index) const;
    void               scale(const int& index, const double& scale);
    void               read(const GFitsTable& hdu);
    void               write(GFitsTable& hdu) const;
    void               add_axis(std::vector<double> axis_lo, 
				std::vector<double> axis_hi,
				std::string name_lo, std::string name_hi,
				std::string unit);    
    void               add_par(std::string name, std::string unit);
};






/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponseTable {
    GCTAResponseTable copy() {
        return (*self);
    }
    double __getitem__(int GCTAResponseTableInx[]) {
        if (GCTAResponseTableInx[0] == 1) {
            return (*self)(GCTAResponseTableInx[1]);
        }
        else {
            return (*self)(GCTAResponseTableInx[1], GCTAResponseTableInx[2]);
        }
    }
    void __setitem__(int GCTAResponseTableInx[], double value) {
        if (GCTAResponseTableInx[0] == 1) {
            (*self)(GCTAResponseTableInx[1]) = value;
        }
        else {
            (*self)(GCTAResponseTableInx[1], GCTAResponseTableInx[2]) = value;
        }
    } 
};
