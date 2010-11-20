/***************************************************************************
 *                     GException.i  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
%module GException

%{
/* Put headers and other declarations here */
#include "GException.hpp"
%}
%include exception.i

%exception {
    try {
        $action
    }
    catch (const GException::mem_alloc& e) {
        SWIG_exception(SWIG_MemoryError, e.what());
    }
    catch (const GException::empty& e) {
        SWIG_exception(SWIG_ValueError, e.what());
    }
    catch (const GException::fits_invalid_type& e) {
        SWIG_exception(SWIG_ValueError, e.what());
    }
    catch (const GException::out_of_range& e) {
        SWIG_exception(SWIG_IndexError, e.what());
    }
    catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}
