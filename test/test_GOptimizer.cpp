/***************************************************************************
 *             test_GOptimizer.cpp  -  test GOptimizer class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <stdlib.h>
#include "test_GOptimizer.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *  Test: Model parameter handling                                         *
 ***************************************************************************/
void MyOptFct::init_members(void)
{
    m_event = GData::iterator();
    m_end   = GData::iterator();
    return;
}
void MyOptFct::copy_members(const MyOptFct& fct)
{
    m_data  = fct.m_data;
    m_event = fct.m_event;
    m_end   = fct.m_end;
    return;
}
MyOptFct& MyOptFct::operator= (const MyOptFct& fct)
{
    if (this != &fct) {
        free_members();
        init_members();
        copy_members(fct);
    }
    return *this;
}
void MyOptFct::first_item(void)
{
    m_event = m_data.begin();
    m_end   = m_data.end();
    return;
}
bool MyOptFct::next_item(void)
{
    m_event++;
    return (m_event == m_end);
}
GVector MyOptFct::get_gradient(void) 
{
    GVector gradient = GVector(10);
    return (gradient);
}
void MyOptFct::set_data(const GData& data)
{
    m_data  = data;
    m_event = m_data.begin();
    m_end   = m_data.end();
    return;
}


/***************************************************************************
 *  Test: Optimizer                                                        *
 ***************************************************************************/
void test_optimizer(void)
{
    // Write header
    std::cout << "Test GOptimizer: ";
    
    // Number of observations in data
    int nobs = 50;

    // Setup GData for optimizing
    GData           data;
    GLATObservation obs;
    try {
        obs.load_binned("data/lat/cntmap.fits.gz", "", "");
        for (int i = 0; i < nobs; ++i)
            data.add(obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup GData for optimizing." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
//std::cout << data << std::endl;


    // Setup optimizer
    MyOptFct opt;
    try {
        opt.set_data(data);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup optimizer." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Time simple GData iterator
    double  t_elapse;
    try {
        clock_t t_start = clock();
        int num = 0;
        int sum = 0;
        for (GData::iterator event = data.begin(); event != data.end(); ++event) {
            num++;
            sum += (int)event->counts();
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GData." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << "." << std::endl;
    std::cout << " - Reference time for GData::iterator: " << t_elapse << std::endl;

    // Iterate over optimizer
    try {
        clock_t t_start = clock();
        int     num   = 0;
        int     sum   = 0;
        int     sum_m = 0;
        for (GOptimizerFunction::iterator item = opt.begin(); item != opt.end(); ++item) {
            num++;
            sum += (int)item.data();
//            double  model    = item.model();
//            GVector gradient = item.gradient();
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
        if (sum != 1359*nobs || num != 200000*nobs) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in MyOptFct::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << " - Time for MyOptFct::iterator doing same operations: " << t_elapse 
              << std::endl;
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "****************************" << std::endl;
    std::cout << "* GOptimizer class testing *" << std::endl;
    std::cout << "****************************" << std::endl;

    // Execute the tests
    test_optimizer();

    // Return
    return 0;
}
