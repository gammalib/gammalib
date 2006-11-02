/***************************************************************************
 *                  test_GVector.cpp  -  test vector class                 *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "test_GVector.hpp"

using namespace std;

int main(void)
{
  // Dump header
  cout << endl;
  cout << "*************************" << endl;
  cout << "* GVector class testing *" << endl;
  cout << "*************************" << endl;

  // Set parameters
  int num = 5;

  // Define vectors
  GVector test(num);
  GVector result(num);
  GVector smaller(num+1);
  GVector bigger(num+1);
  for (int i = 0; i < num; ++i)
    test(i) = (i+1) * 1.1;
  for (int i = 0; i < num-1; ++i)
    smaller(i) = (i+1) * 1.1;
  for (int i = 0; i < num+1; ++i)
    bigger(i) = (i+1) * 1.1;

  // Test 1: Allocate zero vector
  try {
    cout << "Test 1: Allocate zero vector: " << endl;
    GVector test1(0);
	cout << "  Size = " << test1.size() << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }
  
  // Test 2: Allocate too large vector
  try {
    cout << "Test 2: Allocate too large vector: " << endl;
    GVector test2(1000000000);
	cout << "  Size = " << test2.size() << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }

  // Test 3: Assign values
  try {
    cout << "Test 3: Assign values:" << endl;
    //
	cout << "  test3(1) = pi: ";
    GVector test3(3);
	test3(1) = acos(-1.0);
	cout << test3 << endl;
	cout << "  Size = " << test3.size() << endl;
    //
	cout << "  GVector out of range access: ";
    #if defined(G_RANGE_CHECK)
	cout << test3(3) << endl;
	#endif
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }

  // Test 4: Vector copy constructor
  try {
    cout << "Test 4: Define vector using copy constructor:" << endl;
	//
	cout << "  GVector: ";
	cout << test << endl;
	//
	cout << "  GVector test = GVector: ";
	GVector test4 = test;
	cout << test4 << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }

  // Test 5: Vector assignment
  try {
    cout << "Test 5: Vector assignment:" << endl;
	//
	cout << "  GVector: ";
	cout << test << endl;
	//
	cout << "  GVector = GVector: ";
	result = test;
	cout << result << endl;
	//
	cout << "  GVector = GVector (bigger vector): ";
	result = bigger;
	cout << result << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }

  // Test 6: Assignment and arithmetics
  try {
    cout << "Test 6: Assignment and arithmetics:" << endl;
	//
	cout << "  GVector: ";
	cout << test << endl;
	//
	cout << "  GVector += GVector: ";
	result  = test;
	result += test;
	cout << result << endl;
	//
	cout << "  GVector += 2.0: ";
	result  = test;
	result += 2.0;
	cout << result << endl;
	//
	cout << "  GVector -= GVector: ";
	result  = test;
	result -= test;
	cout << result << endl;
    //
	cout << "  GVector -= 2.0: ";	
	result  = test;
	result -= 2.0;
	cout << result << endl;
    //
	cout << "  GVector *= 2.0: ";	
	result  = test;
	result *= 2.0;
	cout << result << endl;
    //
	cout << "  GVector /= 2.0: ";	
	result  = test;
	result /= 2.0;
	cout << result << endl;
    //
	cout << "  GVector = -GVector: ";	
	result = -test;
	cout << result << endl;
    //
    cout << "  Devide by zero: ";
	result  = test;
	result /= 0.0;
	cout << result << endl;
    //
    cout << "  GVector + GVector: ";
	result = test + test;
	cout << result << endl;
    //
    cout << "  GVector + 2.0: ";
	result = test + 2.0;
	cout << result << endl;
    //
    cout << "  GVector + 2: ";
	result = test + 2;
	cout << result << endl;
    //
    cout << "  2.0 + GVector: ";
	result = 2.0 + test;
	cout << result << endl;
    //
    cout << "  2 + GVector: ";
	result = 2 + test;
	cout << result << endl;
    //
    cout << "  GVector - GVector: ";
	result = test - test;
	cout << result << endl;
    //
    cout << "  GVector - 2.0: ";
	result = test - 2.0;
	cout << result << endl;
    //
    cout << "  2.0 - GVector: ";
	result = 2.0 - test;
	cout << result << endl;
    //
	cout << "  Scalar (or dot) product GVector * GVector: ";	
	cout << (test * test) << endl;
    //
    cout << "  GVector * 2.0: ";
	result = test * 2.0;
	cout << result << endl;
    //
    cout << "  2.0 * GVector: ";
	result = 2.0 * test;
	cout << result << endl;
    //
    cout << "  |GVector| (vector norm): ";
	cout << norm(test) << endl;
    //
    cout << "  min(GVector): ";
	cout << min(test) << endl;
    //
    cout << "  max(GVector): ";
	cout << max(test) << endl;
    //
    cout << "  sum(GVector): ";
	cout << sum(test) << endl;
	//
	cout << "  GVector: ";
	cout << test << endl;
    //
    cout << "  acos(GVector/10.0): ";
	cout << acos(test/10.0) << endl;
    //
    cout << "  acosh(GVector): ";
	cout << acosh(test) << endl;
    //
    cout << "  asin(GVector/10.0): ";
	cout << asin(test/10.0) << endl;
    //
    cout << "  asinh(GVector/10.0): ";
	cout << asinh(test/10.0) << endl;
    //
    cout << "  atan(GVector/10.0): ";
	cout << atan(test/10.0) << endl;
    //
    cout << "  atanh(GVector/10.0): ";
	cout << atanh(test/10.0) << endl;
    //
    cout << "  cos(GVector): ";
	cout << cos(test) << endl;
    //
    cout << "  cosh(GVector): ";
	cout << cosh(test) << endl;
    //
    cout << "  exp(GVector): ";
	cout << exp(test) << endl;
    //
    cout << "  fabs(cos(GVector)): ";
	cout << fabs(cos(test)) << endl;
    //
    cout << "  log(GVector): ";
	cout << log(test) << endl;
    //
    cout << "  log10(GVector): ";
	cout << log10(test) << endl;
    //
    cout << "  sin(GVector): ";
	cout << sin(test) << endl;
    //
    cout << "  sinh(GVector): ";
	cout << sinh(test) << endl;
    //
    cout << "  sqrt(GVector): ";
	cout << sqrt(test) << endl;
    //
    cout << "  tan(GVector): ";
	cout << tan(test) << endl;
    //
    cout << "  tanh(GVector): ";
	cout << tanh(test) << endl;
    // 
    cout << "  Incompatible size GVector + GVector: ";
	result = test + bigger;
	cout << result << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }
  try {
    cout << "  cross(a,b) (using 5-dim vectors): ";
	cout << cross(test,test) << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }
  try {
    cout << "  cross(a,b) (using vectors with different dimension): ";
	cout << cross(test,bigger) << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }
  try {
    cout << "  cross(a,b) (using 3-dim vectors): " << endl;
	GVector test_cross_a(3);
	GVector test_cross_b(3);
	test_cross_a(0) = 1.0;
	test_cross_b(1) = 1.0;
	cout << "   a:     " << test_cross_a << endl;
	cout << "   b:     " << test_cross_b << endl;
	cout << "   cross: " << cross(test_cross_a, test_cross_b) << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }

  // Test 7: Vector comparison
  try {
    cout << "Test 7: Comparison:" << endl;
	//
	cout << "  GVector == GVector: ";
	cout << (test == test) << endl;
	//
	cout << "  GVector == GVector(0): ";
	GVector test7(num);
	cout << (test == test7) << endl;
	//
	cout << "  GVector == GVector (bigger): ";
	cout << (test == bigger) << endl;
	//
	cout << "  GVector != GVector: ";
	cout << (test != test) << endl;
	//
	cout << "  GVector != GVector(0): ";
	cout << (test != test7) << endl;
	//
	cout << "  GVector != GVector (bigger): ";
	cout << (test != bigger) << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }


  // Return
  return 0;
}
