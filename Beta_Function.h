/*
  File:         Beta_Function.h
  Version:      0.0.1
  Date:         30-Jan-2019
  Revision:     01-Feb-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Credits:	https://en.wikipedia.org/wiki/Beta_function
		http://www.mymathlib.com/c_source/functions/gamma_beta/beta_function.c

  Beta_Function.h - Library for 'duino
  https://github.com/newEndeavour/Beta_Function

  Beta_Function implements the mathemical Beta Function.  
  
  Copyright (c) 2018-2019 Jerome Drouin  All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Editions:	
	0.0.1	: first version
	
*/


#ifndef BETA_FUNCTION_H
#define BETA_FUNCTION_H

double 			Beta_Function(double a, double b);
long double 		xBeta_Function(long double a, long double b);

double 			Ln_Beta_Function(double a, double b);
long double 		xLn_Beta_Function(long double a, long double b);

double 			Beta_Distribution_Func(double x, double a, double b);
static long double 	xBeta_Distribution_Func(double x, double a, double b);
static long double 	Beta_Continued_Fraction( long double x, long double a, long double b);

double 			Incomplete_Beta_Function(double x, double a,double b);
long double 		xIncomplete_Beta_Function(long double x, long double a, long double b);



#endif
