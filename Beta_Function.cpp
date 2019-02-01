/*
  File:         Beta_Function.cpp
  Version:      0.0.1
  Date:         27-Jan-2019
  Revision:     01-Feb-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:	Please go to Beta_Function for Edition Notes.

  Credits:	https://en.wikipedia.org/wiki/Beta_Function
		http://www.mymathlib.com/functions/gamma_beta.html

  Beta_Function - Library for 'duino
  https://github.com/newEndeavour/Beta_Function

  Beta_Function implements the Beta function. 

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
  along double with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

// Includes
#include <Arduino.h>   // required for Serial.print()
#include <math.h>      // required for powl(), sinl(), fabsl() and ldexpl().
#include <float.h>     // required for DBL_MAX and LDBL_MAX and LDBL_EPSILON
#include <Beta_Function.h>
#include <Gamma_Function.h>


// Library definitions
double 			Beta_Function(double a, double b);
long double 		xBeta_Function(long double a, long double b);

double 			Ln_Beta_Function(double a, double b);
long double 		xLn_Beta_Function(long double a, long double b);

double 			Beta_Distribution_Func(double x, double a, double b);
static long double 	xBeta_Distribution_Func(double x, double a, double b);
static long double 	Beta_Continued_Fraction( long double x, long double a, long double b);

double 			Incomplete_Beta_Function(double x, double a,double b);
long double 		xIncomplete_Beta_Function(long double x, long double a, long double b);


// Constants
static const long double ln_LDBL_MAX =  1.13565234062941435e+4L;



// Active /////////////////////////////////////////////////////////////////
//Beta_Function
double Beta_Function(double a, double b)
{
long double beta = xBeta_Function( (long double) a, (long double) b);

	return (beta < DBL_MAX) ? (double) beta : DBL_MAX;
}


//xBeta_Function
long double xBeta_Function(long double a, long double b)
{
long double lnbeta;

	//If (a + b) <= Gamma_Function_Max_Arg() then simply return
     	//gamma(a)*gamma(b) / gamma(a+b).                         

   	if ( (a + b) <= Gamma_Function_Max_Arg() )
      		return xGamma_Function(a) / (xGamma_Function(a + b) / xGamma_Function(b));

	//If (a + b) > Gamma_Function_Max_Arg() then simply return
     	//exp(lngamma(a) + lngamma(b) - lngamma(a+b) ).           
	lnbeta = xLn_Gamma_Function(a) + xLn_Gamma_Function(b)
                                                 - xLn_Gamma_Function(a + b);
   
	return (lnbeta > ln_LDBL_MAX) ? (long double) LDBL_MAX : expl(lnbeta);

}


//Ln_Beta_Function
double Ln_Beta_Function(double a, double b)
{
	return (double) xLn_Beta_Function( (long double) a, (long double) b );
}


//xLn_Beta_Function
long double xLn_Beta_Function(long double a, long double b)
{
	// If (a + b) <= Gamma_Function_Max_Arg() return
     	//  log(gamma(a)*gamma(b) / gamma(a+b))
	if ( (a + b) <= (long double) Gamma_Function_Max_Arg() )
      		if ( a == 1.0L && b == 1.0L ) 
			return 0.0L;
      		else 
			return logl( xGamma_Function(a) /
                             	( xGamma_Function(a + b) / xGamma_Function(b) ));

     	// If (a + b) > Gamma_Function_Max_Arg() return
     	// lngamma(a) + lngamma(b) - lngamma(a+b)
	return xLn_Gamma_Function(a) + xLn_Gamma_Function(b)
				     - xLn_Gamma_Function(a+b);
}


//Beta_Distribution_Func
double Beta_Distribution_Func(double x, double a, double b)
{
	if ( x <= 0.0 ) 
		return 0.0;

   	if ( x >= 1.0 ) 
		return 1.0;

	return (double) xBeta_Distribution_Func(x, a, b);

}


//xBeta_Distribution_Func
static long double xBeta_Distribution_Func(double xx, double aa, double bb) 
{
long double x = (long double) xx;
long double a = (long double) aa;
long double b = (long double) bb;

	// Both shape parameters are strictly greater than 1.

   	if ( aa > 1.0 && bb > 1.0 )
      		if ( x <= (a - 1.0L) / ( a + b - 2.0L ) )
         		return Beta_Continued_Fraction(x, a, b);
      		else
         		return 1.0L - Beta_Continued_Fraction( 1.0L - x, b, a );
  
        // Both shape parameters are strictly less than 1. 
   	if ( aa < 1.0 && bb < 1.0 )  
      		return (a * xBeta_Distribution_Func(xx, aa + 1.0, bb) 
                      		+ b * xBeta_Distribution_Func(xx, aa, bb + 1.0) ) / (a + b); 
   
        // One of the shape parameters exactly equals 1.
   	if ( aa == 1.0 )
      		return 1.0L - powl(1.0L - x, b) / ( b * xBeta_Function(a,b) );

   	if ( bb == 1.0 ) 
		return powl(x, a) / ( a * xBeta_Function(a,b) );

      	// Exactly one of the shape parameters is strictly less than 1. 
	if ( aa < 1.0 )  
      		return xBeta_Distribution_Func(xx, aa + 1.0, bb)
            		+ powl(x, a) * powl(1.0L - x, b) / ( a * xBeta_Function(a,b) );
 
        // The remaining condition is b < 1.0 */
   	return xBeta_Distribution_Func(xx, aa, bb + 1.0)
        	    - powl(x, a) * powl(1.0L - x, b) / ( b * xBeta_Function(a,b) );

}


//Beta_Continued_Fraction
static long double Beta_Continued_Fraction( long double x, long double a,
                                                                 long double b)
{
long double Am1 = 1.0L;
long double A0 = 0.0L;
long double Bm1 = 0.0L;
long double B_0 = 1.0L;		// WARNING !!! : Renamed B0 -> B_0: Possible name Conflict with Variable define BO in binary.h (Arduino.h) //
long double e = 1.0L;
long double Ap1 = A0 + e * Am1;
long double Bp1 = B_0 + e * Bm1;
long double f_less = Ap1 / Bp1;
long double f_greater = 0.0L;
long double aj = a;
long double am = a;
static long double eps = 10.0L * LDBL_EPSILON;
int j = 0;
int m = 0;
int k = 1;


	if ( x == 0.0L ) 
		return 0.0L;
   
	while ( (2.0L * fabsl(f_greater - f_less) > eps * fabsl(f_greater + f_less)) ) {
      		Am1 	= A0;
      		A0 	= Ap1;
      		Bm1 	= B_0;
      		B_0 	= Bp1;
      		am 	= a + m;
      		e 	= - am * (am + b) * x / ( (aj + 1.0L) * aj );
      		Ap1 	= A0 + e * Am1;
      		Bp1 	= B_0 + e * Bm1;
      		k 	= (k + 1) & 3;
      		
		if (k == 1) 
			f_less = Ap1/Bp1;
      		else 
			if (k == 3) 
				f_greater = Ap1/Bp1;
      
		if ( fabsl(Bp1) > 1.0L) {
         		Am1 	= A0 / Bp1;
         		A0 	= Ap1 / Bp1;
         		Bm1 	= B_0 / Bp1;
         		B_0 	= 1.0;
      		} else {
         		Am1 	= A0;
         		A0 	= Ap1;
         		Bm1 	= B_0;
         		B_0 	= Bp1;
      		}
      
		m++;
      		j 	+= 2;
      		aj 	= a + j;
      		e 	= m * ( b - m ) * x / ( ( aj - 1.0L) * aj  );
      		Ap1 	= A0 + e * Am1;
      		Bp1 	= B_0 + e * Bm1;
      		k 	= (k + 1) & 3;
      		
		if (k == 1) 
			f_less = Ap1/Bp1;
      		else 
			if (k == 3) 
				f_greater = Ap1/Bp1;
   	}
	
   	return expl( a * logl(x) + b * logl(1.0L - x) + logl(Ap1 / Bp1) ) /
                                                ( a * xBeta_Function(a,b) );

}


//Incomplete_Beta_Function
double Incomplete_Beta_Function(double x, double a, double b)
{
long double beta = xIncomplete_Beta_Function( (long double) x, (long double) a, (long double) b);

	return (beta > DBL_MAX) ? DBL_MAX : (double) beta; 

}


//xIncomplete_Beta_Function
long double xIncomplete_Beta_Function(long double x, long double a, long double b) 
{

	// Both shape parameters are strictly greater than 1.
   	if ( a > 1.0L && b > 1.0L )
      	if ( x <= (a - 1.0L) / ( a + b - 2.0L ) )
        	return Beta_Continued_Fraction(x, a, b);
      	else
        	return xBeta_Function( (double) a, (double) b )
                 	                  - Beta_Continued_Fraction( 1.0L - x, b, a );
  
        // Both shape parameters are strictly less than 1.
    	if ( a < 1.0L && b < 1.0L )  
      		return xIncomplete_Beta_Function(x, a + 1.0L, b) 
                	                  + xIncomplete_Beta_Function(x, a, b + 1.0L); 
   
        // One of the shape parameters exactly equals 1.
   	if ( a == 1.0L ) 
		return ( 1.0L - powl(1.0L - x, b) ) / b;

	if ( b == 1.0L ) 
		return powl(x, a) / a;
      	
	// Exactly one of the shape parameters is strictly less than 1.
   	if ( a < 1.0L )  
      		return ( (a + b) * xIncomplete_Beta_Function(x, a + 1.0L, b)
                	                       + powl(x, a) * powl(1.0L - x, b) ) / a;
 
        // The remaining condition is b < 1.
   	return ( (a + b) * xIncomplete_Beta_Function(x, a, b + 1.0L)
        	                               - powl(x, a) * powl(1.0L - x, b) ) / b;

}


// Inactive /////////////////////////////////////////////////////////////////




//{Endoffile}


