#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#include "SolverWing.h"
float coeff[pps*6];
float XX[pps];

////////////////////////////////////////////////////////////////////////////////
int main(void) {
	float param[12] = {0.015025628551640,
   			0.006616870933959,
			0.535701875370922,
			0.038128511401401,
			-1.998369922315678,
			0.413720925926081,
			-0.058963126687472,
			1.631935344288543,
			0.007495015102874,
			0.002270498904051,
			4.717485467228419,
			12.423983727690421};
	float objs[3];
//	float a[12];
//	int npoints = 80;
//	long int npointslong;
//	float coeff[npoints*6], XX[npoints], X[2*npoints-1], Y[2*npoints-1];
//	float Alfa = 5.f, Re = 500000.f;
//	parsec2pointsPrep( XX, coeff, npoints);
	
//	parsec(a, param);
//	printf("a=\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n%10.6f,\n",
//		a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11]);
	
//	parsec2points( X, Y, a, npoints, XX, coeff);
//	npointslong = (long int) (2*npoints-1);
//	xfoilfun_(X, Y, &Re, &Alfa, &npointslong);
	parsec2pointsPrep( XX, coeff, pps);
	f(objs, param);
	printf("The objectives are:%f,%f,%f\n",objs[0],objs[1],objs[2]);
}
