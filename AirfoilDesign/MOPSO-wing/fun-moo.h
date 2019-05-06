/***************************************************************************
                          fun-moo.h  -  description
                             -------------------
    begin                : Tue May 15 2001
    copyright            : (C) 2001 by Max Salazar
    email                : max@maxnet.cc
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
using namespace std;

 /////////////////////////////////////////////////////////////////////////////////////// 1
double f1_1_mo(double *x) {
	double r;
	
	if (x[0] <= 1.0)
		r = -x[0];
	if ((1.0 < x[0]) && (x[0] <= 3.0))
		r = -2.0 + x[0];
	if ((3.0 < x[0]) && (x[0] <= 4.0))
		r = 4.0 - x[0];
	if (x[0] > 4.0)
		r = -4.0 + x[0];
		
	return r;
}

double f2_1_mo(double *x) {
	return pow(x[0] - 5.0, 2);
}

/////////////////////////////////////////////////////////////////////////////////////// 2
double f1_2_mo(double *x) {
	double r = 0.0;
	unsigned int i;
		
	for(i = 0; i < 2; i++)
		r += -10.0 * exp(-0.2 * sqrt(pow(x[i], 2) + pow(x[i + 1], 2)) );
	
	return r;
}

double f2_2_mo(double *x) {
	double r = 0.0;
	unsigned int i;
		
	for(i = 0; i < 3; i++)
		r += pow(fabs(x[i]), 0.8) + 5.0 * sin(pow(x[i], 3));
	return r;
}

/////////////////////////////////////////////////////////////////////////////////////// 3
double f1_3_mo(double *x) {
	double r = 0.0;
	unsigned int i;
		
	for(i = 0; i < 3; i++)
		r += pow(x[i] - (1.0 / sqrt(2)), 2);
	
	return 1.0 - exp(-r);
}

double f2_3_mo(double *x) {
	double r = 0.0;
	unsigned int i;
		
	for(i = 0; i < 3; i++)
		r += pow(x[i] + (1.0 / sqrt(2)), 2);
	
	return 1.0 - exp(-r);
}

/////////////////////////////////////////////////////////////////////////////////////// 4
double f1_4_mo(double *x) {
	return pow( x[0], 2);
}

double f2_4_mo(double *x) {
	return pow( x[0] - 2.0, 2);
}

/////////////////////////////////////////////////////////////////////////////////////// 5
double f1_5_mo(double *x) {
	double A1, A2, B1, B2;
	
	A1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0);
	A2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0);
	B1 = 0.5 * sin(x[0]) - 2.0 * cos(x[0]) + sin(x[1]) - 1.5 * cos(x[1]);
	B2 = 1.5 * sin(x[0]) - cos(x[0]) + 2.0 * sin(x[1]) - 0.5 * cos(x[1]);
	
	return 1.0 * -(1.0 + pow(A1 - B1, 2) + pow(A2 - B2, 2) );
}

double f2_5_mo(double *x) {
	return 1.0 * -(pow(x[0] + 3.0, 2) + pow(x[1] + 1.0, 2) );
}

/////////////////////////////////////////////////////////////////////////////////////// 6
double f1_6_mo(double *x) {
	return 0.5 * (pow(x[0], 2) + pow(x[1], 2)) + sin(pow(x[0], 2) + pow(x[1], 2));
}

double f2_6_mo(double *x) {
	return (pow(3.0 * x[0] - 2.0 * x[1] + 4.0, 2) / 8.0) + (pow(x[0] - x[1] + 1.0, 2) / 27.0) + 15.0;
}

double f3_6_mo(double *x) {
	return (1 / (pow(x[0], 2) + pow(x[1], 2) + 1.0)) - 1.1 * exp(-pow(x[0], 2) - pow(x[1], 2));
}

/////////////////////////////////////////////////////////////////////////////////////// 7
double f1_7_mo(double *x) {
	return x[0];
}

double f2_7_mo(double *x) {
double a, q;
	a = 2.0;
	q = 4.0;
	
	return (1.0 + 10.0 * x[1]) * (1.0 - pow((x[0] / (1.0 + 10.0 * x[1])), a) - (x[0] / (1.0 + 10.0 * x[1])) * sin(2.0 * PHI * q * x[0]));
}

/////////////////////////////////////////////////////////////////////////////////////// 8
double f1_8_mo(double *x) {
	return (pow(x[0] - 2.0, 2) / 2.0) + (pow(x[1] + 1.0, 2) / 13.0) + 3.0;
}

double f2_8_mo(double *x) {
	return (pow(x[0] + x[1] - 3.0, 2) / 36.0) + (pow(-x[0] + x[1] + 2.0, 2) / 8.0) - 17.0;
}

double f3_8_mo(double *x) {
	return (pow(x[0] + 2 * x[1] - 1.0, 2) / 175.0) + (pow(2.0 * x[1] - x[0], 2) / 17.0) - 13.0;
}

/////////////////////////////////////////////////////////////////////////////////////// 9
double f1_9_mo(double *x) {
	return x[0];
}

double f2_9_mo(double *x) {
	double g, h;
	g = 11.0 + pow(x[1], 2) - 10.0 * cos(2.0 * PHI * x[1]);
	
	if (x[0] <= g)
		h = 1.0 - sqrt(x[0]/g);
	else
	 h=0.0;
	
	 return g*h;
}
/////////////////////////////////////////////////////////////////////////////////////// 10
double f1_10_mo(double *x) {
	return x[0];
}

double f2_10_mo(double *x) {
	double g;
	g = 2.0 - exp(-pow((x[1] - 0.2) / 0.004, 2)) - 0.8 * exp(-pow((x[1] - 0.6) / 0.4, 2));
	
	 return g/x[0];
}


///////////////////////////////////////////////////////////////////////////////////////
double f1_10100_mo(double *x) {
	  return((-x[0]*x[0]+x[1]));
}

double f2_10100_mo(double *x) {
	return((x[0]/2.0 + x[1] + 1));
}

double f1_10200_mo(double *x) {
	  return (pow(x[0]-2,2) + pow(x[1] - 1,2) + 2);
}

double f2_10200_mo(double *x) {
	return (9*x[0] - pow(x[1] - 1,2));
}


double f1_10300_mo(double *x) {// four bar plane truss
  double F = 10;
  double E = 2*pow(10,5);
  double L = 200;

  return( L*(2*x[0]+sqrt(2)*x[1]+sqrt(x[2])+x[3]) );
}

double f2_10300_mo(double *x) {
  double F = 10;
  double E = 2*pow(10,5);
  double L = 200;

  return( (F*L)/E * ( 2/x[0] + (2*sqrt(2))/x[1] - (2*sqrt(2))/x[2] + 2/x[3]) );


}
double f1_10400_mo(double *x) {//deb2 multiobjective genetic algorithms: problem dificulties and construction of test problems
  return (x[0]);

}

double f2_10400_mo(double *x) {//deb2 multiobjective genetic algorithms: problem dificulties and construction of test problems
  double g=2.0-exp(-pow(((x[1]-0.2)/0.004),2))-0.8*exp(-pow(((x[1]-0.6)/0.4),2));
  return ((double)g/x[0]);

}

double f1_10500_mo(double *x) {//deb
  return (x[0]);
}

double f2_10500_mo(double *x) {
  double f=x[0];
  double h;
  double g=11+pow(x[1],2)-10*cos(2*3.1415926*x[1]);
  if(f<=g)
    h=1-pow((f/g),(double)(1.0/2.0));
  else
    h=0;

  return(g*h);

 
}
////////////////////////////////////////////////////////////////////////////////
// ZDT3-10D
////////////////////////////////////////////////////////////////////////////////
double f1_ZDT3_10D(double *x) {
	return x[0];
}

double f2_ZDT3_10D(double *x) {
	float r=0, PI = 3.14159265;
	int D = 10;
	for (int i = 1; i < D; i++) r +=  x[i];
	return (1+9*r/(D-1)) * (1-pow(x[0]/(1+9*r/(D-1)),(float)0.5) - 
		(x[0]/(1+9*r/(D-1))) * sin(10*PI*x[0]));
}

void ZDT3_10D(double *fitness, double *x) {
	float r=0, PI = 3.14159265;
	int D = 10;
	fitness[0] = x[0];
	for (int i = 1; i < D; i++) r +=  x[i];
	fitness[1] = (1+9*r/(D-1)) * (1-pow(x[0]/(1+9*r/(D-1)),(float)0.5) - 
		(x[0]/(1+9*r/(D-1))) * sin(10*PI*x[0]));
}


///////////////////////////////////////////////////////////////////////////////////////
/*
double f1__mo(double *x) {
	return x[0];
}

double f2__mo(double *x) {
	double r = 0.0;
		
		r = 2.0 - exp(pow((x[1] - 0.2) / 0.004, 2)) - 0.8 * exp(pow((x[1] - 0.6) / 0.4, 2));
	
	return r / x[0];
}
*/
