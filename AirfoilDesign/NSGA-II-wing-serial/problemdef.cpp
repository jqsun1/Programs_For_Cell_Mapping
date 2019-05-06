/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <cstring>

# include "global.h"
# include "rand.h"
#include "DTLZ.h"
#include "ExampleProblems.h"
#include "TransFunctions.h"
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
#include "SolverWing.h"



extern char    strTestInstance[256];
extern int nreal;
extern int nobj;




extern"C" {
	// everything is passed by pointer
	void mysub_(double * OBJS, long int * NOBJ, double * VARS, long int * NVAR);
}

void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
	vector<double> x_var(xreal, xreal+nreal);
	vector<double> y_obj(obj, obj+nobj);
	if(!strcmp(strTestInstance,"DTLZ1"))
	{
		y_obj = DTLZ1(x_var, y_obj);
	}
	if(!strcmp(strTestInstance,"DTLZ2"))
	{
		y_obj = DTLZ2(x_var, y_obj);
	}
	if(!strcmp(strTestInstance,"DTLZ3"))
	{
		y_obj = DTLZ3(x_var, y_obj);
	}
	if(!strcmp(strTestInstance,"DTLZ4"))
	{
		y_obj = DTLZ4(x_var, y_obj);
	}
	if(!strcmp(strTestInstance,"wfg1"))
	{
		y_obj = Problems::WFG1(x_var, 4, y_obj.size() );
	}
	if(!strcmp(strTestInstance,"wfg2"))
	{
		y_obj = Problems::WFG2(x_var, 4, y_obj.size() );
	}
	if(!strcmp(strTestInstance,"wfg3"))
	{
		y_obj = Problems::WFG3(x_var, 4, y_obj.size() );
	}
	if(!strcmp(strTestInstance,"wfg4"))
	{
		y_obj = Problems::WFG4(x_var, 4, y_obj.size() );
	}
	if(!strcmp(strTestInstance,"ZDT1"))
	{
	    double f1, f2, g, h;
		int i;
		f1 = x_var[0];
		g = 0.0;
		for (i=1; i<x_var.size(); i++)
		{
		    g += x_var[i];
		}
		g = 9.0*g/(x_var.size()-1);
		g += 1.0;
		h = 1.0 - sqrt(f1/g);
		f2 = g*h;
		y_obj[0] = f1;
		y_obj[1] = f2;
	}
	if(!strcmp(strTestInstance,"ZDT2"))
	{
		double g = 0;
		for(int n=1;n<x_var.size();n++)
			g+= xreal[n];
		g = 1 + 9*g/(x_var.size()-1);
		y_obj[0] = xreal[0];
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));
	}
	if(!strcmp(strTestInstance,"ZDT3"))
	{
		double g = 0;
		for(int n=1;n<x_var.size();n++)
			g+= xreal[n];
		g = 1 + 9*g/(x_var.size()-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(x_var[0]/g) - x_var[0]*sin(10*PI*x_var[0])/g);

	}
	if(!strcmp(strTestInstance,"ZDT4"))
	{
		double g = 0;
		for(int n=1;n<x_var.size();n++)
		{
			double x = 10*(x_var[n] - 0.5);
			g+= x*x - 10*cos(4*PI*x);
		}
		g = 1 + 10*(x_var.size()-1) + g;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));
	}
	if(!strcmp(strTestInstance,"ZDT6"))
	{
		double g = 0;
		for(int n=1;n<x_var.size();n++)
			g+= x_var[n]/(x_var.size() - 1);
		g = 1 + 9*pow(g,0.25) ;

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*PI*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));
	}
	if(!strcmp(strTestInstance,"ZDT3FORTRAN"))
	{
		long int NOBJ = 2;
		long int NVAR = (long int)x_var.size();
		double OBJS [NOBJ];
		double VARS [NVAR];
		
		for(int n=0;n<x_var.size();n++)
			VARS[n]= xreal[n];
		mysub_(OBJS, &NOBJ, VARS, &NVAR);

		y_obj[0] = OBJS[0];
		y_obj[1] = OBJS[1];

	}
	if(!strcmp(strTestInstance,"airfoil"))
	{
		int nvar = x_var.size();
		float x[nvar];
		float y[nobj];
		for(int n=0;n<nvar;n++)
			x[n]= (float)xreal[n];
		f(y, x);
		for(int n=0;n<nobj;n++)
			y_obj[n]= (double)y[n];
	}
	
//	obj = &y_obj[0];
	for(int i=0; i<nobj; i++){
		obj[i] = y_obj[i];
	}
	return;
}


