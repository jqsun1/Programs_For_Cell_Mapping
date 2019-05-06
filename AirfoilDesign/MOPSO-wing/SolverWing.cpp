#include <iostream>
#include <fstream>      // std::ofstream
#include <iomanip>      // std::setw
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "SolverWing.h"










extern"C" {
	// everything is passed by pointer
	void xfoilfun_(bool * CONV, double * CL, double * CD, double * CM, double *XX, double *YY, double * Re, double * Alfa, long int * npointslong, long int * samples);
}
////////////////////////////////////////////////////////////////////////////////
void matLU(double *A, double *L, double *U, int dim){
	//LU decomposition of a square matrix using Doolittle algorithm
	//matrices are treated as 1D pointers here
	int i,j,k,r;
	double sum;
	//const int dim = npp1;

	//initialization
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			L[i*dim+j] = 0.f;
			U[i*dim+j] = 0.f;
		}
	}

	//assign the first row of U and first column of L
	for(j=0;j<dim;j++)
		U[j] = A[j];
	for(i=0;i<dim;i++)
		L[i*dim] = A[i*dim]/U[0];

	//assign the rest rows/cols
	for(k=1;k<dim;k++){
		//kth row of U
		for(j=k;j<dim;j++){
			sum = 0.f;
			for(r=0;r<=k-1;r++)
				sum+=L[k*dim+r]*U[r*dim+j];
			U[k*dim+j] = A[k*dim+j]-sum;
		}
		//kth col of L
		L[k*dim+k] = 1.f; //diagnoal elements are 1
		for(i=k+1;i<dim;i++){
			sum = 0.f;
			for(r=0;r<=k-1;r++)
				sum+=L[i*dim+r]*U[r*dim+k];
			L[i*dim+k] = (A[i*dim+k]-sum)/U[k*dim+k];
		}
	}
}
////////////////////////////////////////////////////////////////////////////////
void matsolve(double *A, double *b, double *x, double* L, double* U, double *y, int n){
	//solve linear system Ax=b using LU decomposition
	int k,r;
	//const int n = npp1;
	//double L[n][n],U[n][n],y[n];
	double sum;

	matLU(A,L,U,n); //LU decomposition

	//Ly=b
	for(k=0;k<n;k++){
		sum = 0.f;
		for(r=0;r<=k-1;r++)
			sum+=L[k*n+r]*y[r];
		y[k] = b[k]-sum;
	}
	//Ux=y
	for(k=n;k!=0;k--){
		sum = 0.f;
		for(r=k+1;r<=n;r++)
			sum+=U[(k-1)*n+r-1]*x[r-1];
		x[k-1] = (y[k-1]-sum)/U[(k-1)*n+k-1];
	}
}
////////////////////////////////////////////////////////////////////////////////


void parsec(double * a, double * param) {
	//   the function returns the PARSEC coefficients for for an airfoil using 
	//	   the 11 coefficients.
	//	   the arguments are:
	//	   Input: 
	//	       param : the vector of length 12 which includes:
	//	           Rleu Upper leading edge radius
	//	           Rlel   lower leading edge radius
	//	           Xup    position of upper crest
	//	           Zup    upper crest point
	//	           Zxxup  upper crest curvature
	//	           Xlo    position of lower crest
	//	           Zlo    lower crest point
	//	           Zxxlo  lower crest curvature
	//	           Tte    trailing edge thickness
	//	           Toff   trailing edge offset
	//	           Ate    trailing edge direction angle
	//	           Bte    trailing edge wedge angle
	//	   Output: 
	//	       a : 12 coefficients. the first 6 are for the upper surface and 
	//				the second 6 are for the lower surface 
	//	           curves of airfoil in the form:
	//	           Z = a1 * x^(1/2) + a2 * x^(3/2) + a3 * x^(5/2) +...
	//	               a4 * x^(7/2) + a5 * x^(9/2) + a6 * x^(11/2)
	//	   Author: Yousef Naranjani
	//	   Oct 05 2015
	double pi = 4.f * atan(1.f);
	double L[5][5], U[5][5], Y[5];
	// The 12 input parameters
	double Rleu = param[0];
	double Rlel = param[1];
	double Xup  = param[2];
	double Zup  = param[3];
	double Zxxup= param[4];
	double Xlo  = param[5];
	double Zlo  = param[6];
	double Zxxlo= param[7];
	double Tte  = param[8];
	double Toff = param[9];
	double Ate  = param[10] * pi / 180.f;
	double Bte  = param[11] * pi / 180.f;
	//	% for the upper surface:
	double A1up = pow(2*Rleu,.5f);
	double XupN15 = pow(Xup,-1.5f);
	double XupN5 = pow(Xup,-.5f);
	double Xup5 = pow(Xup,.5f);
	double Xup15 = pow(Xup,1.5f);
	double Xup25 = pow(Xup,2.5f);
	double Xup35 = pow(Xup,3.5f);
	double Xup45 = pow(Xup,4.5f);
	double Xup55 = pow(Xup,5.5f);
	
	double A[5][5] = {
		{1.f			,1.f			,1.f			,1.f			,1.f},
		{Xup15			,Xup25			,Xup35			,Xup45			,Xup55},
		{1.5f			,2.5f			,3.5f			,4.5f			,5.5f},
		{1.5f*Xup5		,2.5f*Xup15		,3.5f*Xup25		,4.5f*Xup35		,5.5f*Xup45},
		{3.f/4.f*XupN5	,15.f/4.f*Xup5	,35.f/4.f*Xup15	,64.f/3.f*Xup25	,99.f/4.f*Xup35}
		};

	double RHS[5] = {Toff+.5f*Tte		- A1up,
					Zup 				- A1up * Xup5,
					tan(-Ate-.5f*Bte) 	- A1up * .5f,
					0 					- A1up * .5f * XupN5,
					Zxxup 				+ A1up * .25f * XupN15};
	matsolve((double*)A,RHS,a+1,(double*)L,(double*)U,Y,5);
	a[0] = A1up;
	
	//	% for the lower surface: (Warning: we are reusing some of the variables!)
	A1up = -pow(2*Rlel,.5f);
	XupN15 = pow(Xlo,-1.5f);
	XupN5 = pow(Xlo,-.5f);
	Xup5 = pow(Xlo,.5f);
	Xup15 = pow(Xlo,1.5f);
	Xup25 = pow(Xlo,2.5f);
	Xup35 = pow(Xlo,3.5f);
	Xup45 = pow(Xlo,4.5f);
	Xup55 = pow(Xlo,5.5f);
	
	 double A2[5][5] = {
		{1.f			,1.f			,1.f			,1.f			,1.f},
		{Xup15			,Xup25			,Xup35			,Xup45			,Xup55},
		{1.5f			,2.5f			,3.5f			,4.5f			,5.5f},
		{1.5f*Xup5		,2.5f*Xup15		,3.5f*Xup25		,4.5f*Xup35		,5.5f*Xup45},
		{3.f/4.f*XupN5	,15.f/4.f*Xup5	,35.f/4.f*Xup15	,64.f/3.f*Xup25	,99.f/4.f*Xup35}
		};

	double RHS2[5] = {Toff-.5f*Tte		- A1up,
					Zlo 				- A1up * Xup5,
					tan(-Ate+.5f*Bte) 	- A1up * .5f,
					0 					- A1up * .5f * XupN5,
					Zxxlo 				+ A1up * .25f * XupN15};
	matsolve((double*)A2,RHS2,a+7,(double*)L,(double*)U,Y,5);
	a[6] = A1up;
	
}
////////////////////////////////////////////////////////////////////////////////
void parsec2pointsPrep(double * XX, double * coeff, int npoints) {
	// parsec2pointsPrep prepares X points and its half powers (so that Z values
	// can be evaluated using these points in another function)
	// Input:
	//		npoints	the number of points for each side of airfoil. recommended:80
	// Output:
	//		XX		the x points. XX[0] = 0, XX[npoints-1] = 1
	//				these points are more dense on the sides to capture the 
	//				change in the profile more accurately (refer to naca.f in Xfoil)
	//		coeff	the array of half powers of X.(Used late in the form Z=X.a)
	//				the format: 
	//				coeff =[XX[0]^(1/2)	XX[0]^(3/2) 	XX[0]^(5/2) XX[0]^(7/2)...
	//						XX[0]^(9/2)	XX[0]^(11/2)	XX[1]^(1/2) XX[1]^(3/2) ... 
	//						...			XX[npoints-1]^(9/2)	XX[npoints-1]^(11/2)];
	int npowers = 6; 			// total # of half powers used
	double AN = 1.5, ANP, frac; 	// distribution of points are done using the method in Xfoil
	ANP = AN + 1;
	for (int i = 0; i < npoints; i++) {
		frac = (double)i / (double) (npoints-1);
		XX[i] = 1.f - ANP*frac*pow(1.f-frac, AN) - pow(1.f-frac, ANP);
		if (i == npoints-1) XX[i] = 1;
		for (int j = 0; j < npowers; j++) {
			coeff[i*npowers+j] = pow(XX[i],(j*2.f+1.f)/2.f)	;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////
void parsec2points(double * X, double * Y, double * a, int & npoints, double * XX, double * coeff) {
	// this function creates the airfoil profile with npoints on each side using
	// the parsec parameters a and will return the coordinates X and Y
	//	Input:
	//			a		PARSEC coefficients of length 12. The first 6 are to be used
	//					for top sirface and a[6:11] to be used for lower surface
	//			npoints	The number of points to be generated on each side of airfoil
	//					recommended: 80
	//			coeff	the prepared array of half powers of X in: Z=X.a
	//					it is calculated separately because it is a repeated
	//					the format: 
	//					X = [x[0]^(1/2)	x[0]^(3/2) x[0]^(5/2) x[0]^(7/2)...
	//						x[0]^(9/2)	x[0]^(11/2)	x[1]^(1/2) x[1]^(3/2) ... 
	//						...			x[npoints-1]^(9/2)	x[npoints-1]^(11/2)];
	//	Outpt:
	//			X,Y		X,Y coordinated of the arfoil surface sarting from the 
	//					tail (1.0,0) and doing an CCW turn to end at tail. The
	//					length of X,Y are 2*npoints.
	//	Author: Yousef Naranjani
	//	Oct 06 2015
	int npowers = 6; 			// total # of half powers used
	int i, cnt = -1;
	for (i = npoints-1; i>=0; i--) { //Upper surface
		cnt++;
		X[cnt] = XX[i];
		Y[cnt] = coeff[i*npowers+0] * a[0] +
				coeff[i*npowers+1] * a[1] +
				coeff[i*npowers+2] * a[2] +
				coeff[i*npowers+3] * a[3] +
				coeff[i*npowers+4] * a[4] +
				coeff[i*npowers+5] * a[5]; 
		//printf("%4d,  %10.5f,  %10.5f\n",cnt,X[cnt],Y[cnt]);
	}
	//printf("The lower surface:\n");
	for (i = 1; i<npoints; i++) {	 //Lower surface
		cnt++;
		X[cnt] = XX[i];
		Y[cnt] = coeff[i*npowers+0] * a[6] +
				coeff[i*npowers+1] * a[7] +
				coeff[i*npowers+2] * a[8] +
				coeff[i*npowers+3] * a[9] +
				coeff[i*npowers+4] * a[10] +
				coeff[i*npowers+5] * a[11];
		//printf("%4d,  %10.5f,  %10.5f\n", cnt, X[cnt],Y[cnt]); 
	}
}
////////////////////////////////////////////////////////////////////////////////
bool fileRead(double * xtmp, double * ytmp, int & i, std::string & airfoil, std::string filename) {
	double x=-1,y=-1;
	std::ifstream myfile(filename.c_str());
	getline(myfile, airfoil);
	//cout << "TheNAME:" <<airfoil << endl;
	i = 0;
	while (myfile >> x >> y) {
		xtmp[i] = x; 
		ytmp[i] = y;
		i++;
		//cout << i << ", "<< x << ", " << y << endl;
	}
	if (i == 0) return true;
	else		return false;
	
}
////////////////////////////////////////////////////////////////////////////////
void f(double * y, double * x) {
	// This solver recieved the PARSEC parameters as input then transform them
	// into coordinates of airfoil and pass them to the XFOIL solver for the 
	// aerodynamic properties to be calculated. The output of the function 
	// include all the opjectives. The list of arguments are:
	// Input:
	// 			x	PARSEC parameters array
	// Output:
	// 			y 	array of objectives including
	//				- Cl / Cd 
	//				- dCl / dAlfa 
	//				Cm^2		
	double a[12], ClCd = 0;
	int npoints = 80;
	int i, prevInd=-1;//prevInd=he index of the previous successful converged case
	long int samples = 60, length = (long int)2*npoints-1;
	bool CONV[samples];
	double X[length], Y[length];
	double Re = 0.5e6, Alfa[samples], CL[samples], CD[samples], CM[samples];
	y[0] = y[1] = y[2] = inf;
	parsec(a, x);
	//printf("x=\n");
	//for (i=0;i<D;i++) printf("%8.5f, ",x[i]);
	//printf("\n");
	//printf("a=\n");
	//for (i=0;i<12;i++) printf("%8.5f, ",a[i]);
	//printf("\n");
	parsec2points( X, Y, a, npoints, XX, coeff);
	y[0] = y[1] = y[2] = inf;
	Alfa[0] = CL[0] = CD[0] = CM[0] = 0.f;
	CONV[0] = false;
	for (i = 1; i < samples; i++) {
		Alfa[i] = Alfa[i-1] + 0.25;
		CONV[i] = false;
		CL[i] = CD[i] = CM[i] = 0;
	}
	
	xfoilfun_(CONV, CL, CD, CM, X, Y, &Re, Alfa, &length, &samples);
	//
	for (i = 0; i < samples; i++) {
		//printf("CONV=%d,A=%5.2f,CL=%8.6f,CD=%8.6f,CL/CD=%8.6f,CM=%8.6f\n",CONV[i],Alfa[i],CL[i],CD[i],CL[i]/CD[i],CM[i]);
		if (CONV[i]==true && CD[i] != 0 && CL[i]/CD[i] > ClCd) {
			ClCd = CL[i] / CD[i];
			y[0] = -ClCd;
			if (prevInd != -1) {
				y[1] = - (CL[i]-CL[prevInd]) / (Alfa[i]-Alfa[prevInd]);
			} 
			else if (Alfa[i]==0) {
				y[1] = 0;
			}
			else {
				y[1] = - (CL[i]) / (Alfa[i]);
			}
			y[2] = powf(CM[i],2.f);
			//printf(" objectives :%f,%f,%f\n",y[0],y[1],y[2]);
			prevInd = i;
		}
		else if (CONV[i]==true) {
			prevInd = i;
		}
	}
	//printf("%8.6f,%8.6f,%8.6f\n",y[0],y[1],y[2]);
}


