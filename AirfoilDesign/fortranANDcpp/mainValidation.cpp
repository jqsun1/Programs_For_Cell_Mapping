#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
using namespace std;
#include "SolverWing.h"

float coeff[pps*6];
float XX[pps];

int main(void) {
	std::string filename = "./input/naca1408.dat";
	std::string airfoil;
	int tmplen = 900, i, length;
	long int len=0, samples = 48;
	bool err;
	float xtmp[tmplen], ytmp[tmplen];
	float Re = 0.5e6, Alfa[samples], CL[samples], CD[samples], CM[samples];
	bool CONV[samples];
	
	///////////////
	//coeff = new float [pps*6];
	//XX = new float [pps];
	parsec2pointsPrep( XX, coeff, pps);
	///////////////
	err = fileRead(xtmp, ytmp, length, airfoil, filename);
	if (err == true) {
		cout << "Error reading the file.\nAborting.\n";
		return 0;
	}
	float * x, * y;
	x = new float [len];
	y = new float [len];
	Alfa[0] = 0.f;
	// Alfa[1] = 2.f; Alfa[2] = 3.f; Alfa[3] = 4.f; Alfa[4] = 4.5f;	Alfa[5] = 5;
	CL[0] = CD[0] = CM[0] = 0;
	CONV[0] = true;
	//for (i = 0; i < 6; i++) {
	//	CONV[i] = true;
	//	CL[i] = CD[i] = CM[i] = 0;
	//}
	for (i = 1; i < 48; i++) {
		Alfa[i] = Alfa[i-1] + 0.25;
		CONV[i] = true;
		CL[i] = CD[i] = CM[i] = 0;
	}
	memcpy(x,xtmp,len * sizeof(float));
	memcpy(y,ytmp,len * sizeof(float));	
	len = (long int) length;
	cout << "Airfoil name: " << airfoil << endl;
	printf("Alfa		CL		CD		CM	\n");
	//for (Alfa = 8.75f; Alfa<9.f; Alfa=Alfa+0.25) {
		xfoilfun_(CONV, CL, CD, CM, xtmp, ytmp, &Re, Alfa, &len, &samples);
		//cout << "C++ Results: \nAlfa = " << Alfa << "\nCL = " << CL
		//	<< "\nCD = " << CD << "\nCM = " << CM << "\nCOnvergence: "
		//	<< CONV << endl;
		if (CONV) {
			printf(" \n");
		}
		else {
			
		}
	//}
	return 0;
}
