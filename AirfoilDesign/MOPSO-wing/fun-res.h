/***************************************************************************
                          function.h  -  description
                             -------------------
    begin                : Tue Feb 13 2001
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

double PHI = 3.141592653589793238462;

//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g01 IEEE		*
 *											*
 ************************/

double g01_IEEE(double *x) {
	unsigned int i;
	double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
	
	for(i = 0; i < 4; i++) {
		sum1 += x[i];
		sum2 += pow(x[i], 2);
	}
	
	for(i = 4; i < 13; i++)
		sum3 += x[i];
	
	return ((5.0 * sum1) - (5.0 * sum2) - sum3);
}

double g01_r1(double *x) {
	return ((2.0 * x[0]) + (2.0 * x[1]) + x[9] + x[10] - 10.0);
}

double g01_r2(double *x) {
	return ((2.0 * x[0]) + (2.0 * x[2]) + x[9] + x[11] - 10.0);
}

double g01_r3(double *x) {
	return ((2.0 * x[1]) + (2.0 * x[2]) + x[10] + x[11] - 10.0);
}

double g01_r4(double *x) {
	return (-(8.0 * x[0]) + x[9]);
}

double g01_r5(double *x) {
	return (-(8.0 * x[1]) + x[10]);
}

double g01_r6(double *x) {
	return (-(8.0 * x[2]) + x[11]);
}

double g01_r7(double *x) {
	return (-(2.0 * x[3]) - x[4] + x[9]);
}

double g01_r8(double *x) {
	return (-(2.0 * x[5]) - x[6] + x[10]);
}

double g01_r9(double *x) {
	return (-(2.0 * x[7]) - x[8] + x[11]);
}


//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g02 IEEE		*
 *											*
 ************************/

double g02_IEEE(double *x) {
	unsigned int i;
	double sum1 = 0.0, sum2 = 0.0, prod = 1.0;
	
	for(i = 0; i < 20; i++) {
		sum1 += pow(cos(x[i]), 4);
		sum2 += (double)(i + 1) * pow(x[i], 2);
		prod = prod * pow(cos(x[i]), 2);
	}
	
	return -(fabs((sum1 - (2.0 * prod)) / sqrt(sum2)));
}

double g02_r1(double *x) {
	unsigned int i;
	double prod = 1.0;
		
	for(i = 0; i < 20; i++)
		prod = prod * x[i];
		
	return 0.75 - prod;
}

double g02_r2(double *x) {
	unsigned int i;
	double sum = 0.0;
		
	for(i = 0; i < 20; i++)
		sum += x[i];
		
	return sum - (7.5 * 20.0);
}


//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g03 IEEE		*
 *											*
 ************************/

double g03_IEEE(double *x) {
	unsigned int i;
	double prod = 1.0;
	
	for(i = 0; i < 10; i++)
		prod *= x[i];
	
	return -(pow(sqrt(10), 10) * prod);
}

double g03_r1(double *x) {
	unsigned int i;
	double sum = 0.0;
	
	for(i = 0; i < 10; i++)
		sum += pow(x[i], 2);
	
	return sum - 1.0;
}


//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g04 IEEE		*
 *											*
 ************************/

double g04_IEEE(double *x) {
	return (5.3578547 * pow(x[2], 2)) + (0.8356891 * x[0] * x[4]) + (37.293239 * x[0]) - 40792.141;
}

double g04_r1(double *x) {
	return 85.334407 + (0.0056858 * x[1] * x[4]) + (0.0006262 * x[0] * x[3]) - (0.0022053 * x[2] * x[4]) - 92.0;
}

double g04_r2(double *x) {
	return -85.334407 - (0.0056858 * x[1] * x[4]) - (0.0006262 * x[0] * x[3]) + (0.0022053 * x[2] * x[4]);
}

double g04_r3(double *x) {
	return 80.51249 + (0.0071317 * x[1] * x[4]) + (0.0029955 * x[0] * x[1]) + (0.0021813 * pow(x[2], 2)) - 110.0;
}

double g04_r4(double *x) {
	return -80.51249 - (0.0071317 * x[1] * x[4]) - (0.0029955 * x[0] * x[1]) - (0.0021813 * pow(x[2], 2)) + 90.0;
}

double g04_r5(double *x) {
	return 9.300961 + (0.0047026 * x[2] * x[4]) + (0.0012547 * x[0] * x[2]) + (0.0019085 * x[2] * x[3]) - 25.0;
}

double g04_r6(double *x) {
	return -9.300961 - (0.0047026 * x[2] * x[4]) - (0.0012547 * x[0] * x[2]) - (0.0019085 * x[2] * x[3]) + 20.0;
}


//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g06 IEEE		*
 *											*
 ************************/

double g06_IEEE(double *x) {
	return pow(x[0] - 10.0, 3) + pow(x[1] - 20.0, 3);
}

double g06_r1(double *x) {
	return - pow(x[0] - 5.0, 2) - pow(x[1] - 5.0, 2) + 100.0;
}

double g06_r2(double *x) {
	return pow(x[0] - 6.0, 2) + pow(x[1] - 5.0, 2) - 82.81;
}

//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g07 IEEE		*
 *											*
 ************************/

double g07_IEEE(double *x) {
	return 	pow(x[0], 2) + pow(x[1], 2) + (x[0] * x[1]) - (14.0 * x[0]) - (16 * x[1]) + pow(x[2] - 10.0, 2) +
					(4 * pow(x[3] - 5.0, 2)) + pow(x[4] - 3.0, 2) +
					(2 * pow(x[5] - 1.0, 2)) + (5 * pow(x[6], 2)) + (7 * pow(x[7] - 11, 2)) + (2 * pow(x[8] - 10, 2)) +
					pow(x[9] -7, 2) + 45.0;
}

double g07_r1(double *x) {
	return -105.0 + (4.0 * x[0]) + (5.0 * x[1]) - (3.0 * x[6]) + (9.0 * x[7]);
}

double g07_r2(double *x) {
	return (10.0 * x[0]) - (8.0 * x[1]) - (17.0 * x[6]) + (2 * x[7]);
}

double g07_r3(double *x) {
	return (-8.0 * x[0]) + (2.0 * x[1]) + (5.0 * x[8]) - (2.0 * x[9]) - 12.0;
}

double g07_r4(double *x) {
	return (3.0 * pow(x[0] - 2.0, 2)) + (4.0 * pow(x[1] - 3.0, 2)) + (2.0 * pow(x[2], 2)) - (7.0 * x[3]) - 120.0;
}

double g07_r5(double *x) {
	return (5.0 * pow(x[0], 2)) + (8 * x[1]) + pow(x[2] - 6.0, 2) - (2.0 * x[3]) - 40.0;
}

double g07_r6(double *x) {
	return pow(x[0], 2) + (2.0 * pow(x[1] - 2.0, 2)) - (2.0 * x[0] * x[1]) + (14.0 * x[4]) - (6.0 * x[5]);
}

double g07_r7(double *x) {
	return (0.5 * pow(x[0] - 8.0, 2)) + (2.0 * pow(x[1] - 4.0, 2)) + (3.0 * pow(x[4], 2)) - x[5] - 30.0;
}

double g07_r8(double *x) {
	return (-3.0 * x[0])  + ( 6.0 * x[1]) + (12.0 * pow(x[8] - 8.0, 2)) - (7.0 * x[9]);
}


//////////////////////////////////////////////////////////

/************************
 *											*
 *	Problema g08 IEEE		*
 *											*
 ************************/

double g08_IEEE(double *x) {
	return (pow(sin(2.0 * PHI * x[0]), 3) * sin(2.0 * PHI * x[1])) / (pow(x[0], 3) * (x[0] + x[1]));
}

double g08_r1(double *x) {
	return pow(x[0], 2) - x[1] + 1.0;
}

double g08_r2(double *x) {
	return 1.0 - x[0] + pow(x[1] - 4.0, 2);
}
