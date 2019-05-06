/***************************************************************************
                          fun.h  -  description
                             -------------------
    begin                : Wed Nov 15 2000
    copyright            : (C) 2000 by Max Salazar
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

/**********************************
 *																*
 *	Design of a Pressure Vessel		*
 *																*
 **********************************/

//double PHI = 3.141592653589793238462;

//Ecuaciones
double g01_Er1(double *x) {
	return -(floor(x[0]) * 0.0625) + (0.0193 * x[2]);
}

double g01_Er2(double *x) {
	return -(floor(x[1]) * 0.0625) + (0.00954 * x[2]);
}

double g01_Er3(double *x) {
	return (-(PHI * pow(x[2], 2) * x[3])) - (((double)(4.0/3.0)) * PHI * pow(x[2], 3)) + 1296000.0 ;
}

double g01_Er4(double *x) {
	return x[3] - 240.0;
}

//Funcion del Problema 1 Coello
double g01_E(double *x) {
	return (0.6224 * (floor(x[0]) * 0.0625) * x[2] * x[3]) + (1.7781 * (floor(x[1]) * 0.0625) * pow(x[2], 2)) +
				 (3.1661 * pow((floor(x[0]) * 0.0625), 2) * x[3]) + (19.84 * pow((floor(x[0]) * 0.0625), 2) * x[2]);
}

//////////////////////////////////////////////////////////


/************************
 *											*
 *	Welded Beam Design	*
 *											*
 ************************/
float P = 6000.0, L = 14.0, taumax = 13600.0, sigmax = 30000.0, deltamax = 0.25,
			E = 30e6, G = 12e6;

//Ecuaciones
double Mf(double *x) {
	return P * (L + (x[1] / 2.0));
}

double Rf(double *x) {
	return sqrt((pow(x[1], 2) / 4.0) + pow(((x[0] + x[2]) / 2.0), 2));
}

double Jf(double *x) {
	return 2.0 * ( (sqrt(2) * x[0] * x[1]) * ((pow(x[1], 2) / 12.0) + pow(((x[0] + x[2]) / 2.0), 2) ) );
}

double taup(double *x) {
	return P / (sqrt(2) * x[0] * x[1]);
}

double taupp(double *x) {
	return (Mf(x) * Rf(x)) / Jf(x);
}

double tau(double *x) {
	return sqrt( pow(taup(x), 2 ) +
							(2.0 * taup(x) * taupp(x) * (x[1] / (2.0 * Rf(x))) ) +
							 pow(taupp(x), 2) );
}

double sig(double *x) {
	return ((6.0 * P * L) / (x[3] * pow(x[2], 2)));
}

double delta(double *x) {
	return ((4.0 * P * pow(L, 3)) / (E * pow(x[2], 3) * x[3]));
}

double Pcf(double *x) {
	return ((4.013 * E * sqrt((pow(x[2], 2) * pow(x[3], 6)) / 36.0)) / pow(L, 2)) *
					(1.0 - ((x[2] / (2.0 * L)) * sqrt(E / (4.0 * G))) );
}

//Restricciones
double g02_Er1(double *x) {
	return tau(x) - taumax;
}

double g02_Er2(double *x) {
	return sig(x) - sigmax;
}

double g02_Er3(double *x) {
	return x[0] - x[3];
}

double g02_Er4(double *x) {
	return (0.10471 * pow(x[0], 2)) + ((0.04811 * x[2] * x[3]) * (14.0 + x[1])) - 5.0;
}

double g02_Er5(double *x) {
	return 0.125 - x[0];
}

double g02_Er6(double *x) {
	return delta(x) - deltamax;
}

double g02_Er7(double *x) {
	return P - Pcf(x);
}

//Funcion del Problema 2 Coello
double g02_E(double *x) {
	return (1.10471 * pow(x[0], 2) * x[1]) + (0.04811 * x[2] * x[3] * (14.0 + x[1]));
}
