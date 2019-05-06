/***************************************************************************
                          fun-SR.h  -  description
                             -------------------
    begin                : Fri May 11 2001
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

double g01_SR(double *x) {
        return ((4.0 - (2.1 * pow(x[0], 2)) + (pow(x[0], 4) / 3.0)) * pow(x[0], 2)) + (x[0] * x[1]) + ((-4.0 + (4.0 * pow(x[1], 2))) * pow(x[1], 2));
}

double g02_SR(double *x) {
        return /*(*/0.65 - (0.75 / (1.0 + pow(x[0], 2))) - (0.65 * x[0] * atan(1.0 / x[0]));//) * -1.0;
}

double g03_SR(double *x) {
        int i, j, a[2][25];
        double aux, sum_f1, sum_f = 0.0, K = 500.0;

        //Crea la matriz a
        for(i = 0; i < 5; i++)
                for(j = 0; j < 5; j++) {
                        a[0][j + (5 * i)] = -32 + (16 * j);
                        a[1][j + (5 * i)] = -32 + (16 * i);
                }

        for(j = 0; j < 25; j++) {
                sum_f1 = 0.0;
                for(i = 0; i < 2; i++) {
                        aux = x[i] - a[i][j];
                        sum_f1 = sum_f1 + pow(aux, 6);
                }
                sum_f1 = sum_f1 + (j + 1.0);
                sum_f = sum_f + (1.0 / sum_f1);
        }

        return 1.0 / ((1.0 / K) + sum_f);
}
