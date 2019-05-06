/***************************************************************************
                          variables.h  -  description
                             -------------------
    begin                : Mon Jan 29 2001
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
#include "SolverWing.h"
unsigned int num_fun(unsigned int fun){
	
	unsigned int num = 1;
	
	switch (fun) {
 		case 100: num = 2;	break;
 		case 200: num = 2;	break;
 		case 300: num = 2;	break;
 		case 400: num = 2;	break;
 		case 500: num = 2; break;
 		case 600: num = 3; break;
 		case 700: num = 2; break;
 		case 800: num = 3; break;
 		case 900: num = 2; break;
 		case 1000: num = 2; break;
	        case 10100:num =2;break;
	        case 10200:num =2;break;
	        case 10300:num =2;break;
	        case 10400:num =2;break;
	        case 10500:num =2;break;
	        case 20100:num = 2;break; // ZDT3-10D
	        case 20200:num = 2;break; // ZDT3-10D
	        case 20300:num = 3;break; //airfoil
	}
	
	return num;
	
}

unsigned int num_dim(unsigned int fun){
	
	unsigned int num = 0;
	
	switch (fun) {
 		case 1: num = 13;	break;
 		case 2: num = 20; break;
 		case 3: num = 5;	break;
 		case 4: num = 2;	break;
 		case 5: num = 10;	break;
 		case 6: num = 2;	break;
 		case 7: num = 4;	break;
 		case 8: num = 4;	break;
 		case 9: num = 2;	break;
 		case 10: num = 1;	break;
 		case 11: num = 2;	break;
 		case 100: num = 1; break;
 		case 200: num = 3; break;
 		case 300: num = 3; break;
 		case 400: num = 1; break;
 		case 500: num = 2; break;
 		case 600: num = 2; break;
 		case 700: num = 2; break;
 		case 800: num = 2; break;
 		case 900: num = 2; break;
 		case 1000: num = 2; break;
	        case 10100:num=2;break;
	        case 10200:num=2;break;
	        case 10300:num=4;break;
	        case 10400:num=2;break;
	        case 10500:num=2;break;
	        case 20100:num=10;break; // ZDT3-10D
	        case 20200:num=10;break; // ZDT3-10D
	        case 20300:num=12;break; // airfoil

	}
	
	return num;
	
}


unsigned int num_rest(unsigned int fun){
	
	unsigned int num = 0;
	
	switch (fun) {
 		case 1: num = 9;	break;
 		case 2: num = 2;	break;
 		case 3: num = 6;	break;
 		case 4: num = 2;	break;
 		case 5: num = 8;	break;
 		case 6: num = 2;	break;
 		case 7: num = 4;	break;
 		case 8: num = 7;	break;
	}
	
	return num;
	
}

void initpop(double *pop, unsigned int fun, unsigned int M, unsigned int D) {
  double F = 10.0;
  double E = 2*pow(10,5);
  double L = 200.0;
  double r = 10.0;
	
	unsigned int i, j;

	switch (fun) {
  	case 1:	for(i = 0; i < M; i++)
	  					for(j = 0; j < 9; j++)
  							pop[(D * i) + j] = RandomDouble(0.0, 1.0);
  					for(i = 0; i < M; i++)
  						for(j = 9; j < 12; j++)
  							pop[(D * i) + j] = RandomDouble(0.0, 100.0);
 						for(i = 0; i < M; i++)
	 						pop[(D * i) + 12] = RandomDouble(0.0, 1.0);
  					break;
  	case 2:	for(i = 0; i < M; i++)
  						for(j = 0; j < 20; j++)
  							pop[(D * i) + j] = RandomDouble(0.0, 10.0);
  					break;
  	case 3:	for(i = 0; i < M; i++)
	  					pop[(D * i) + 0] = RandomDouble(78.0, 102.0);
  					for(i = 0; i < M; i++)
  						pop[(D * i) + 1] = RandomDouble(33.0, 45.0);
	  				for(i = 0; i < M; i++)
	  					for(j = 2; j < 5; j++)	
		  					pop[(D * i) + j] = RandomDouble(27.0, 45.0);
  					break;
  	case 4: for(i = 0; i < M; i++)
	  					pop[(D * i) + 0] = RandomDouble(13.0, 100.0);
  					for(i = 0; i < M; i++)
  						pop[(D * i) + 1] = RandomDouble(0.0, 100.0);
  					break;
  	case 5:	for(i = 0; i < M; i++)
  						for(j = 0; j < 10; j++)
  							pop[(D * i) + j] = RandomDouble(-10.0, 10.0);
  					break;
  	case 6:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
  							pop[(D * i) + j] = RandomDouble(0.0, 10.0);
  					break;
  	case 7:	for(i = 0; i < M; i++) {
  						for(j = 0; j < 2; j++)
  							pop[(D * i) + j] = RandomDouble(1.0, 99.0);
  						for(j = 2; j < 4; j++)
  							pop[(D * i) + j] = RandomDouble(10.0, 200.0);  					
  					}
  					break;
  	case 8:	for(i = 0; i < M; i++) {
 							pop[(D * i) + 0] = RandomDouble(0.1, 2.0);
  						for(j = 1; j < 3; j++)
  							pop[(D * i) + j] = RandomDouble(0.1, 10.0);
 							pop[(D * i) + 3] = RandomDouble(0.1, 2.0);
  					}
  					break;
  	case 9:	for(i = 0; i < M; i++)
 							pop[(D * i) + 0] = RandomDouble(-3.0, 3.0);
 						for(i = 0; i < M; i++)
 							pop[(D * i) + 1] = RandomDouble(-2.0, 2.0);
  					break;
  	case 10: 	for(i = 0; i < M; i++)
 								pop[(D * i) + 0] = RandomDouble(-5.0, 5.0);
  						break;
  	case 11: 	for(i = 0; i < M; i++)
  							for(j = 0; j < 2; j++)
  								pop[(D * i) + j] = RandomDouble(-65.536, 65.536);
  						break;
	case 100:	for(i = 0; i < M; i++)
 						pop[(D * i) + 0] = RandomDouble(-5.0, 10.0);
  					break;
  	case 200:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(-5.0, 5.0);
  					break;
  	case 300:	for(i = 0; i < M; i++)
  						for(j = 0; j < 3; j++)
 							pop[(D * i) + j] = RandomDouble(-4.0, 4.0);
  					break;	
  	case 400:	for(i = 0; i < M; i++)
 						pop[(D * i) + 0] = RandomDouble(pow(-10.0, 5), pow(10.0, 5));
  					break;	
  	case 500:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(-3.1416, 3.1416);
  					break;
  	case 600:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(-30.0, 30.0);
  					break;
  	case 700:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(0.0, 1.0);
  					break;
	case 800:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(-400.0, 400.0);
  					break;
	case 900:	for(i = 0; i < M; i++)
 							pop[(D * i) + 0] = RandomDouble(0.0, 1.0);
 						for(i = 0; i < M; i++)
 							pop[(D * i) + 1] = RandomDouble(-30.0, 30.0);
  					break;
	case 1000:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(0.1, 1.0);
  					break;
	case 10100:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(0.0, 7.0);
  					break;
	case 10200:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
 							pop[(D * i) + j] = RandomDouble(-20.0, 20.0);
  					break;
	case 10300:
	  for(i = 0; i < M; i++){
	    pop[(D * i) + 0]=RandomDouble((double)(F/r),(double)3.0*(F/r));
	    pop[(D * i) + 1]=RandomDouble((double)sqrt(2.0)*(F/r),(double)3.0*(F/r));
	    pop[(D * i) + 2]=RandomDouble((double)sqrt(2.0)*(F/r),(double)3.0*(F/r));
	    pop[(D * i) + 3]=RandomDouble((double)(F/r),(double)3.0*(F/r));
	    
	  }
	  break;
	case 10400:	for(i = 0; i < M; i++)
  						for(j = 0; j < 2; j++)
						  pop[(D * i) + j] = RandomDouble(0.1, 1.0);
	break;
	case 10500:	
	  for(i = 0; i < M; i++){
	    pop[(D * i) + 0] = RandomDouble(0.0, 1.0);
	    pop[(D * i) + 1] = RandomDouble(-30.0, 30.0);

	  }
	break;
	case 20100:	for(i = 0; i < M; i++) // ZDT3-10D
  		for(j = 0; j < 10; j++)
 			pop[(D * i) + j] = RandomDouble(0.0, 1.0);
  	break;
  	case 20200:	for(i = 0; i < M; i++) // ZDT3-10D
  		for(j = 0; j < 10; j++)
 			pop[(D * i) + j] = RandomDouble(0.0, 1.0);
  	break;
	
	case 20300:	for(i = 0; i < M; i++) // airfoil
  		for(j = 0; j < D; j++)
 			pop[(D * i) + j] = RandomDouble(lb[j], ub[j]);
  	break;
  	}
}

unsigned int  constraints(double *x, unsigned int fun) {
	
  unsigned int i, j,violadas=0;
  switch(fun){  
  /*
  case 10100: 
    violadas=0;
    //    cout<<"                                                       cotz"<<endl;
    if( 0 < (double)(x[0]/6.0 + x[1] - 13.0/2.0 )) violadas ++;
    if( 0 < (double)(x[0]/2.0 + x[1] - 15.0/2.0) ) violadas ++;
    if( 0 < (double)(5.0*x[0] + x[1] - 30.0) ) violadas ++;
    //cout<<violadas<<endl;
    return violadas;
  case 10200: 
    //    return 0;
    //cout<<"       
    violadas++;
    if( 0 < pow(x[0],2) + pow(x[1],2) - 225 ) violadas++;
    if( 0 < x[0] - 3*x[1] + 10 ) violadas++;
    //cout<<"hola"<<endl;
    return violadas;
    break;
    */
  default:
    return 0;
  }
}

// constraints2 resturns the total number of constraint violations
unsigned int  constraints2(double *x, double *fitness, unsigned int fun) {
	unsigned int violations = 0;
	switch(fun){
		case 20200: 
			//if (fitness[1] > 0.6 && fitness[0] < 0.7) violations++;
			if (x[1] > 0.02 && x[1] < 0.03) violations++;
			return violations;
		case 20300:
			if (fitness[0] == inf) violations++;
			if (fitness[1] == inf) violations++;
			if (fitness[2] == inf) violations++;
			return violations;
		dafault:
			return 0;
	}
}



void chfun(double *fitness, double *pop, unsigned int fun) {
  //CONT_FUN++;
	switch (fun) {
/*
  	case 1:	return g01_IEEE(pop);
  	case 2:	return g02_IEEE(pop);
  	case 3:	return g04_IEEE(pop);
  	case 4:	return g06_IEEE(pop);
  	case 5:	return g07_IEEE(pop);
  	case 6:	return g08_IEEE(pop);
  	case 7:	return g01_E(pop);
  	case 8:	return g02_E(pop);
  	case 9:	return g01_SR(pop);
  	case 10:	return g02_SR(pop);
  	case 11:	return g03_SR(pop);
  	case 100:	 return f1_1_mo(pop);
  	case 101: return f2_1_mo(pop);
  	case 200:	return f1_2_mo(pop);
  	case 201: return f2_2_mo(pop);
  	case 300:	return f1_3_mo(pop);
  	case 301: return f2_3_mo(pop);
  	case 400:	return f1_4_mo(pop);
  	case 401: return f2_4_mo(pop);
  	case 500: return f1_5_mo(pop);
  	case 501: return f2_5_mo(pop);
  	case 600: return f1_6_mo(pop);
  	case 601: return f2_6_mo(pop);
  	case 602: return f3_6_mo(pop);
  	case 700: return f1_7_mo(pop);
  	case 701: return f2_7_mo(pop);
  	case 800: return f1_8_mo(pop);
  	case 801: return f2_8_mo(pop);
  	case 802: return f3_8_mo(pop);
  	case 900: return f1_9_mo(pop);
  	case 901: return f2_9_mo(pop);
  	case 1000: return f1_10_mo(pop);
  	case 1001: return f2_10_mo(pop);
  	case 10100: return f1_10100_mo(pop);
  	case 10101: return f2_10100_mo(pop);
  	case 10200: return f1_10200_mo(pop);
  	case 10201: return f2_10200_mo(pop);
  	case 10300: return f1_10300_mo(pop);
  	case 10301: return f2_10300_mo(pop);
  	case 10400: return f1_10400_mo(pop);
  	case 10401: return f2_10400_mo(pop);
  	case 10500: return f1_10500_mo(pop);
  	case 10501: return f2_10500_mo(pop);
  	case 20100: return f1_ZDT3_10D(pop); // ZDT3-10D
  	case 20101: return f2_ZDT3_10D(pop); // ZDT3-10D
*/
	case 20200: 
				ZDT3_10D(fitness, pop); // ZDT3-10D
				break;
	case 20300: 
				return f(fitness, pop); // airfoil
				break;
 	}

//return;
}

double chconst(double *pop, unsigned int fun, unsigned int cons) {
	
	switch (fun) {
	case 1:	switch(cons) {
  					case 1: return g01_r1(pop);
  					case 2: return g01_r2(pop);
  					case 3: return g01_r3(pop);
  					case 4: return g01_r4(pop);
  					case 5: return g01_r5(pop);
  					case 6: return g01_r6(pop);
  					case 7: return g01_r7(pop);
  					case 8: return g01_r8(pop);
  					case 9: return g01_r9(pop);
					}
	case 2:	switch(cons) {
  					case 1: return g02_r1(pop);
  					case 2: return g02_r2(pop);
					}
	case 3:	switch(cons) {
  					case 1: return g04_r1(pop);
  					case 2: return g04_r2(pop);
  					case 3: return g04_r3(pop);
  					case 4: return g04_r4(pop);
  					case 5: return g04_r5(pop);
  					case 6: return g04_r6(pop);
					}
	case 4:	switch(cons) {
  					case 1: return g06_r1(pop);
  					case 2: return g06_r2(pop);
					}
	case 5:	switch(cons) {
  					case 1: return g07_r1(pop);
  					case 2: return g07_r2(pop);
  					case 3: return g07_r3(pop);
  					case 4: return g07_r4(pop);
  					case 5: return g07_r5(pop);
  					case 6: return g07_r6(pop);
  					case 7: return g07_r7(pop);
  					case 8: return g07_r8(pop);
					}
	case 6:	switch(cons) {
  					case 1: return g08_r1(pop);
  					case 2: return g08_r2(pop);
					}
	case 7:	switch(cons) {
  					case 1: return g01_Er1(pop);
  					case 2: return g01_Er2(pop);
  					case 3: return g01_Er3(pop);
  					case 4: return g01_Er4(pop);
					}
	case 8:	switch(cons) {
  					case 1: return g02_Er1(pop);
  					case 2: return g02_Er2(pop);
  					case 3: return g02_Er3(pop);
  					case 4: return g02_Er4(pop);
  					case 5: return g02_Er5(pop);
  					case 6: return g02_Er6(pop);
  					case 7: return g02_Er7(pop);
					}
	}

	return 0;
}


void keepin(double *pop, double *vel, unsigned int fun, unsigned int M, unsigned int D) {
	unsigned int i, j;
	double linf[D], lsup[D];
	  double F = 10.0;
	  double E = 2*pow(10,5);
	  double L = 200.0;
	  double r = 10.0;
	
	switch (fun) {
    case 1:	for(i = 0; i < 9; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 1.0;
          	}
          	for(i = 9; i < 12; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 100.0;
          	}
     				linf[12] = 0.0;
           	lsup[12] = 1.0;
          	break;
    case 2:	for(i = 0; i < D; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 10.0;
          	}
          	break;
  	case 3: linf[0] = 78.0;
          	linf[1] = 33.0;
          	linf[2] = 27.0;
          	linf[3] = 27.0;
          	linf[4] = 27.0;
          	
          	lsup[0] = 102.0;
          	lsup[1] = 45.0;
          	lsup[2] = 45.0;
          	lsup[3] = 45.0;
          	lsup[4] = 45.0;
          	break;
		case 4: linf[0] = 13.0;
          	linf[1] = 0.0;
          	
          	lsup[0] = 100.0;
          	lsup[1] = 100.0;
          	break;
    case 5:	for(i = 0; i < D; i++) {
      				linf[i] = -10.0;
            	lsup[i] = 10.0;
          	}
          	break;
    case 6:	for(i = 0; i < D; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 10.0;
          	}
          	break;
    case 7: linf[0] = 1.0;
          	linf[1] = 1.0;
          	linf[2] = 10.0;
          	linf[3] = 10.0;
          	
          	lsup[0] = 99.0;
          	lsup[1] = 99.0;
          	lsup[2] = 200.0;
          	lsup[3] = 200.0;
          	break;
		case 8:	for(i = 0; i < D; i++) {
      				linf[i] = 0.1;
          	}
          	lsup[0] = 2.0;
          	lsup[1] = 10.0;
          	lsup[3] = 10.0;
          	lsup[4] = 2.0;
          	break;
    case 9: linf[0] = -3.0;
          	linf[1] = -2.0;
          	
          	lsup[0] = 3.0;
          	lsup[1] = 2.0;
          	break;
		case 10:	linf[0] = -5.0;
          	
          		lsup[0] = 5.0;
          	break;
    case 11:	for(i = 0; i < D; i++) {
      					linf[i] = -65.536;
            		lsup[i] = 65.536;
          		}
          	break;
    case 100:	linf[0] = -5.0;
          	
          		lsup[0] = 10.0;
          	break;
    case 200:	for(i = 0; i < D; i++) {
    						linf[i] = -5.0;
          		lsup[i] = 5.0;
          		}
          	break; 	
    case 300:	for(i = 0; i < D; i++) {
    					linf[i] = -4.0;
          		lsup[i] = 4.0;
          		}
          	break;
    case 400:	linf[0] = pow(-10.0, 5);
          		lsup[0] = pow(10.0, 5);
          	break;
    case 500:	for(i = 0; i < D; i++) {
    					linf[i] = -3.1416;
          		lsup[i] = 3.1416;
          		}
          	break;
     case 600:	for(i = 0; i < D; i++) {
    					linf[i] = -30.0;
          		lsup[i] = 30.0;
          		}
          	break;
	case 700:	for(i = 0; i < D; i++) {
    					linf[i] = 0.0;
          		lsup[i] = 1.0;
          		}
          	break;
	case 800:	for(i = 0; i < D; i++) {
    					linf[i] = -400.0;
          		lsup[i] = 400.0;
          		}
          	break;
     case 900: linf[0] = 0.0;
          	linf[1] = -30.0;
          	
          	lsup[0] = 1.0;
          	lsup[1] = 30.0;
          	break;
	case 1000:	for(i = 0; i < D; i++) {
    					linf[i] = 0.1;
          		lsup[i] = 1.0;
          		}
          	break;

 	
	case 10100:	for(i = 0; i < D; i++) {
    					linf[i] = 0.0;
          		lsup[i] = 7.0;
          		}
          	break;

	case 10200:	for(i = 0; i < D; i++) {
    					linf[i] = -20.0;
          		lsup[i] = 20.0;
          		}
          	break;

	case 10300:
	  
	  linf[0]=(double)(F/r);//0.05
	  lsup[0]=(double)3.0*(F/r);//0.15
	  linf[1]=(double)sqrt(2)*(F/r);//0.070710678
	  lsup[1]=(double)3.0*(F/r);//0.15
	  linf[2]=(double)sqrt(2)*(F/r);//0.070710678
	  lsup[2]=(double)3.0*(F/r);//0.15
	  linf[3]=(double)(F/r);//0.05
	  lsup[3]=(double)3.0*(F/r);//0.015

	  break;

	case 10400:	for(i = 0; i < D; i++) {
    					linf[i] = 0.1;
          		lsup[i] = 1.0;
          		}
          	break;

	case 10500:	
	  linf[0]=0.0;
	  linf[1]=-30.0;
	  lsup[0]=1.0;
	  lsup[1]=30.0;

	  break;
	case 20100:	for(i = 0; i < D; i++) { // ZDT3-10D
    				linf[i] = 0.0;
          			lsup[i] = 1.0;
          		}
	break;
	
	case 20200:	for(i = 0; i < D; i++) { // ZDT3-10D
    				linf[i] = 0.0;
          			lsup[i] = 1.0;
          		}
	break;
	
	case 20300:	for(i = 0; i < D; i++) { // airfoil
    				linf[i] = lb[i];
          			lsup[i] = ub[i];
          		}
	break;


 	}
  	
 	for(i = 0; i < M; i++) {
 		for(j = 0; j < D; j++) {
 			if (pop[(D * i) + j] < linf[j]) {
 				pop[(D * i) + j] = linf[j];
 				vel[(D * i) + j] = -vel[(D * i) + j];
 			}
 			if (pop[(D * i) + j] > lsup[j]) {
 				pop[(D * i) + j] = lsup[j];
 				vel[(D * i) + j] = -vel[(D * i) + j];
 			}
  			
 		}
 	}
 	
}



void ranges(unsigned int fun,int D,double *linf,double * lsup){
	unsigned int i, j;
	  double F = 10.0;
	  double E = 2*pow(10,5);
	  double L = 200.0;
	  double r = 10.0;
	
	switch (fun) {
    case 1:	for(i = 0; i < 9; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 1.0;
          	}
          	for(i = 9; i < 12; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 100.0;
          	}
     				linf[12] = 0.0;
           	lsup[12] = 1.0;
          	break;
    case 2:	for(i = 0; i < D; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 10.0;
          	}
          	break;
  	case 3: linf[0] = 78.0;
          	linf[1] = 33.0;
          	linf[2] = 27.0;
          	linf[3] = 7.0;
          	linf[4] = 27.0;
          	
          	lsup[0] = 102.0;
          	lsup[1] = 45.0;
          	lsup[2] = 45.0;
          	lsup[3] = 45.0;
          	lsup[4] = 45.0;
          	break;
		case 4: linf[0] = 13.0;
          	linf[1] = 0.0;
          	
          	lsup[0] = 100.0;
          	lsup[1] = 100.0;
          	break;
    case 5:	for(i = 0; i < D; i++) {
      				linf[i] = -10.0;
            	lsup[i] = 10.0;
          	}
          	break;
    case 6:	for(i = 0; i < D; i++) {
      				linf[i] = 0.0;
            	lsup[i] = 10.0;
          	}
          	break;
    case 7: linf[0] = 1.0;
          	linf[1] = 1.0;
          	linf[2] = 10.0;
          	linf[3] = 10.0;
          	
          	lsup[0] = 99.0;
          	lsup[1] = 99.0;
          	lsup[2] = 200.0;
          	lsup[3] = 200.0;
          	break;
		case 8:	for(i = 0; i < D; i++) {
      				linf[i] = 0.1;
          	}
          	lsup[0] = 2.0;
          	lsup[1] = 10.0;
          	lsup[3] = 10.0;
          	lsup[4] = 2.0;
          	break;
    case 9: linf[0] = -3.0;
          	linf[1] = -2.0;
          	
          	lsup[0] = 3.0;
          	lsup[1] = 2.0;
          	break;
		case 10:	linf[0] = -5.0;
          	
          		lsup[0] = 5.0;
          	break;
    case 11:	for(i = 0; i < D; i++) {
      					linf[i] = -65.536;
            		lsup[i] = 65.536;
          		}
          	break;
    case 100:	linf[0] = -5.0;
          	
          		lsup[0] = 10.0;
          	break;
    case 200:	for(i = 0; i < D; i++) {
    						linf[i] = -5.0;
          		lsup[i] = 5.0;
          		}
          	break; 	
    case 300:	for(i = 0; i < D; i++) {
    					linf[i] = -4.0;
          		lsup[i] = 4.0;
          		}
          	break;
    case 400:	linf[0] = pow(-10.0, 5);
          		lsup[0] = pow(10.0, 5);
          	break;
	case 500:
	  for(i = 0; i < D; i++) {
	    linf[i] = -3.1416;
	    lsup[i] = 3.1416;
	  }
          	break;
     case 600:	for(i = 0; i < D; i++) {
    					linf[i] = -30.0;
          		lsup[i] = 30.0;
          		}
          	break;
	case 700:	for(i = 0; i < D; i++) {
    					linf[i] = 0.0;
          		lsup[i] = 1.0;
          		}
          	break;
	case 800:	for(i = 0; i < D; i++) {
    					linf[i] = -400.0;
          		lsup[i] = 400.0;
          		}
          	break;
     case 900: linf[0] = 0.0;
          	linf[1] = -30.0;
          	
          	lsup[0] = 1.0;
          	lsup[1] = 30.0;
          	break;
	case 1000:	for(i = 0; i < D; i++) {
    					linf[i] = 0.1;
          		lsup[i] = 1.0;
          		}
          	break;
	case 10100:	for(i = 0; i < D; i++) {
    					linf[i] = 0.0;
          		lsup[i] = 7.0;
          		}
          	break;
	case 10200:	for(i = 0; i < D; i++) {
    					linf[i] = -20.0;
          		lsup[i] = 20.0;
          		}
          	break;
	case 10300:
	  
	  linf[0]=(double)(F/r);//0.05
	  lsup[0]=(double)3.0*(F/r);//0.15
	  linf[1]=(double)sqrt(2)*(F/r);//0.070710678
	  lsup[1]=(double)3.0*(F/r);//0.15
	  linf[2]=(double)sqrt(2)*(F/r);//0.070710678
	  lsup[2]=(double)3.0*(F/r);//0.15
	  linf[3]=(double)(F/r);//0.05
	  lsup[3]=(double)3.0*(F/r);//0.015

	  break;

	case 10400:	for(i = 0; i < D; i++) {
    					linf[i] = 0.1;
          		lsup[i] = 1.0;
          		}
          	break;

	case 10500:	
	  linf[0]=0.0;
	  linf[1]=-30.0;
	  lsup[0]=1.0;
	  lsup[1]=30.0;

	  break;
	case 20100:	for(i = 0; i < D; i++) { // ZDT3-10D
    				linf[i] = 0.0;
          			lsup[i] = 1.0;
          		}
	break;
	case 20200:	for(i = 0; i < D; i++) { // ZDT3-10D
    				linf[i] = 0.0;
          			lsup[i] = 1.0;
          		}
	break;
	
	case 20300:	for(i = 0; i < D; i++) { // airfoil
    				linf[i] = lb[i];
          			lsup[i] = ub[i];
          		}
	break;

 	}
}
