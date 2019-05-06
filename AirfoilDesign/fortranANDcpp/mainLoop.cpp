#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <omp.h> 
#include <fstream>		// std::cout, std::endl
#include <iomanip>		// std::setw
#include <string.h>
#include <sstream>
#include <iostream>
using namespace std;
#define cellLen D*sizeof(int)


//put your problem definition and setting header file here
//#include "input/ZDT1-5D.h"

#include "SolverWing.h"

#include "EASCMheaderLoop.h"

//export OMP_STACKSIZE=512M;

float coeff[pps*6];
float XX[pps];

int main(void) {
int loopCnt;

	char loopCntStr[2];
	
	double T1, Ttmp, Tparallel;										/*The timers*/
	
	long long int i,j,k, ind; 								/*Loop Counters*/
	string filename;
	
	


	omp_set_num_threads(20);
	

	float h[D];												/*Cell width*/
	
	int * gloCell = new int[maxArrLen*D]; 					/*All the cells*/
	
	long long int gloCellInd = -1; 							/*gloCell Index*/
	float * gloCellObj = new float[maxArrLen*O];			/*All the cellObjs*/
	
	bool * gloCellFlag = new bool[maxArrLen];				/*Flag for studied cells*/
	
	
	int * solCell = new int[maxArrLen*D]; 					/*All the solution cells*/
	
	long long int solCellInd = 0;
	float * solCellObj = new float[maxArrLen*O];			/*All the solCellObjs*/
	
	int * tmpSolCell = new int[maxTempLen*D]; 				/*Temp sol cells of iteration*/
	float * tmpSolCellObj = new float[maxTempLen*O];
		
	/* Variables to be used in the parallel loop */
	int * todoCell = new int[maxTempLen*D];					/*Cells to be studied*/
	
	long long int todoLen = 0;								/* # cells to be studied*/ 
	bool * solFlag = new bool[maxTempLen];					/* to be added to solCells*/
	int * neiCell = new int[maxTempLen*2*D*D];				/*todoCell Neighbor cells*/
	float * neiCellObj = new float[maxTempLen*2*D*O];		/*Neighbors Objs*/
	bool * destCellFlag = new bool[maxTempLen*2*D];
	
	int subDivInd = 0;
	bool subDivFlag = true;
	time_t t = time(0);
	struct tm * now = localtime( & t );
	int Nrecord[D];
	for (i = 0; i < D; i++) Nrecord[i] = N[i];
	
	// correcting the bounds to include the upper and lower margin.
	//bound_correction(lb, ub, N, subDivSize);
	
	
	int fin_result_dim = 0;
	if (!interm_result) { //in case we are not interested in intermediate results 
		for (i = 0; i<D; i++) {
			if (subDivSize[i]!=0) fin_result_dim=i+1;
		}
	}
	parsec2pointsPrep( XX, coeff, pps);
for (loopCnt = 0; loopCnt <= 0; loopCnt++) {
	for (i = 0; i < D; i++) N[i] = Nrecord[i];
	sprintf(loopCntStr,"%d",loopCnt);
	T1 = omp_get_wtime();
	Ttmp = omp_get_wtime();
	//filename = "mkdir " + fcn;
	//system(filename.c_str());
	ofstream report, log;
	if (log_gen) {
		filename = "./" + fcn + "/"  + fcn + "_log.dat";
		log.open(filename.c_str());
		log << "problem initialized" << endl;
	}
	cout << filename << endl;
	filename = "./" + fcn  + "/" + fcn + "_report.dat";
	report.open(filename.c_str());
	report << "Function name:\n" << fcn << endl;
	report << "Simulation Date:" << endl;
	
		/*initial info in the report file*/
		t = time(0); //get time now
		
		report << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-'
    	     <<  now->tm_mday << endl;
		report << "D\n" << D << "\n" << "O\n" << O << "\nN" << endl;
		for (i = 0; i < D; i++) report << N[i] << ", ";
		report << "\nSubDiv" << endl;
		for (i = 0; i < D; i++) report << subDivSize[i] << ", ";
		report << "\nlb" << endl;
		for (i = 0; i < D; i++) report << lb[i] << ", ";
		report << "\nub" << endl;
		for (i = 0; i < D; i++) report << ub[i] << ", ";
		report << "\n Max array length (maxArrLen):\n" << maxArrLen <<endl;
		report << "Max temporary array length (maxTempLen):\n" << maxTempLen << endl;
	memset(gloCell, -1, maxArrLen*D*sizeof(int));
	gloCellInd = -1;
	memset(gloCellObj,0,maxArrLen*O*sizeof(float));
	memset(gloCellFlag, 0, maxArrLen*sizeof(bool));
	memset(solCell, -1, maxArrLen*D*sizeof(int));
	solCellInd = 0;
	memset(solCellObj,0,maxArrLen*O*sizeof(float));
	memset(todoCell, -1, maxTempLen*D*sizeof(int));
	todoLen = 0;
	subDivInd = 0;
	subDivFlag = true;
	
	
	
	
	
	for (i = 0; i < D; i++)		h[i] = (ub[i]-lb[i])/((float)N[i]);
	Ttmp = omp_get_wtime() - Ttmp;
	if (log_gen) log << "Variables allocated. (" << Ttmp << " s)" << endl;
	
	
	Ttmp = omp_get_wtime();
	
	if (readFromFile) {
		initFromFile(todoCell, todoLen, 
			gloCell, gloCellObj, gloCellFlag, 
			gloCellInd, lb, h, N, fcn, loopCntStr,log_gen, log);
	}
	else {
		if (log_gen) log << "Assigning " << initCellLen << " initial cells..." << endl;
		cout << "Assigning " << initCellLen << " initial cells..." << endl;
		initialize(todoCell, todoLen, gloCell, gloCellObj, gloCellFlag, 
			gloCellInd, initCellLen, lb,h);
	}
	//printf("printing out\n");
	//for (i = 0; i < todoLen; i++) {
	//	printf("todoCell[");
	//	for (j = 0; j < D; j++) {
	//		printf("%4d,",todoCell[i*D+j]);
	//	}
	//	printf("]\n");
	//}
	Ttmp = omp_get_wtime() - Ttmp;
	if (log_gen) log << "Done. (" << Ttmp << " s)" << endl;

	
	//for (i = 0; i < todoLen; i++) 
	//	printf("gloCell[%4d,%4d,%4d],f=[%4.1f,%4.1f,%4.1f]\n",gloCell[i*D],gloCell[i*D+1],gloCell[i*D+2],gloCellObj[i*D],gloCellObj[i*D+1],gloCellObj[i*D+2]);

	
	while (subDivFlag == true) {
		Tparallel = omp_get_wtime();	
		if (log_gen) log << "\nSubDivInd = " << subDivInd << ". Cells list length = "
						<< todoLen << "\n" << endl; 
		cout << "\nSubDivInd = " << subDivInd << ". Cells list length = "
			<< todoLen << "\n" << endl; 
		do {
			Ttmp = omp_get_wtime();	
			
			if (log_gen) log << "**Parallel loop started with " << todoLen 
							<< " cells and their " << todoLen*2*D <<" neighbors..." << endl;
			if (log_gen) cout << "Parallel loop started with " << todoLen 
							<< " cells and their " << todoLen*2*D <<" neighbors..." << endl;
			memset(neiCell, -1, maxTempLen*2*D*D*sizeof(int));
			memset(neiCellObj, 0, maxTempLen*2*D*O*sizeof(float));
			memset(destCellFlag, 0, maxTempLen*2*D*sizeof(bool));
			memset(solFlag, 0, maxTempLen);
			#pragma omp parallel default(none) private(i, j, k, ind) shared(N, lb, h, todoLen, todoCell, gloCell, gloCellObj, gloCellInd, neiCell, neiCellObj, solFlag, destCellFlag, gloCellFlag, XX, coeff) //schedule(guided,2) 
			{
				int * cell = new int[D];						/*Cell Coordinates*/
				float * cellObj = new float[O];					/*Cell Objective*/
				#pragma omp for schedule(dynamic)
				for (i = 0; i < todoLen; i++) {
					//printf("Thread: %d\n",omp_get_thread_num());
		
					memcpy(cell, todoCell+i*D, D*sizeof(int));
					
					funEval(cellObj, cell, (long long int) 1, gloCell, gloCellObj, gloCellInd,lb, h);
					neighFind(neiCell+i*2*D*D, cell, N);
					funEval(neiCellObj+i*2*D*O, neiCell+i*2*D*D, (long long int) (2*D), gloCell, gloCellObj, gloCellInd,lb, h);
					//printf("cell=[%2d,%2d,%2d],funEval=[%4.1f,%4.1f,%4.1f],Thread=%d\nNei=[%2d,%2d,%2d],[%2d,%2d,%2d],...,[%2d,%2d,%2d]\n",cell[0],cell[1],cell[2],cellObj[0],cellObj[1],cellObj[2],omp_get_thread_num(),neiCell[i*2*D*D],neiCell[i*2*D*D+1],neiCell[i*2*D*D+2],neiCell[i*2*D*D+3],neiCell[i*2*D*D+4],neiCell[i*2*D*D+5],neiCell[(i+1)*2*D*D-3],neiCell[(i+1)*2*D*D-2],neiCell[(i+1)*2*D*D-1]);	
					explore(destCellFlag+i*2*D, solFlag+i, cell, cellObj, neiCell+i*2*D*D, neiCellObj+i*2*D*O, N, h, lb);
					ind = binSrch(gloCell, cell, gloCellInd);
					gloCellFlag[ind] = true;
				}
				#pragma omp barrier
				delete[] cell;
				delete[] cellObj;
			}
			Ttmp = omp_get_wtime() - Ttmp;
			if (log_gen) log << "Done. (" << Ttmp << " s)" << endl;	
			/*
			////test print
			int * cell = new int[D];						//Cell Coordinates
			float * cellObj = new float[O];					//Cell Objective
			for (i = 0; i < todoLen; i++) {
				memcpy(cell, todoCell+i*D, D*sizeof(int));	
				funEval(cellObj, cell, (long long int) 1, gloCell, gloCellObj, gloCellInd,lb, h);
				printf("cell=[%2d,%2d,%2d],funEval=[%5.2f,%5.2f,%5.2f],Thread=%d\n",cell[0],cell[1],cell[2],cellObj[0],cellObj[1],cellObj[2],omp_get_thread_num());
				printf("solFlag=%dNeighbors:\n",solFlag[i]);
			delete[] cell;
			delete[] cellObj;
			
				for (j=0; j<2*D; j++) {
					printf("[%2d,%2d,%2d],[%5.2f,%5.2f,%5.2f],Flag=%d\n",neiCell[i*2*D*D+j*D],neiCell[i*2*D*D+j*D+1],neiCell[i*2*D*D+j*D+2],neiCellObj[i*2*D*O+j*O],neiCellObj[i*2*D*O+j*O+1],neiCellObj[i*2*D*O+j*O+2],destCellFlag[i*2*D+j]);
				}
			}
			*/
			////test print
			
			/*inserting funEvals into gloCell*/
			if (log_gen) log << "Updating "<< gloCellInd <<" global cells and their Objs..." << endl; 
			Ttmp = omp_get_wtime();
			insert(gloCell, gloCellInd, gloCellObj, gloCellFlag, neiCell, 
				(long long int) todoLen*2*D*D, neiCellObj);
			Ttmp = omp_get_wtime() - Ttmp;
			if (log_gen) log << "Done.(gloCell pop=" << gloCellInd <<")("<< Ttmp <<" s)" << endl;
			//for (i = 0; i <= gloCellInd; i++) {
			//	printf("gloCell[%4d,%4d,%4d],f=[%4.1f,%4.1f,%4.1f]\n",gloCell[i*D],gloCell[i*D+1],gloCell[i*D+2],gloCellObj[i*D],gloCellObj[i*D+1],gloCellObj[i*D+2]); }
			
			
			/* updating the solution*/
			if (log_gen) log << "Updating the solution and dominancy check (sol pop=" 
							<< solCellInd <<")..." << endl;
			Ttmp = omp_get_wtime();
			updateSol(solCell, solCellObj, solCellInd, tmpSolCell, tmpSolCellObj, 
				todoCell, todoLen, solFlag, gloCell, gloCellInd, gloCellObj);
			Ttmp = omp_get_wtime() - Ttmp;
			if (log_gen) log << "Done.(sol pop=" << solCellInd <<")("<< Ttmp <<" s)" << endl;
			
			//for (i = 0; i < solCellInd; i++) {
			//	printf("solCell[%4d,%4d,%4d],f=[%4.1f,%4.1f,%4.1f]\n",solCell[i*D],solCell[i*D+1],solCell[i*D+2],solCellObj[i*D],solCellObj[i*D+1],solCellObj[i*D+2]); }
			
			
			/*creat a new list for the next iteration*/
			updateToDoCell(todoCell, todoLen, destCellFlag, neiCell, gloCell, 
				gloCellInd, gloCellFlag);
			
			
			/* just for test*///qs(todoCell, (long long int) 0,(long long int) (todoLen-1));
			//printf("todoLen=%lld\n",todoLen);
			//for (i = 0; i < todoLen; i++) 
			//	printf("todoCell[%4d,%4d,%4d]\n",todoCell[i*D],todoCell[i*D+1],todoCell[i*D+2]); 
		} while (todoLen>0);
		Tparallel = omp_get_wtime() - Tparallel;
		if (log_gen) log << "Done. (" << Tparallel <<" s)" << endl;
		cout << "Done. (" << Tparallel <<" s)" << endl;
		
		/* Exporting the solution */
		if (interm_result) {
			if (log_gen) log << "Exporting "<< solCellInd <<" solutions to file..." << endl;
			cout << "Exporting "<< solCellInd <<" solutions to file..." << endl;
			Ttmp = omp_get_wtime();
			exportSol(fcn, loopCntStr, solCell, solCellInd, solCellObj, lb, h, subDivInd);
			Ttmp = omp_get_wtime() - Ttmp;
			if (log_gen) log << "Done.("<< Ttmp <<" s)" << endl;
			cout << "Done.("<< Ttmp <<" s)" << endl;
		} 
		else if (fin_result_dim == subDivInd){
			if (log_gen) log << "Exporting "<< solCellInd <<" solutions to file..." << endl;
			cout << "Exporting "<< solCellInd <<" solutions to file..." << endl;
			Ttmp = omp_get_wtime();
			exportSol(fcn, loopCntStr, solCell, solCellInd, solCellObj, lb, h, subDivInd);
			Ttmp = omp_get_wtime() - Ttmp;
			if (log_gen) log << "Done.("<< Ttmp <<" s)" << endl;
			cout << "Done.("<< Ttmp <<" s)" << endl;
		}
		
		
		/* SubDivision */
		while (subDivInd < D && subDivSize[subDivInd] == 0) {
			subDivInd++;
		}
		if (subDivInd == D) 
			//break;
			subDivFlag = false;
		else {
			if (log_gen) log << "Doing subDiv on Dim=" << subDivInd << " with " 
							<< gloCellInd << " total cells..." << endl;
			Ttmp = omp_get_wtime();
			subDiv(todoCell, todoLen, gloCell, gloCellInd, 
				solCell, solCellInd, solCellObj,
				N, h, ub, lb, subDivInd, subDivSize);
			Ttmp = omp_get_wtime() - Ttmp;
			if (log_gen) {
				log << "Done.(" << Ttmp <<" s)\nN=[" << endl;
				for (i=0; i<D; i++) log<<N[i]<<", "; log<<"]"<<endl;
			}
			subDivInd++;
			//printf("N=[%d,%d,%d],h=[%f,%f,%f]\n",N[0],N[1],N[2],h[0],h[1],h[2]);
		}
	}
	report << "Total number of function evaluations:\n" << gloCellInd <<
		"\nTotal number of solution points:\n" << solCellInd << endl;
	T1 = omp_get_wtime() - T1;
	report << "Total time (s):\n" << T1 <<" s"<< endl;
	log.close();
	report.close();
	
}	


	delete[] gloCell;
	delete[] gloCellObj;
	delete[] gloCellFlag;
	delete[] solCell;
	delete[] solCellObj;
	delete[] tmpSolCell;
	delete[] tmpSolCellObj;
	delete[] todoCell;
	delete[] solFlag;
	delete[] neiCell;
	delete[] neiCellObj;
	delete[] destCellFlag;
	

}
