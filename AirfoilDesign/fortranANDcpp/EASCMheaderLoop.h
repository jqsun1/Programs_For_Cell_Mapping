#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

// This file includes all the necessary subroutines in the EA+SCM algorithm

////////////////////////////////////////////////////////////////////////////////
 /* Random number generator */
float rnd() {
	return ((float)rand()/(float)RAND_MAX);
}
float rnd(float min, float max) {
	float out;
	out = ((float)rand()/(float)RAND_MAX*((float)max-(float)min))+(float)min;
	return out;
}
////////////////////////////////////////////////////////////////////////////////
void bound_correction(float * lb, float * ub, int * N, int * subDiv) {
	float x;
	for (int i = 0; i<D; i++) {
		x = ((double)ub[i] - (double)lb[i]) / (double)(N[i]*(2*subDiv[i]+1)*2-2);
		ub[i] += x;
		lb[i] -= x;
	}
}
////////////////////////////////////////////////////////////////////////////////
/* Z to X mapping */
void z2x(float * x, int * z, float * lb, float * h){
//printf("The lb value is:%f,%f,%f\n",lb[0],lb[1],lb[2]);
	for(int i = 0; i < D; i++)
		x[i] = lb[i]+h[i]*z[i]+0.5f*h[i];
}
////////////////////////////////////////////////////////////////////////////////
/* X to Z mapping */
void x2z(int * z, float * x, float * lb, float * h, int * N) {
	for (int i = 0; i<D; i++) {
		z[i] = (int)floor((x[i]-lb[i])/h[i]);
	}
	for (int i = 0; i<D; i++) {
		if (z[i] < 0 || z[i] >= N[i]) memset(z, -1, D*sizeof(int));
	}
}
////////////////////////////////////////////////////////////////////////////////
// comparison of two coordinates a and b in the array
// if a>b returns 1
// if a<b returns -1
// if a=b returns 0
int cmp(int * p1, int * p2) {
	int i=0; //counter
	int out = 0; // output
	
	while (i < D) {
		if (p1[i] > p2[i]) {
			out = 1;
			break;
		}
		else if (p1[i] < p2[i]) {
			out = -1;
			break;
		}
		else i++;
	}
	//printf("comparing [%5d,%5d],[%5d,%5d] returned: %5d\n",p1[0],p1[1],p2[0],p2[1],out);
	return out;
}
////////////////////////////////////////////////////////////////////////////////
// The parallel quick sort algorithm for the arrays.
void qs(int *v, long long int first,long long int last) {
  long long int start[2], end[2], pivotInd;
  int * pivot = new int[D];
  int i, tmp[D];
	if (first < last) {
		start[1] = first;
		end[0] = last;
		pivotInd = (first + last) / 2;
		memcpy(pivot, v+(pivotInd*D), cellLen);
		while (start[1] <= end[0]) {

			while (cmp(pivot, v+start[1]*D)>0) //v[start[1]] < pivot
				start[1]++;
			while (cmp(v + end[0]*D, pivot)>0)//pivot < v[end[0]]
				end[0]--;
			if (start[1] <= end[0]) {
				memcpy(tmp,v+ D*start[1], cellLen);
				memcpy(v+ D*start[1],v+ D*end[0] , cellLen);
				memcpy(v+ D*end[0], tmp , cellLen);
				start[1]++;
				end[0]--;
			}
		}
		start[0] = first; 
		end[1]   = last; 
		#pragma omp parallel 
		{
			#pragma omp for nowait
			for(i = 0; i <= 1; i++) {
				qs(v, start[i], end[i]);
			}
		}
	}
	delete[] pivot;
}
void qs(int *v, long long int first,long long int last, float * funEval) {
  long long int start[2], end[2], pivotInd;
  int * pivot = new int[D];
  int i, tmp[D];
  float tmpFunVal[O];
	if (first < last) {
		start[1] = first;
		end[0] = last;
		pivotInd = (first + last) / 2;
		memcpy(pivot, v+(pivotInd*D), cellLen);
		while (start[1] <= end[0]) {

			while (cmp(pivot, v+start[1]*D)>0)
				start[1]++;
			while (cmp(v + end[0]*D, pivot)>0)
				end[0]--;
			if (start[1] <= end[0]) {
				
				memcpy(tmp,v+ D*start[1], cellLen);
				memcpy(v+ D*start[1],v+ D*end[0] , cellLen);
				memcpy(v+ D*end[0], tmp , cellLen);
				
				memcpy(tmpFunVal, funEval + O* start[1], O*sizeof(float));
				memcpy(funEval + O* start[1], funEval + O*end[0], O*sizeof(float));
				memcpy(funEval + O*end[0], tmpFunVal, O*sizeof(float));
				
				start[1]++;
				end[0]--;
			}
		}
		start[0] = first; 
		end[1]   = last; 
		#pragma omp parallel 
		{
			#pragma omp for nowait
			for(i = 0; i <= 1; i++) {
				qs(v, start[i], end[i], funEval);
			}
		}
	}
	delete[] pivot;
}
void qs(int *v, long long int first,long long int last, float * funEval, bool * flag) {
  long long int start[2], end[2], pivotInd;
  int * pivot = new int[D];
  int i, tmp[D];
  float tmpFunVal[O];
  bool tmpFlag;
	if (first < last) {
		start[1] = first;
		end[0] = last;
		pivotInd = (first + last) / 2;
		memcpy(pivot, v+(pivotInd*D), cellLen);
		while (start[1] <= end[0]) {

			while (cmp(pivot, v+start[1]*D)>0)
				start[1]++;
			while (cmp(v + end[0]*D, pivot)>0)
				end[0]--;
			if (start[1] <= end[0]) {
				
				memcpy(tmp,v+ D*start[1], cellLen);
				memcpy(v+ D*start[1],v+ D*end[0] , cellLen);
				memcpy(v+ D*end[0], tmp , cellLen);
				
				memcpy(tmpFunVal, funEval + O* start[1], O*sizeof(float));
				memcpy(funEval + O* start[1], funEval + O*end[0], O*sizeof(float));
				memcpy(funEval + O*end[0], tmpFunVal, O*sizeof(float));
				
				tmpFlag = flag[start[1]];
				flag[start[1]] = flag[end[0]];
				flag[end[0]] = tmpFlag;
				
				start[1]++;
				end[0]--;
			}
		}
		start[0] = first; 
		end[1]   = last; 
		#pragma omp parallel 
		{
			#pragma omp for nowait
			for(i = 0; i <= 1; i++) {
				qs(v, start[i], end[i], funEval, flag);
			}
		}
	}
	delete[] pivot;
}
////////////////////////////////////////////////////////////////////////////////
long long int binSrch(int * arr, int * x, long long int maxInd) {
	long long int first, last, middle, Ind = -1;
	if (maxInd == (long long int)-1 || cmp(arr,x)>0 || cmp(x,arr+maxInd*D)>0) {
		return -1;
	}
	else {
		first = 0;
		last = maxInd;
		middle = (first+last)/2;
		while (first <= last) {
			if (cmp(x,arr+middle*D)>0)
				first = middle + 1;    
			else if (cmp(x,arr+middle*D)==0)
				return middle;
			else
				last = middle - 1;
 
			middle = (first + last)/2;
		}
		if (first > last)
			return -1;   
	}
}
////////////////////////////////////////////////////////////////////////////////
// This function finds the orthogonal neighbors of a certain cell
void neighFind(int * neighbors, int * cell, int * N) {
	int * z = new int[D];
	for (int i = 0; i  < D; i++) {
		memcpy(z,cell,cellLen);//dest,source, amount
		z[i] = z[i]-1;
		if (z[i] < 0 || z[i] >= N[i]) {
			memset(neighbors+2*i*D,-1,cellLen);
		}
		else {
			memcpy(neighbors+2*i*D, z, cellLen);
		}
		z[i] = z[i]+2;
		if (z[i] < 0 || z[i] >= N[i]) {
			memset(neighbors+(2*i+1)*D,-1,cellLen);
		}
		else {
			memcpy(neighbors+(2*i+1)*D, z, cellLen);
		}
	}
	delete[] z;
}
////////////////////////////////////////////////////////////////////////////////
//This function exclude specific members from a list that are marked by a flag. The input can be a list of cell (with or without correcponding FunEval) and a list of flags. The program removes the items carrying a True flag from the list
// Input:	list:		The list of cells to be truncated
//			listInd:	The list length
//			flags:		The flags. if flag[i]=true => list[i] should be removed
// Output:	list:		The updated list (shorter)
//			listInd:	Length of the updated list
void listShrink (int * list, long long int & listInd, bool * flags) {
long long int tmp = listInd, flagInd = 0;

for (long long int i = 0; i < tmp; i++) {
	if (flags[i] == 0) {
		memmove(list+flagInd*D,list+i*D,cellLen);
		flagInd++;
	}
	else
		listInd--;
}
//memset(list+(listInd+1)*D, -1, (maxTempLen-listInd-1)*cellLen);
}
// Overloaded function for accepting the funEvaluations as well.
void listShrink (int * list, long long int & listInd, float * funEvals, bool * flags) {
long long int tmp = listInd, flagInd = 0;

for (long long int i = 0; i < tmp; i++) {
	if (flags[i] == 0) {
		memmove(list+flagInd*D,list+i*D,cellLen);
		memmove(funEvals+flagInd*O,funEvals+i*O,O*sizeof(float));
		flagInd++;
	}
	else
		listInd--;
}

//memset(list+(listInd+1)*D, -1, (maxTempLen-listInd-1)*cellLen);
//memset(funEvals+(listInd+1)*O, 0, (maxTempLen-listInd-1)*O*sizeof(float));

}

////////////////////////////////////////////////////////////////////////////////
void insert(int * cellsGlo,long long int & cellsGloInd, float * funEvalGlo, bool * cellsGloFlag, int * cells, long long int cellsLen, float * funEval) {
	// steps: 	1) sort the cells
	//			2) flag repetitions and the (-1)s and the ones existing in cellsGlo
	//			3) shrink cells to remove flagged cells
	//			4) put it after the cellsGlo with funEvals when shrinked
	//			5) qs the cellsGlo with the funEvals and flags
	long long int i=0, ind=0;
	bool * remFlag = new bool [cellsLen];
	memset(remFlag,0,cellsLen*sizeof(bool));
	int * tmpCell = new int[D];
	memset(tmpCell, -1, cellLen);
	int * cellsTmp = new int[cellsLen*D];
	memcpy(cellsTmp, cells, cellsLen*D*sizeof(int));
	float * funEvalTmp = new float[cellsLen*O];
	memcpy(funEvalTmp, funEval, cellsLen*O*sizeof(float));
	
	//qs(cells, 0, cellsLen-1, funEval);
	qs(cellsTmp, 0, cellsLen-1, funEvalTmp);
	for (i = 0; i < cellsLen; i++) {
		//flag the cells to be removed in here
		if (cellsTmp[i*D]==-1 || cmp(cellsTmp+i*D, tmpCell)==0) {
			remFlag[i] = true;
		}
		else {
			memcpy(tmpCell, cellsTmp+i*D, cellLen);
			ind = binSrch(cellsGlo, tmpCell, cellsGloInd);
			if (ind != (long long int)-1) remFlag[i] = true;
		}
	}
	
	//for (i = 0; i < cellsLen; i++) {
	//	printf("i=%4lld, cell=%3d,%3d,%3d, flag=%d\n",i,cells[i*D],cells[i*D+1],cells[i*D+2], remFlag[i]);
	//}
	// Shrink here
	listShrink (cellsTmp, cellsLen, funEvalTmp, remFlag);
	
	memcpy(cellsGlo+(cellsGloInd+1)*D,cellsTmp, cellsLen*cellLen);
	memcpy(funEvalGlo+O*(cellsGloInd+1),funEvalTmp,cellsLen*O*sizeof(float));
	memset(cellsGloFlag+cellsGloInd+1, 0, cellsLen * sizeof(bool));
	
	cellsGloInd += cellsLen;
	
	qs(cellsGlo, 0, cellsGloInd, funEvalGlo, cellsGloFlag);
	delete[] remFlag;
	delete[] tmpCell;
	delete[] cellsTmp;
	delete[] funEvalTmp;
}
////////////////////////////////////////////////////////////////////////////////
// This functions is used to find the objective values for one or a set of cells
// Input:	cells 		The cell numbers to be evaluated
// 			nCells		the length of cells (%lld)
//			funEvalGlo	function evaluation database
//			cellsGlo	the cell numbers for function evaluation database
//			cellGloInd	the length of cellsGlo
// Output: 	function evals	The corresponding objective values for the 'cells' 
void funEval(float * functionEvals, int * cells, long long int nCells, int * cellsGlo, float * funEvalGlo, long long int cellsGloInd, float * lb, float * h){
	long long int insInd;
	float y[O], x[D];
	for (long long int i = 0; i < nCells; i++) {
		insInd = binSrch(cellsGlo, cells+i*D, cellsGloInd);
		if (cells[i*D] == -1) { // not a feasble neighbor
			memset(functionEvals+i*O, -1, O*sizeof(float));
		}
		else {
			if (insInd >= 0) { //the cell already existed in the list
				memcpy(functionEvals+i*O,funEvalGlo+insInd*O,O*sizeof(float));
				//printf("cell [%3d,%3d,%3d] found at index: %lld, funEval=%f,%f,%f\n",cells[i*D],cells[i*D+1],cells[i*D+2],insInd,functionEvals[i*O], functionEvals[i*O+1],functionEvals[i*O+2]);
			}
			else {
				z2x( x, cells+i*D, lb, h);
				f(y,x); 
				memcpy(functionEvals+i*O,y,O*sizeof(float));
			}
		}
	}
	//delete[] cell;
} 
////////////////////////////////////////////////////////////////////////////////
void initialize(int * todoCell, long long int & todoLen, int * gloCell, float * gloCellObj, bool * gloCellFlag, long long int & gloCellInd, int initCellLen, float * lb, float * h) {
	int z[D];
	float x[D],y[O];
	int * tmpCell = new int[initCellLen*D];
	float * tmpCellObj = new float[initCellLen*O];
	for (int i = 0; i < initCellLen; i++) {
		for (int k = 0; k < D; k++) {
			z[k] = (int)rnd(0,N[k]);
		}	
		z2x(x,z,lb, h);
		f(y,x);
		memcpy(tmpCell+i*D,z,D*sizeof(int));
		memcpy(tmpCellObj+i*O,y,O*sizeof(float));
		//
	}
	insert(gloCell, gloCellInd, gloCellObj, gloCellFlag, tmpCell, (long long int)initCellLen, tmpCellObj);
	memcpy(todoCell, gloCell, (gloCellInd+1)*D*sizeof(int));
	todoLen = gloCellInd+1;
	delete[] tmpCell;
	delete[] tmpCellObj;
}
////////////////////////////////////////////////////////////////////////////////
void initFromFile(int * todoCell, long long int & todoLen, 
			int * gloCell, float * gloCellObj, bool * gloCellFlag, 
			long long int & gloCellInd, float * lb, float * h, int * N,
			std::string fcn, char * loopCntStr,bool log_gen,
			std::ofstream& log)
{
	string fName;
	int len = 15; /*lenght of each # in the input file plus the space after it */
	long long int todoCnt = 0, fileCnt = 0, i;
	//char cLine[(len+2)*D];
	//char valC[len+1];
	bool * remFlag = new bool[maxTempLen];
	int tmpCell[D];
	//valC[len] = '\0';
	float valF[D];
	ifstream file;
	fName = "./" + fcn + "/"  + fcn + "_coord.txt";
	file.open(fName.c_str());
	string value;
	if(!file) {
    	cout << "The file '(Problem name)_coord.txt' does not exist." << endl;
    } 
	else {
		while (getline(file, value)) {
			std::stringstream ss(value); // must #include <sstream>
			fileCnt++;
			//cout << "the Line:" << value << endl;
			for (i = 0; i<D; i++) {
				//strncpy(valC, cLine+i*len, len);
				//valF[i] = (float)atof(valC);
				ss >> valF[i];
				//printf("%8.2e,",valF[i]);
			}	
			//printf("\n");
			x2z(todoCell + todoCnt*D, valF, lb, h, N);
			//printf("==[%d,%d]\n",todoCell[todoCnt*D],todoCell[todoCnt*D+1]);
			todoCnt++;
		}
		qs(todoCell, 0, todoCnt-1);
		memset(tmpCell, -1, D*sizeof(int));
		for (i = 0; i < todoCnt; i++) {
			if (todoCell[i*D]==-1 || cmp(todoCell+i*D, tmpCell)==0) {
				remFlag[i] = true;
			}
			else {
				memcpy(tmpCell, todoCell+i*D, D*sizeof(int));
			}
		}
		listShrink (todoCell, todoCnt, remFlag);
		printf("Read from file completed!\n");
		if (log_gen) log << "Read from file completed. " << fileCnt <<" GA point resulted in "
						<< todoCnt << " initial cells." << endl;
		cout << "Read from file completed. " << fileCnt <<" GA point resulted in "
			<< todoCnt << " initial cells." << endl;
	}	
	todoLen = todoCnt;
	file.close();
	delete[] remFlag;
	
}

/*
void initFromFile(int * todoCell, long long int & todoLen, 
			int * gloCell, float * gloCellObj, bool * gloCellFlag, 
			long long int & gloCellInd, float * lb, float * h, int * N,
			std::string fcn, char * loopCntStr,bool log_gen,
			std::ofstream& log)
{
	string fName;
	int len = 15; //lenght of each # in the input file plus the space after it 
	long long int todoCnt = 0, fileCnt = 0, i;
	char cLine[(len+2)*D];
	char valC[len+1];
	bool * remFlag = new bool[maxTempLen];
	int tmpCell[D];
	valC[len] = '\0';
	float valF[D];
	fName = "./" + fcn + "/"  + fcn + "_coord.txt";
	FILE * coord;
	coord = fopen(fName.c_str(),"r");
	if (coord == (FILE*)NULL)
	{
		cout << "The file '(Problem name)_coord.txt' does not exist." << endl;
	}
	else {
		//todoLen = 0;
		while (!feof(coord)) {
			fgets(cLine, (len+2)*D, coord);
			fileCnt++;
			//cout << cLine << endl;
			for (i = 0; i<D; i++) {
				strncpy(valC, cLine+i*len, len);
				valF[i] = (float)atof(valC);
				//printf("%8.2e,",valF[i]);
			}	
			//printf("\n");
			x2z(todoCell + todoCnt*D, valF, lb, h, N);
			//printf("==[%d,%d]\n",todoCell[todoCnt*D],todoCell[todoCnt*D+1]);
			todoCnt++;
		}
		qs(todoCell, 0, todoCnt-1);
		memset(tmpCell, -1, D*sizeof(int));
		for (i = 0; i < todoCnt; i++) {
			if (todoCell[i*D]==-1 || cmp(todoCell+i*D, tmpCell)==0) {
				remFlag[i] = true;
			}
			else {
				memcpy(tmpCell, todoCell+i*D, D*sizeof(int));
			}
		}
		listShrink (todoCell, todoCnt, remFlag);
		printf("Read from file completed!\n");
		if (log_gen) log << "Read from file completed. " << fileCnt <<" GA point resulted in "
						<< todoCnt << " initial cells." << endl;
		cout << "Read from file completed. " << fileCnt <<" GA point resulted in "
			<< todoCnt << " initial cells." << endl;
	}	
	todoLen = todoCnt;
	fclose(coord);
	delete[] remFlag;
	
}
*/
////////////////////////////////////////////////////////////////////////////////
//This function performs a dominance check on objective function values for two specific points
// Input: Objective function values for two points (f(x),f(y))
// Output: -1 if the input1 is dominated by input2: f(x)>f(y)
//			0 if the are non-dominant
//			1 if input one dominates input2 f(x)<f(y)
//			2 if equal
int domChkSingle(float * x, float * y) {
	bool dominant, dominated, equal;
	dominant = dominated = equal = true;
	for (int k = 0; k < O; k++) { // for all the objectives
		if (y[k]-x[k]>1e-5) {//x[k] < y[k]
			dominated = false;
			equal = false;
		}
		else if (x[k]-y[k]>1e-5) {//y[k] < x[k]
			dominant = false;
			equal = false;
		}
	}
	if (dominated == true && equal == false) {
		// the cell is dominated by input2 (input2 is better)
		return -1;
	}
	else if (dominant == true && equal == false) 
		// input1 is dominant to input2 (input1 is better)
		return 1;
	else if (equal == true)
		return 2;
	else
		return 0;
}
////////////////////////////////////////////////////////////////////////////////
// This function finds the descent slope for 2 adjacent neighbors

float distance(int * cell1, float * y1, int * cell2, float * y2, int * N, float * h, float * lb) {
	float x1[D],x2[D], dx, dy, sum = 0;
	z2x(x1, cell1,lb, h);
	z2x(x2, cell2,lb, h);
	for (int i = 0; i < D; i++) {
		sum += pow(x1[i]-x2[i],2);
	}
	dx = pow(sum,0.5);
	sum = 0;
	for (int i = 0; i < O; i++) {
		sum += pow(y1[i]-y2[i],2);
	}
	dy = pow(sum,0.5);
	return dy/dx;
}
////////////////////////////////////////////////////////////////////////////////
void explore(bool * destCellFlag, bool * solFlag, int * cell, float * cellObj, int * neiCell, float * neiCellObj, int * N, float * h, float * lb) {
	float goodness = 0, dist;
	//int dest[D];								/*Destination Cell*/
	//memset(dest, -1, D*sizeof(int));
	int dest = -1;
	int domVal[2*D];							/*Dominancy comparison value*/
	int i;
	for (i = 0; i < 2*D; i++) {
		if (neiCell[i*D] != -1) {//if the neighbor exists
			domVal[i] = domChkSingle(cellObj, neiCellObj+i*O);
			if (domVal[i] == -1) { 				/*better neighbor found*/
				dist = distance(cell, cellObj, neiCell+i*D, neiCellObj+i*O, N, h, lb);
				if (dist > goodness || (dist == goodness && rnd()>0.5)) { //the best neighbor so far
					goodness = dist;
					//memcpy(dest, neiCell+i*D, D*sizeof(int));
					dest = i;
				}
			}
		}
		else domVal[i] = 1;
	}
	if (goodness == 0) { // no better neighbor is found. We might be on the solution!
		solFlag[0] = true; //make sure to pass the solFlag+i to this function
		for (i = 0; i < 2*D; i++) {
			if ((domVal[i] == 0 || domVal[i] == 2) && neiCell[i*D] != -1) {
				//z2x(x, cell,  lb,  h);
				//printf("x=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14] );
				destCellFlag[i] = true;
				//this part is only for testing
				
				//printf("cen=[%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d],[%5.2f,%5.2f,%5.2f]\nnei=[%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d],[%5.2f,%5.2f,%5.2f]\n"				,cell[0],cell[1],cell[2],cell[3],cell[4],cell[5],cell[6],cell[7],cell[8],cell[9],cell[10],cell[11],cell[12],cell[13],cell[14],cellObj[0], cellObj[1],cellObj[2], neiCell[i*D],neiCell[i*D+1],neiCell[i*D+2],neiCell[i*D+3],neiCell[i*D+4],neiCell[i*D+5],neiCell[i*D+6],neiCell[i*D+7],neiCell[i*D+8],neiCell[i*D+9],neiCell[i*D+10],neiCell[i*D+11],neiCell[i*D+12],neiCell[i*D+13],neiCell[i*D+14],neiCellObj[i*O],neiCellObj[i*O+1],neiCellObj[i*O+2]);
			}
		}
	}
	else {
		destCellFlag[dest] = true;
	}
}
////////////////////////////////////////////////////////////////////////////////
// This fucntion receives a single cell or an array (list) and then does the dominancy check with the solCells and insert it to the solution if needed. It will keep the solCells always updated with the dominancy check
void domChkIns(int * solCells, float * solFunEval, long long int & solCellsInd, int * list, float * listFunEval, long long int & listInd) {
bool * remSolFlag = new bool[solCellsInd+1]; //A shared memory
memset(remSolFlag,0,(solCellsInd+1)*sizeof(bool));
bool * listFlag = new bool[listInd+1]; //A shared memory array
memset(listFlag,0,(listInd+1)*sizeof(bool));
long long int cellInd, i, j;
int k;
float cellFunValue[O];
float otherCellFunValue[O];
int domVal;

// we need to first check the list and make sure it is a non-dominating set. Then we can compare it to the solCells
#pragma omp parallel default(none) private(i, j, k, cellFunValue, otherCellFunValue, domVal) shared(list, listInd, listFunEval, listFlag, remSolFlag, solCells, solCellsInd, solFunEval)
{
	#pragma omp for
	for (i = 0; i < listInd; i++) {
		//listFlag value will be 1 for the cells that should be removed from the set
		memcpy(cellFunValue, listFunEval+i*O, O*sizeof(float));
		for (j = 0; j < listInd; j++) {
			memcpy(otherCellFunValue, listFunEval+j*O, O*sizeof(float));
			// Now that the FunValues are available, we can compare for dominancy
			domVal = domChkSingle(cellFunValue, otherCellFunValue);
			if (domVal == -1) {// input 1 is dominated
				listFlag[i] = true; // list[i] should be removed
				break; // no need to compare this one with the rest since the other cell is better
			}
			else if (domVal == 1) //input one dominates
				listFlag[j] = true; // list j should be removed
			else if (domVal == 2 && i<j && list[i] == list[j]) {//equal case //if the cell #s are equal,one should be removed
				listFlag[i] = true;
			}
		}
	}
}

//should run serially

listShrink (list, listInd, listFunEval, listFlag);

memset(listFlag,0,listInd*sizeof(bool));

#pragma omp parallel default(none) private(i, j, k, cellFunValue, otherCellFunValue, domVal) shared(list, listInd, listFunEval, listFlag, remSolFlag, solCells, solCellsInd, solFunEval)
{
	// DO this loop in parallel
	// in this parallel loop we compare each member in the list with all the members in the solFunEval. There are flags for the solCells to be remoeved because they are dominated by list and there are flags for the list cells that need to be inserted to the solCells
	#pragma omp for
	for (long long int i = 0; i < listInd; i++) { // for all the cells in the list
	//	memset(remSolFlagTemp,0,solCellsInd*sizeof(bool));
		memcpy(cellFunValue, listFunEval+i*O, O*sizeof(float));
		for (long long int j = 0; j < solCellsInd; j++) { // compare to sol cells
			memcpy(otherCellFunValue, solFunEval+j*O, O*sizeof(float));
			// Now that the FunValues are available, we can compare for dominancy
			domVal = domChkSingle(cellFunValue, otherCellFunValue);
			if (domVal == -1) {// input 1 is dominated
				listFlag[i] = true; // list[i] should be removed
				break;
			}
			else if (domVal == 1) {//input one dominates
				remSolFlag[j] = true; // list j should be removed
				//printf("this case!\n");
			}
			else if (domVal == 2 && list[i] == solCells[j]) {//equal case //if the cells are equal,one should be removed
				listFlag[i] = true;
			}
		}
	}
}

// now that we have all the remove flags on the both lists (solCells and list) we should update and shrink the lists and then insert the list into the sol cells
	listShrink (list, listInd, listFunEval, listFlag);
	listShrink (solCells, solCellsInd, solFunEval, remSolFlag);



	// it is now time to insert list into the solCells. This call is NOT parallel
	memcpy(solCells+solCellsInd*D, list, listInd*cellLen);
	memcpy(solFunEval+solCellsInd*O, listFunEval, listInd*O*sizeof(float));
	solCellsInd += listInd;



delete[] listFlag;
delete[] remSolFlag;
}
////////////////////////////////////////////////////////////////////////////////
void updateSol(int * solCell, 
				float * solCellObj, 
				long long int & solCellInd,
				int * tmpSolCell, 
				float * tmpSolCellObj,
				int * todoCell,
				long long int todoLen,
				bool * solFlag,
				int * gloCell,
				long long int gloCellInd,
				float * gloCellObj
	) {
	long long int i, tmpSolCnt=0, ind;
	float * obj = new float [O];
	for (i = 0; i < todoLen; i++) {
		if (solFlag[i] == true) { 			/*that cell might belong to sol*/
			memcpy(tmpSolCell+tmpSolCnt*D, todoCell+i*D, D*sizeof(int));
			ind = binSrch(gloCell, tmpSolCell+tmpSolCnt*D, gloCellInd);
			if (ind == -1) {
				printf("Update Sol: binary search error: cell not found in the gloCell\n");
				memset(tmpSolCellObj+tmpSolCnt*O, 0, O*sizeof(float));
			}
			else {
				memcpy(tmpSolCellObj+tmpSolCnt*O, gloCellObj+ind*O, O*sizeof(float));
			}
			tmpSolCnt++;
		}
	}
	//printf("%lld cells on the tmpSol list\n",tmpSolCnt);
	domChkIns(solCell, solCellObj, solCellInd, tmpSolCell, tmpSolCellObj, tmpSolCnt);
	
	
	delete[] obj;
}
////////////////////////////////////////////////////////////////////////////////
void updateToDoCell(int * todoCell, 
					long long int & todoLen, 
					bool * destCellFlag,
					int * neiCell,
					int * gloCell,
					long long int & gloCellInd,
					bool * gloCellFlag
	) {
	long long int i, todoCellCnt = 0, ind;
	for (i = 0; i < todoLen*2*D; i++) {
			//printf("Now working on neighbor:[%2d,%2d,%2d],destCellFlag=%d\n",neiCell[i*D],neiCell[i*D+1],neiCell[i*D+2],destCellFlag[i]);
		if (destCellFlag[i] == true) {
			
			ind = binSrch(gloCell, neiCell+i*D, gloCellInd);
			if (ind != -1 && gloCellFlag[ind] == false) {
				memcpy(todoCell+todoCellCnt*D, neiCell+i*D, D*sizeof(int));
				gloCellFlag[ind] = true;
				todoCellCnt++;
				if (todoCellCnt >= maxTempLen) { 		/*The list got full*/
					printf("Warning: Max length of the todoCell reached\n");
				}
			}
			else if (ind == -1) {
				printf("updateToDoCell Warning: the nei cell not found in the gloCell\n");
			}
		}
	}
	todoLen = todoCellCnt;
}
////////////////////////////////////////////////////////////////////////////////
void refine(int * cells, long long int cellsInd, int * N, int * div, int divInd) {
// we can do this in parallel
long long int i;
int * z = new int[D];
int Nnew[D], division;
division = div[divInd] * 2 + 1;
memcpy(Nnew,N,D*sizeof(int));
Nnew[divInd] *= division;
//#pragma omp parallel default(none) private(i, cell, z) shared(cells, N, Nnew, divInd, cellsInd, division)
{
	//#pragma omp for
	for (i = 0; i <= cellsInd; i++) {
		//cell = cells[i];
		//memcpy(z, cells+i*D, cellLen);
		//c2z(z, cell, N);
		//z[divInd] = z[divInd] * division + division/2;
		cells[i*D+divInd] = cells[i*D+divInd] * division + division/2; 
		//cell = z2c(z, Nnew);
		//cells[i] = cell;
	}
}
delete[] z;
//N[divInd] *= division; should be done in the main code
}
////////////////////////////////////////////////////////////////////////////////
void subDiv(int * todoCell, long long int & todoLen,
			int * gloCell, long long int & gloCellInd,
			int * solCell, long long int & solCellInd, float * solCellObj,
			int * N, float * h, float * ub, float * lb,
			int & subDivInd, int * subDivSize) {
	long long int i;
	int div = subDivSize[subDivInd]*2 + 1;
	for (i = 0; i <= gloCellInd; i++) {
		gloCell[i*D+subDivInd] = gloCell[i*D+subDivInd] * div + div/2; 
	}
	for (i = 0; i < solCellInd; i++) {
		solCell[i*D+subDivInd] = solCell[i*D+subDivInd] * div + div/2; 
	}
	// the new todoCell is made of the solCell
	memset(todoCell, -1, D*maxTempLen*sizeof(int));
	memcpy(todoCell, solCell, D*solCellInd*sizeof(int));
	todoLen = solCellInd;
	
	memset(solCell, -1, D*maxArrLen*sizeof(int));
	memset(solCellObj, 0, O*maxArrLen*sizeof(float));
	solCellInd = 0;
	
	N[subDivInd] *= div;
	h[subDivInd] = (ub[subDivInd]-lb[subDivInd])/((float)N[subDivInd]);
}
////////////////////////////////////////////////////////////////////////////////
void exportSol(std::string fcn,char * loopCntStr, int * solCell, long long int solCellInd, float * solCellObj, float * lb, float * h, int & subDivInd) {

	std::string filename;
	float * x = new float[D];
	ofstream ps, pf;
	long long int i, j;
	char str[3];
	sprintf(str,"%d",subDivInd);
	filename = "./" + fcn + "/" + fcn + "_PS_sub" + str + ".dat";
	ps.open(filename.c_str());
	filename = "./" + fcn + "/" + fcn + "_PF_sub" + str + ".dat";
	pf.open(filename.c_str());
	for (i = 0; i < solCellInd; i++) {
		//c2x(x, solCells[i], N, h, lb);
		z2x(x, solCell+i*D, lb, h);
		for (j = 0 ; j < D; j++) ps << setw(15) << scientific << x[j];
		ps << endl;
		for (j = 0 ; j < O; j++) pf << setw(15) << scientific << solCellObj[i*O+j];
		pf << endl;
	}
	ps.close();
	pf.close();
	delete[] x;
}

