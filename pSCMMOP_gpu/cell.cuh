#include <math.h>
#include "cuda_runtime.h"

#define O 2 // number of output
#define D 5 // number of input
#define pop 243 // number of neighbours (3^D)
#define C 1 // number of constrains

//MOP
__device__ void f(float * r,float * x){
	/*
	//ZDT-2
	r[0]=(powf(x[0]-1,4)+powf(x[1]-1,2))*(.5f)+(powf(x[0]-1,2)+powf(x[1]-1,2))*.5f;
	r[1]=(powf(x[0]+1,2)+powf(x[1]+1,2))*(.5f)+(powf(x[0]+1,2)+powf(x[1]+1,2))*.5f;
	*/

	////twocircles
	//r[0]=powf(x[0],2)+powf(x[1],2);
	//r[1]=powf(x[0]-10,2)+powf(x[1],2);

	//2obj
	int i;
	float w[D],fun;
	r[0]=0;
	for(i=0;i<D;i++)
		r[0]+=x[i];
	for(i=0;i<D;i++){
		if(i==0 || i==1)
			w[i]=0.01f*expf(-powf(x[i]/20,2.5f));
		else
			w[i]=0.01f*expf(-x[i]/15);
	}
	fun=1;
	for(i=0;i<D;i++)
		fun*=1-w[i];
	r[1] = 1-fun;

	////threecircles
	//r[0] = powf(x[0]+5,2)+powf(x[1],2);
	//r[1] = powf(x[0]-5,2)+powf(x[1],2);
	//r[2] = powf(x[0],2)+powf(x[1]-5,2);

}

//Constraints
__device__ void g(float* r,float *x){
	//set r[0]=-1 if there is no constraints
	//default inequalities are less than zero

	r[0]=-1;
	
	//r[0]=(powf(x[0]-5,2)+powf(x[1],2)-5); //inside
	//r[0]=-(powf(x[0]-5,2)+powf(x[1],2)-5); //outside
}

/*
__device__ int roundf(float x){
	int y;
	if(floorf(x)==floorf(x+0.5f))
		y = x;
	else
		y = floorf(x) + 1;
	return y;
}
*/

__device__ __host__ int ztocell(const int * z, const int * N){
	int i;
	int ncell;
	int temp[D];

	for(i=0;i!=D;i++)
		temp[i]=z[i];
	ncell = temp[0];
	
	int b = N[0];
	for(i=0;i!=D;i++)
		temp[i]-=1;
	for(i=1;i!=D;i++){
		ncell+=temp[i]*b;
		b*=N[i];
	}
	//ncell+=1;
	return ncell;
}

__device__ void xtoz(int * z, const float *x, const float *h, const float * lb){
	int i;
	for(i=0;i!=D;i++)
		z[i] = roundf((x[i]-lb[i])/h[i] + 0.5f);
}

__device__ void ztox(float * xd, const int * z, const float * h, const float * lb){
	int i;
	for(i=0;i!=D;i++)
		xd[i] = lb[i]+h[i]*z[i]-0.5f*h[i];
}

__device__ __host__ void celltoz(int * z, const int cs, const int * N){
	int coord[D];
	int i;
	int cell=cs;
	cell-=1;
	for(i=0;i!=D;i++){
		coord[i] = cell%N[i] + 1;
		cell = floorf(float(cell)/float(N[i]));
		z[i] = coord[i];
	}
}

__device__ __host__ void neighbour_finder(int * neighbour, const int * N, int cs){
	/*
	Input arguments:
	neighbour:    intialized arary storing neighbour cells
	N:            cell space partition
	cs:           central cell number
	pop:          maximum number of possible neighbours (dim^3)

	Output arguments:
	neighbour:    neighbour cells with zero as invalid cells
	*/
	int z_center[D];
	int z_current[D];
	int i;
	int j;
	bool goodcell=true;
	celltoz(z_center,cs,N);
	for(i=0;i!=D;i++)
		z_current[i]=z_center[i]-1;
	for(i=0;i!=pop;i++){
		for(j=0;j!=D;j++){
			if(z_current[j]>N[j] || z_current[j]<1){
				goodcell=false;
				break;
			}
		}
		if(goodcell){
			neighbour[i] = ztocell(z_current,N);
		}
		else
			neighbour[i]=0;
		goodcell=true;
		z_current[0]+=1;
		for(j=0;j!=D-1;j++){
			if(z_current[j]>z_center[j]+1){
				z_current[j]=z_center[j]-1;
				z_current[j+1]=z_current[j+1]+1;
			}
		}
	}
	//rule out the central cell
	for(i=0;i!=pop;i++){
		if(neighbour[i]==cs)
			neighbour[i]=0;
	}
}

__device__ __host__ void cartprod(int * X, const int * div){
	//cartesin product of subdivision vectors (-div[i]:div[i])
	/*
	Input arguments:
	div:          subdivision vector with (2*div+1) division on each cell
	sizeThisSet:  number of each dimension's subcells (2*div+1)

	Output arguments:
	X:            allcomb of integer coordinate incrementation in 1D array
	*/

	int i,j,dim;
	int ixVect[D],
		sizeThisSet[D];

	dim = 1;
	for(i=0;i<D;i++){
		sizeThisSet[i] = 2*div[i]+1;
		dim*=sizeThisSet[i];
	}

	for(i=0;i<dim;i++){
		celltoz(ixVect,i+1,sizeThisSet);
		for(j=0;j<D;j++){
			X[i*D+j] = -div[j]+ixVect[j]-1;
		}
	}
}

__device__ __host__ void subdivision(int * rcells, const int cell, const int * N, int * N_new, 
									 const int dim, const int * div, const int * allcomb){
	//subdivide the given cell set with 2*div+1 smaller cells
	/*
	Input arguments:
	cell:      single coarse cell remain to be refined
	N:         coarse cell space partition
	div:       subdivision times (2*div+1)
	dim:       number of subdivision cells in a cell (2*div+1)^N

	Output arguments:
	N:         newly refined cell space partition
	rcells:    cell number of smaller cells
	*/
	int i,j,k;
	int cs,cr;
	int zcs[D],zinc[D];

	// refinement partition
	for(i=0;i<D;i++)
		N_new[i]= N[i]*(2*div[i]+1);

	cs = cell;
	celltoz(zcs,cs,N);
	for(j=0;j<D;j++)
		zcs[j] = zcs[j]*(2*div[j]+1)-div[j]; // new central cell coord
		
	for(j=0;j<dim;j++){
		for(k=0;k<D;k++)
			zinc[k] = zcs[k]+allcomb[j*D+k];
		cr = ztocell(zinc,N_new); // newly refined cell
		rcells[j] = cr;
	}
}