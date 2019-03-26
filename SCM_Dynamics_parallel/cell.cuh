#include <math.h>
#include "cuda_runtime.h"

#define D 2 // number of input
#define O 2 // number of output
#define SamNum 15 //number of random sampling points within each cell
#define pop 27 //3^D number of neighbours

//test problems of the potential field from Rn->R1
__device__ __host__ float fp(float *x){
	////2d potential field
	//float A[4] = {-200.f,-100.f,-170.f,-15.f};
	//float a[4] = {-1.f,-1.f,-6.5f,-0.7f};
	//float b[4] = {0.f,0.f,11.f,0.6f};
	//float c[4] = {-10.f,-10.f,-6.5f,0.7f};
	//float x0[4] = {1.f,0.f,-0.5f,-1.f};
	//float y0[4] = {0.f,0.5f,1.5f,1.f};
	//float sum = 0.f;
	//for(int i=0;i<4;i++)
	//	sum+=A[i]*expf(a[i]*powf(x[0]-x0[i],2.f)+b[i]*(x[0]-x0[i])*(x[1]-y0[i])+c[i]*powf(x[1]-y0[i],2.f));
	//return sum;

	//3d potential field
	float v12,v13,v23;
	float r12,r13,r23;
	float x2,x3,y3;
	x2 = x[0]; x3 = x[1]; y3 = x[2];
	r12 = x2;
	r13 = sqrtf(x3*x3+y3*y3);
	r23 = sqrtf((x3-x2)*(x3-x2)+y3*y3);
	v12 = powf(1.f/r12,12.f)-2.f*powf(1.f/r12,6.f);
	v13 = powf(1.f/r13,12.f)-2.f*powf(1.f/r13,6.f);
	v23 = powf(1.f/r23,12.f)-2.f*powf(1.f/r23,6.f);
	return v12+v13+v23;

	////2d test function for saddle
	//return sinf(0.5f*x[0]*x[0]-.25f*powf(x[1],3.f)+3.f)*cosf(2.f*x[0]+1.f-expf(x[1]));
}

//test problems of f(x)=0 with Rn->Rn
__device__ __host__ void f(float * r,const float * x){
	////2d 9 point
	//r[0] = 4*x[0]*(powf(x[0],2)+x[1]-11)+2*(x[0]+powf(x[1],2)-7);
	//r[1] = 2*(powf(x[0],2)+x[1]-11)+4*x[1]*(x[0]+powf(x[1],2)-7);

	//2d gradient of Beale's function
	r[0] = (1.5f-x[0]+x[0]*x[1])*(-1.f+x[1])+(2.25f-x[0]+x[0]*x[1]*x[1])*
		(-1.f+x[1]*x[1])+(2.625f-x[0]+x[0]*powf(x[1],3.f))*(-1.f+powf(x[1],3.f));
	r[1] = (1.5f-x[0]+x[0]*x[1])*x[0]+(2.25f-x[0]+x[0]*x[1]*x[1])*2.f*x[0]*x[1]+
		(2.625f-x[0]+x[0]*powf(x[1],3))*3.f*x[0]*x[1]*x[1];

	////2d potential field
	//float A[4] = {-200.f,-100.f,-170.f,-15.f};
	//float a[4] = {-1.f,-1.f,-6.5f,-0.7f};
	//float b[4] = {0.f,0.f,11.f,0.6f};
	//float c[4] = {-10.f,-10.f,-6.5f,0.7f};
	//float x0[4] = {1.f,0.f,-0.5f,-1.f};
	//float y0[4] = {0.f,0.5f,1.5f,1.f};
	//float sum1=0.f, sum2=0.f, P;
	//int i;
	//for(i=0;i<4;i++){
	//	P = expf(a[i]*powf(x[0]-x0[i],2.f)+b[i]*(x[0]-x0[i])*(x[1]-y0[i])+c[i]*powf(x[1]-y0[i],2.f));
	//	sum1+=A[i]*P*(2.f*a[i]*(x[0]-x0[i])+b[i]*(x[1]-y0[i]));
	//	sum2+=A[i]*P*(b[i]*(x[0]-x0[i])+2.f*c[i]*(x[1]-y0[i]));
	//}
	//r[0] = sum1;
	//r[1] = sum2;

	////3d atom potential field
	//float x2,x3,y3;
	//float r12,r13,r23;
	//x2 = x[0];
	//x3 = x[1];
	//y3 = x[2];
	//r12 = x2;
	//r13 = sqrtf(x3*x3+y3*y3);
	//r23 = sqrtf((x3-x2)*(x3-x2)+y3*y3);
	//r[0] = 12.f/powf(r12,8.f)*(1.f-powf(1.f/r12,6.f))*x2+12.f/powf(r23,8.f)*(1-powf(1.f/r23,6.f))*x2;
	//r[1] = 12.f/powf(r13,8.f)*(1.f-powf(1.f/r13,6.f))*x3+12.f/powf(r23,8.f)*(1-powf(1.f/r23,6.f))*x3;
	//r[2] = 12.f/powf(r13,8.f)*(1.f-powf(1.f/r13,6.f))*y3+12.f/powf(r23,8.f)*(1-powf(1.f/r23,6.f))*y3;

	////gradient of 2d test saddle function
	//r[0] = 2.f*sinf(x[0]*x[0]/2.f-powf(x[1],3.f)/4.f+3.f)+
	//	x[0]*cosf(x[0]*x[0]/2.f-powf(x[1],3.f)/4.f+3.f)*
	//	cosf(expf(x[1])-2.f*x[0]-1.f);
	//r[1] = -sinf(x[0]*x[0]/2.f-powf(x[1],3.f)/4.f+3.f)*expf(x[1])*
	//	sinf(expf(x[1])-2.f*x[0]-1.f)-3.f/4.f*
	//	cosf(x[0]*x[0]/2.f-powf(x[1],3.f)/4.f+3.f)*
	//	cosf(expf(x[1])-2.f*x[0]-1.f);

	//int i;
	//float sum = 0.f;
	//for(i=0;i<D;i++)
	//	sum+=x[i];
	//for(i=0;i<D;i++)
	//	r[i] = -sum+100.f*x[i]+powf(x[i],2)-powf(x[i],3);

	////Neuro
	//r[0] = x[0]*x[0]+x[2]*x[2]-1.f;
	//r[1] = x[1]*x[1]+x[3]*x[3]-1.f;
	//r[2] = x[4]*powf(x[2],3)+x[5]*powf(x[3],3);
	//r[3] = x[4]*powf(x[0],3)+x[5]*powf(x[1],3);
	//r[4] = x[4]*x[0]*x[2]*x[2]+x[5]*x[3]*x[3]*x[1];
	//r[5] = x[4]*x[0]*x[0]*x[2]+x[5]*x[1]*x[1]*x[3];

	////Chemical
	//float R = 10.f;
	//float R5 = 0.193f;
	//float R6 = 0.002597f/sqrtf(40.f);
	//float R7 = 0.003448f/sqrtf(40.f);
	//float R8 = 0.00001799f/sqrtf(40.f);
	//float R9 = 0.0002155f/sqrtf(40.f);
	//float R10 = 0.00003846f/sqrt(40.f);
	//r[0] = x[0]*x[1]+x[0]-3.f*x[4];
	//r[1] = 2.f*x[0]*x[1]+x[0]+x[1]*x[2]*x[2]+R8*x[1]-R*x[4]+
	//	2.f*R10*x[1]*x[1]+R7*x[1]*x[2]+R9*x[1]*x[3];
	//r[2] = 2.f*x[1]*x[2]*x[2]+2*R5*x[2]*x[2]-8.f*x[4]+R6*x[2]+R7*x[1]*x[2];
	//r[3] = R9*x[1]*x[3]+2.f*x[3]*x[3]-4.f*R*x[4];
	//r[4] = x[0]*(x[1]+1.f)+R10*x[1]*x[1]+x[1]*x[2]*x[2]+R8*x[1]+
	//	R5*x[2]*x[2]+x[3]*x[3]-1.f+R6*x[2]+R7*x[1]*x[2]+R9*x[1]*x[3];

	////Interval Arithmetic
	//r[0] = x[0]-0.25428722f-0.18324757f*x[3]*x[2]*x[8];
	//r[1] = x[1]-0.37842197f-0.16275449f*x[0]*x[9]*x[5];
	//r[2] = x[2]-0.27162577f-0.16955071f*x[0]*x[1]*x[9];
	//r[3] = x[3]-0.19807914f-0.15585316f*x[6]*x[0]*x[5];
	//r[4] = x[4]-0.44166728f-0.19950920f*x[6]*x[5]*x[2];
	//r[5] = x[5]-0.14654113f-0.18922793f*x[7]*x[4]*x[9];
	//r[6] = x[6]-0.42937161f-0.21180486f*x[1]*x[4]*x[7];
	//r[7] = x[7]-0.07056438f-0.17081208f*x[0]*x[6]*x[5];
	//r[8] = x[8]-0.34504906f-0.19612740f*x[9]*x[5]*x[7];
	//r[9] = x[9]-0.42651102f-0.21466544f*x[3]*x[7]*x[0];

	////Economics
	//int i,j,k;
	//float sum;
	//for(k=0;k<D-1;k++){
	//	sum = 0.f;
	//	for(i=0;i<D-k-1;i++)
	//		sum+=x[i]*x[i+k];
	//	r[k] = (x[k]+sum)*x[D-1];
	//}
	//sum = 0.f;
	//for(j=0;j<D-1;j++)
	//	sum+=x[j];
	//r[D-1] = sum+1.f;

	////More
	//float sum;
	//int i,j;
	//for(i=0;i<D;i++){
	//	sum = 0.f;
	//	for(j=0;j<D;j++)
	//		sum+=cosf(x[j]);
	//	r[i] = D-sum+(i+1)*(1.f-cosf(x[i]))-sinf(x[i]);
	//}
}

__device__ __host__ void matMul(float *A, float *B, float* C, int rA, int cA, int rB, int cB){
	//matrix multiplication with 1D pointer as input
	int i,j,k;
	for(i=0;i<rA;i++){
		for(j=0;j<cB;j++){
			C[i*cB+j] = 0.f;
			for(k=0;k<cA;k++)
				C[i*cB+j]+=A[i*cA+k]*B[k*cB+j];
		}
	}
}

__device__ __host__ void matTran(float *A, float *B, int rA, int cA){
	//matrix transpose with 1D pointer as input
	int i,j;
	for(i=0;i<rA;i++){
		for(j=0;j<cA;j++)
			B[j*rA+i] = A[i*cA+j];
	}
}

__device__ __host__ void matLU(float *A, float *L, float* U){
	//LU decomposition of a square matrix using Doolittle algorithm
	//matrices are treated as 1D pointers here
	int i,j,k,r;
	float sum;
	const int dim = D;

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

__device__ __host__ void matInv(float * A, float * iA){
	//matrix inversion using LU decomposition, 1D pointer
	//is used to represent the inversion matrix
	const int n = D;
	float L[n][n], U[n][n];
	float x[n], y[n], b[n];
	int i,j,k,r;
	float sum;

	matLU(A,(float*)L,(float*)U);

	//initialization
	for(i=0;i<n;i++){
		b[i] = 0.f;
		x[i] = 0.f;
		y[i] = 0.f;
	}

	//acquire each column of A by solving linear systems
	for(i=0;i<n;i++){
		b[i] = 1.f;
		//Ly=b
		for(k=0;k<n;k++){
			sum = 0.f;
			for(r=0;r<=k-1;r++)
				sum+=L[k][r]*y[r];
			y[k] = b[k]-sum;
		}
		//Ux=y
		for(k=n;k!=0;k--){
			sum = 0.f;
			for(r=k+1;r<=n;r++)
				sum+=U[k-1][r-1]*x[r-1];
			x[k-1] = (y[k-1]-sum)/U[k-1][k-1];
		}
		b[i] = 0.f;
		//ith column of the inversion
		for(j=0;j<n;j++)
			iA[j*n+i] = x[j];
	}
}

//point-to-point mapping xd->xc with one step
__device__ __host__ void p_map_1step(float * xc, const float * xd, const float * lb, const float * ub, const int * N){
	//Newton's method and armijio condition for pointwise dynamical system
	//The system is driven from xd to xc, namely, xd->xc
	float fx[O], fxp[O];
	float J[O][D], iJ[D][O];
	float v[D], h[D], xp[D];
	int i,j;
	float norm, hmax, hmin, t, nf, nfp;
	bool flag;

	for(i=0;i<D;i++){
		h[i] = (ub[i]-lb[i])/N[i];
		xp[i] = xd[i];
	}

	//construct local Jacobian matrix thru perturbation
	f(fx,xd);
	for(j=0;j<D;j++){
		xp[j]+=0.0001f*h[j];
		f(fxp,xp);
		for(i=0;i<O;i++)
			J[i][j] = -(fxp[i]-fx[i])/(0.0001f*h[j]);
		xp[j] = xd[j]; //recover for next directions's perturbation
	}

	//calculate normalized evolution direction
	matInv((float*)J,(float*)iJ);
	matMul((float*)iJ,(float*)fx,(float*)v,D,O,O,1);
	//norm = 0.f;
	//for(i=0;i<D;i++)
	//	norm+=v[i]*v[i];
	//norm = sqrtf(norm);
	//for(i=0;i<D;i++)
	//	v[i]/=norm;

	//calculate step size
	hmax = h[0];
	hmin = h[0];
	for(i=1;i<D;i++){
		if(h[i]>hmax)
			hmax = h[i];
		if(h[i]<hmin)
			hmin = h[i];
	}
	t = hmax;
	////all(|fnew|<|fold|)
	//while(true){
	//	flag = true;
	//	for(i=0;i<D;i++)
	//		xp[i] = xd[i]+v[i]*t;
	//	f(fxp,xp);
	//	for(i=0;i<O;i++){
	//		if(fabsf(fxp[i])>fabsf(fx[i])){
	//			flag = false;
	//			t/=2;
	//			break;
	//		}
	//	}
	//	if(flag || t<=hmin/2)
	//		break;
	//}

	//|fnew|<|fold|
	while(true){
		flag = true;
		for(i=0;i<D;i++)
			xp[i] = xd[i]+v[i]*t;
		f(fxp,xp);
		nf = 0.f;
		nfp = 0.f;
		for(i=0;i<D;i++){
			nf+=fx[i]*fx[i];
			nfp+=fxp[i]*fxp[i];
		}
		if(nfp>nf){
			flag = false;
			t/=1.25f;
		}
		if(flag || t<=hmin/2)
			break;
	}

	//pointwise mapping
	for(i=0;i<D;i++)
		xc[i] = xd[i]+v[i]*t;
}

//n step mapping forward xd->xc
__device__ __host__ void p_map(float * xc, const float * xd, const float * lb, const float * ub, const int * N, int st){
	int i,j;
	float xold[D],xnew[D];

	for(i=0;i<D;i++)
		xold[i] = xd[i];

	//shoot forward st times for point mapping
	for(i=0;i<st;i++){
		p_map_1step(xnew,xold,lb,ub,N);
		for(j=0;j<D;j++)
			xold[j] = xnew[j];
	}

	for(i=0;i<D;i++)
		xc[i] = xnew[i];
}


__device__ __host__ unsigned long long ztocell(const int * z, const int * N){
	int i;
	unsigned long long ncell;
	int temp[D];

	for(i=0;i!=D;i++)
		temp[i]=z[i];
	ncell = temp[0];
	
	unsigned long long b = N[0];
	for(i=0;i!=D;i++)
		temp[i]-=1;
	for(i=1;i!=D;i++){
		ncell+=temp[i]*b;
		b*=N[i];
	}
	//ncell+=1;
	return ncell;
}

__device__ __host__ void xtoz(int * z, const float *x, const float *h, const float * lb){
	int i;
	for(i=0;i!=D;i++)
		z[i] = (int)rintf((x[i]-lb[i])/h[i] + 0.5f);
}

__device__ __host__ void ztox(float * xd, const int * z, const float * h, const float * lb){
	int i;
	for(i=0;i!=D;i++)
		xd[i] = lb[i]+h[i]*z[i]-0.5f*h[i];
}

__device__ __host__ void celltoz(int * z, const unsigned long long cs, const int * N){
	int coord[D];
	int i;
	unsigned long long cell=cs;
	cell-=1;
	for(i=0;i!=D;i++){
		coord[i] = cell%N[i] + 1;
		cell = (int)floorf(float(cell)/float(N[i]));
		z[i] = coord[i];
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

__device__ __host__ void subdivision(unsigned long long *rcells, const unsigned long long cell, const int * N, int * N_new, 
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
	unsigned long long cs,cr;
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

__device__ __host__ void neighbour_finder(unsigned long long * neighbour, const int * N, unsigned long long cs){
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
		if(goodcell)
			neighbour[i] = ztocell(z_current,N);
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

__device__ __host__ void fminsearch(const float *v0, float *vf){
	//Nelder-Mead simplex algorithm for single objective optimization with
	//several variables, this code implements the same function in Matlab
	const int Dp1=D+1;
	float x0[D][Dp1];
	float x[D][Dp1];
	float xshrink[D][Dp1];
	float f[Dp1];
	float temp;
	int index[Dp1];
	float xj[D],r[D],m[D],s[D],c[D],cc[D],xtemp[D];
	float infnorm[D];
	float sum;
	int i,j,itemp;
	float tolx,tolf,min_tolx,min_tolf;
	int iter,maxiter;
	float fx1,fr,fxn,fxnp1,fs,fc,fcc;
	bool shrink;

	//generate the first point set near initial guessing
	for(i=0;i<D;i++)
		for(j=0;j<Dp1;j++)
			x0[i][j] = v0[i]; //every column is identical
	for(j=1;j<Dp1;j++){
		if(v0[j-1]!=0.f)
			x0[j-1][j]*=1.05f;
		else
			x0[j-1][j] = 0.00025f;
	}

	//sort the n+1 points with their function values in ascending order
	for(j=0;j<Dp1;j++){
		for(i=0;i<D;i++)
			xj[i] = x0[i][j];
		f[j] = fp(xj);
	}
	for(i=0;i<Dp1;i++)
		index[i] = i;
	//bubble sort ascending order
	for(i=0;i<Dp1;i++){
		for(j=0;j<Dp1-i-1;j++){
			if(f[j]>f[j+1]){
				temp = f[j];
				f[j] = f[j+1];
				f[j+1] = temp;
				itemp = index[j];
				index[j] = index[j+1];
				index[j+1] = itemp;
			}
		}
	}
	//rearrange x point set (matrix) at ascending order
	for(j=0;j<Dp1;j++)
		for(i=0;i<D;i++)
			x[i][j] = x0[i][index[j]];

	//Main algorithm
	min_tolx = 0.0001f;
	min_tolf = 0.0001f;
	maxiter = 5000;
	iter = 0;
	while(iter<maxiter){
		//function value tolerence
		tolf = fabsf(f[1]-f[0]);
		for(i=2;i<Dp1;i++){
			if(fabsf(f[i]-f[0])>tolf)
				tolf = fabsf(f[i]-f[0]);
		}

		//x tolerence using max inf norm
		for(j=0;j<D;j++){
			infnorm[j]=fabsf(x[0][j+1]-x[0][0]);
			for(i=1;i<D;i++){
				if(fabsf(x[i][j+1]-x[i][0])>infnorm[j])
					infnorm[j] = fabsf(x[i][j+1]-x[i][0]);
			}
		}
		tolx = infnorm[0];
		for(i=1;i<D;i++){
			if(infnorm[i]>tolx)
				tolx = infnorm[i];
		}

		//breaking criteria
		if(tolx<min_tolx && tolf<min_tolf)
			break;

		//generate reflected point
		for(i=0;i<D;i++){
			sum = 0.f;
			for(j=0;j<D;j++)
				sum+=x[i][j];
			m[i] = sum/D;
			r[i] = 2.f*m[i]-x[i][D];
		}
		fr = fp(r);

		for(i=0;i<D;i++)
			xtemp[i] = x[i][0];
		fx1 = fp(xtemp);
		for(i=0;i<D;i++)
			xtemp[i] = x[i][D-1];
		fxn = fp(xtemp);

		shrink = false;

		if(fr<fx1){
			//calcuate expansion point
			for(i=0;i<D;i++)
				s[i] = m[i]+2.f*(m[i]-x[i][D]);
			fs = fp(s);
			if(fs<fr){
				//accept s
				for(i=0;i<D;i++)
					x[i][D] = s[i];
				f[D] = fs;
			}
			else{
				//accept r
				for(i=0;i<D;i++)
					x[i][D] = r[i];
				f[D] = fr;
			}
		}
		else{ //fx1<=fr
			if(fr<fxn){
				//fx1<=fr<fv, accept r
				for(i=0;i<D;i++)
					x[i][D] = r[i];
				f[D] = fr;
			}
			else{ //fxr>=fxn
				//perform contraction
				for(i=0;i<D;i++)
					xtemp[i] = x[i][D];
				fxnp1 = fp(xtemp);
				if(fr<fxnp1){
					//outside contraction
					for(i=0;i<D;i++)
						c[i] = m[i]+0.5f*(r[i]-m[i]);
					fc = fp(c);
					if(fc<fr){
						//accept c
						for(i=0;i<D;i++)
							x[i][D] = c[i];
						f[D] = fc;
					}
					else
						shrink = true;
				}
				else{
					//inside contraction
					for(i=0;i<D;i++)
						cc[i] = m[i]+(x[i][D]-m[i])*0.5f;
					fcc = fp(cc);
					if(fcc<fxnp1){
						//accept cc
						for(i=0;i<D;i++)
							x[i][D] = cc[i];
						f[D] = fcc;
					}
					else
						shrink = true;
				}

				//shrink
				if(shrink){
					for(i=0;i<D;i++)
						xshrink[i][0] = x[i][0];
					for(j=1;j<Dp1;j++)
						for(i=0;i<D;i++)
							xshrink[i][j] = x[i][0]+0.5f*(x[i][j]-x[i][0]);
					for(j=0;j<Dp1;j++)
						for(i=0;i<D;i++)
							x[i][j] = xshrink[i][j];
					for(j=1;j<Dp1;j++){
						for(i=0;i<D;i++)
							xtemp[i] = x[i][j];
						f[j] = fp(xtemp);
					}
				}
			}
		}

		//we now have a new point set x and f, sort them again
		for(i=0;i<D;i++)
			for(j=0;j<Dp1;j++)
				x0[i][j] = x[i][j];
		for(i=0;i<Dp1;i++)
			index[i] = i;
		//bubble sort ascending order
		for(i=0;i<Dp1;i++){
			for(j=0;j<Dp1-i-1;j++){
				if(f[j]>f[j+1]){
					temp = f[j];
					f[j] = f[j+1];
					f[j+1] = temp;
					itemp = index[j];
					index[j] = index[j+1];
					index[j+1] = itemp;
				}
			}
		}
		//rearrange x point set at ascending order
		for(j=0;j<Dp1;j++)
			for(i=0;i<D;i++)
				x[i][j] = x0[i][index[j]];

		//next iteration
		iter++;
	}

	//final result is the first point at the final point set
	for(i=0;i<D;i++)
		vf[i] = x[i][0];
}