/*
SCM boundary searching, SCM is performed
by sweeping the whole domain without
subdivision.
By: Free Xiong; 2015-03-05
*/
#include "cuda_runtime.h"
#include <cstdlib>
#include "cell.cuh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>

using namespace std;

#define NumberOfBlock 256
#define ThreadPerBlock 512

//sop based scm for stability boundary finding
__global__ void scm_sop(long tot_cells, int *N, float *lb, float *ub, unsigned long long *img){
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	long i,j;
	float h[D],xcs[D],xcn[D],ximg[D],gradf[D],norm;
	float fcs,fcn,df;
	int zcs[D],zcn[D],zimg[D];
	unsigned long long neighbours[pop];
	unsigned long long cs,cn;
	bool flag,boundary;
	float t,hmax,hmin;

	for(j=0;j<D;j++)
		h[j] = (ub[j]-lb[j])/(float)N[j];

	hmax = h[0];
	hmin = h[0];
	for(j=1;j<D;j++){
		if(h[j]>hmax)
			hmax = h[j];
		if(h[j]<hmin)
			hmin = h[j];
	}

	i = id;
	while(i<tot_cells){
		cs = i+1;
		celltoz(zcs,cs,N);
		ztox(xcs,zcs,h,lb);
		fcs = fp(xcs);
		flag = true;
		df = 1000.f;

		//examine neighbours, pick the steepest descent as image cell
		neighbour_finder(neighbours,N,cs);
		for(j=0;j<pop;j++){
			if(neighbours[j]==0)
				continue;
			cn = neighbours[j];
			celltoz(zcn,cn,N);
			ztox(xcn,zcn,h,lb);
			fcn = fp(xcn);
			if(fcn<fcs && fcn-fcs<df){
				df = fcn-fcs;
				img[i] = cn;
				flag = false;
			}
		}

		//we take the boundary cell as absorbing cell
		boundary = false;
		for(j=0;j<D;j++){
			if(zcs[j]==1 || zcs[j]==N[j]){
				boundary = true;
				break;
			}
		}
		if(boundary&&flag)
			img[i] = 0; //mark as sink cell
		else if(!boundary&&flag){
		//	//use gradient descend as secondary criteria
		//	f(gradf,xcs);
		//	norm = 0.f;
		//	for(j=0;j<D;j++)
		//		norm+=gradf[j]*gradf[j];
		//	norm = sqrtf(norm);

		//	t = hmax; //step size
		//	while(true){
		//		for(j=0;j<D;j++)
		//			ximg[j] = xcs[j]-t*gradf[j]/norm;
		//		if(fp(ximg)<fcs){
		//			xtoz(zimg,ximg,h,lb);
		//			img[i] = ztocell(zimg,N);
		//			break;
		//		}
		//		else if(t<hmin/2.f){
		//			img[i] = cs; //absorbing
		//			break;
		//		}
		//		else if(fp(ximg)>=fcs && t>=hmin/2.f)
		//			t/=2.f;
		//	}

			//use fminsearch as secondary criteria
			fminsearch(xcs,ximg);
			xtoz(zimg,ximg,h,lb);
			img[i] = ztocell(zimg,N);
		}

		i+=blockDim.x*gridDim.x;
	}
}

//sequential scm search with gr, pe and st arrays
void scm_unravel(unsigned long long *img, long tot_cells, int *gr, int *pe, int*st){
	unsigned long long cell_new,cell_old;
	int i,j,k,m;
	int gr1,pe1,st1;
	int g = 1;
	std::vector<unsigned long long> path;
	bool flag;

	for(i=0;i<tot_cells;i++){
		if(gr[i]!=0 && gr[i]!=-1)
			continue;

		//generate a path from the current virgin cell
		cell_old = i+1;
		flag = true;
		while(flag){
			//process sink cell first
			if(cell_old==0){
				for(j=0;j<path.size();j++){
					gr[path[j]-1] = 1;
					st[path[j]-1] = 0;
					pe[path[j]-1] = 1;
				}
				path.clear();
				break;
			}

			switch(gr[cell_old-1]){
			case 0:
				//virgin cell ahead, keep exploring
				gr[cell_old-1]=-1;
				path.push_back(cell_old);
				cell_new = img[cell_old-1];
				cell_old = cell_new;
				break;
			case -1:
				//cyclic structure, new group found
				g++;
				for(j=path.size()-1,k=1;j>=0;j--,k++){
					if(cell_new!=path[j])
						st[path[j]-1] = 0; //zero step number for cyclic cells
					else
						break;
				}
				for(j=0;j<path.size();j++){
					gr[path[j]-1] = g;
					pe[path[j]-1] = k;
				}
				for(j=path.size()-k-1,m=1;j>=0;j--,m++)
					st[path[j]-1] = m;
				
				path.clear();
				flag = false;
				break;
			default:
				//merge to another group
				gr1 = gr[cell_new-1];
				st1 = st[cell_new-1];
				pe1 = pe[cell_new-1];
				for(j=path.size()-1,k=1;j>=0;j--,k++){
					gr[path[j]-1] = gr1;
					pe[path[j]-1] = pe1;
					st[path[j]-1] = st1+k;
				}
				
				path.clear();
				flag = false;
				break;
			}
		}
	}
}

//extract boundary in cell space from a saddle cell
void boundary(int *gr, int *N, int *bn_old, long tot_cells, unsigned long long cell){
	int *tgt_cells;
	int *bn_new;
	long i,j,k;
	unsigned long long ncells[pop], cs, nc;

	bn_new = new int[tot_cells];
	tgt_cells = new int[tot_cells];
	memset(bn_new,0,tot_cells*sizeof(int));
	memset(tgt_cells,0,tot_cells*sizeof(int));
	
	tgt_cells[cell-1] = 1;
	bn_new[cell-1] = 1;
	bn_old[cell-1] = 1; //bn_old must be initialized with all zeros

	//continuation like extraction
	while(true){
		for(i=0;i<tot_cells;i++){
			if(tgt_cells[i]==0)
				continue;
			cs = i+1;
			neighbour_finder(ncells,N,cs);
			for(j=0;j<pop;j++){
				nc = ncells[j];
				if(nc==0)
					continue;
				if(gr[cs-1]!=gr[nc-1] && bn_old[nc-1]==0){
					//new boundary cell brought in
					tgt_cells[nc-1] = 1;
					bn_new[nc-1] = 1;
				}
			}
			//remove cs from the target cell set since it's already been processed
			tgt_cells[cs-1] = 0;
		}

		//check whether steay state reaches
		for(i=0;i<tot_cells;i++)
			if(bn_new[i]!=bn_old[i])
				break;
		
		if(i==tot_cells)
			break;
		else
			memcpy(bn_old,bn_new,tot_cells*sizeof(int));
	}

	delete[] tgt_cells;
	delete[] bn_new;
}

void saveSCM(int *gr, int *pe, int *st, int *bd, long tot_cells, int*N){
	long i;
	ofstream outData;
	//outData.open("SCMproperties.dat"); //store gr, pe, st
	//for(i=0;i<tot_cells;i++)
	//	outData<<gr[i]<<"   "<<pe[i]<<"   "<<st[i]<<"   "<<bd[i]<<endl;
	//outData.close();
	outData.open("boundary_cells.dat");
	for(i=0;i<tot_cells;i++)
		if(bd[i]==1)
			outData<<i+1<<endl;
	outData.close();
	outData.open("SCMpartition.dat");
	for(i=0;i<D;i++)
		outData<<N[i]<<endl;
	outData.close();
}

void pscm_bd(int *N, float *lb, float *ub){
	/*
	parallel analysis of simple cell mapping, the mapping construction
	is with sop approach, which is implemented in parallel. scm unravelling
	is conducted with Prof. Hsu's traditional sequential approach
	*/
	int *device_N,num_sd,*bd,*bd_new;
	fstream file;
	vector<vector<float>> saddles;
	vector<unsigned long long> mcells;
	vector <float> rowVector(D);
	float *device_lb, *device_ub;
	int *gr,*pe,*st,z[D];
	float x[D],h[D];
	unsigned long long *img, *device_img,*scells;
	long tot_cells=1,i,j;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);
	clock_t cpu_time_start, cpu_time_end;
	float cpu_time;

	for(i=0;i<D;i++)
		tot_cells*=N[i];

	cout<<"Performing SCM sweeping..."<<endl;

	cudaMalloc(&device_N, D*sizeof(int));
	cudaMalloc(&device_lb, D*sizeof(float));
	cudaMalloc(&device_ub, D*sizeof(float));	
	cudaMemcpy(device_N, N, D*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_lb, lb, D*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_ub, ub, D*sizeof(int), cudaMemcpyHostToDevice);

	//---------------------------build scm in parallel------------------------------------
	cudaMalloc(&device_img, tot_cells*sizeof(unsigned long long));
	scm_sop<<<NumberOfBlock,ThreadPerBlock>>>(tot_cells,device_N,device_lb,device_ub,device_img);
	cudaDeviceSynchronize();
	img = new unsigned long long[tot_cells];
	cudaMemcpy(img,device_img,tot_cells*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
	cudaFree(device_img);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float gpu_time;
	cudaEventElapsedTime(&gpu_time, start, stop);

	//------------------------sequential scm unravelling----------------------------------
	cpu_time_start = clock();
	gr = new int[tot_cells];
	pe = new int[tot_cells];
	st = new int[tot_cells];
	memset(gr,0,tot_cells*sizeof(int));
	memset(pe,0,tot_cells*sizeof(int));
	memset(st,0,tot_cells*sizeof(int));
	scm_unravel(img,tot_cells,gr,pe,st);

	//-------------------------extract boundary from scm----------------------------------
	//read the saddle points from matlab generated file
	//read cells captured by ga from the file
	num_sd = 0;
	file.open("saddles.dat");
	if(file.is_open()){
		while(file.good()){
			saddles.push_back(rowVector);
			for(i=0;i<D;i++)
				file >> saddles[num_sd][i];
			num_sd++;
		}
	}
	saddles.erase(saddles.end()-1);
	
	num_sd--;
	scells = (unsigned long long*)malloc(num_sd*sizeof(unsigned long long));
	for(i=0;i<D;i++)
		h[i] = (ub[i]-lb[i])/N[i];
	for(i=0;i<num_sd;i++){
		for(j=0;j<D;j++)
			x[j] = saddles[i][j];
		xtoz(z,x,h,lb);
		scells[i] = ztocell(z,N);
	}
	saddles.clear();

	cout<<endl;
	cout<<"Extracting boundaries..."<<endl;
	bd = new int[tot_cells];
	bd_new = new int[tot_cells];
	memset(bd,0,tot_cells*sizeof(int));
	memset(bd_new,0,tot_cells*sizeof(int));
	
	//find boundary start from all saddle cells
	for(i=0;i<num_sd;i++){
		boundary(gr,N,bd_new,tot_cells,scells[i]);
		//elementwise operation
		for(j=0;j<tot_cells;j++)
			if(bd_new[j]==1 && bd[j]==0)
				bd[j] = 1;
		memset(bd_new,0,tot_cells*sizeof(int));
	}

	cpu_time_end = clock();
	cpu_time = float(cpu_time_end-cpu_time_start)/CLOCKS_PER_SEC*1000;
	cout<<endl;
	cout<<"Total SCM analysis runtime is: "<<setw(8)<<gpu_time+cpu_time<<setprecision(6)<<" ms"<<endl;

	//---------------save global properties for post processing-------------------
	saveSCM(gr,pe,st,bd,tot_cells,N);

	delete[] bd;
	delete[] bd_new;
	delete[] scells;
	delete[] img;
	delete[] gr;
	delete[] pe;
	delete[] st;

	cudaFree(device_N);
	cudaFree(device_lb);
	cudaFree(device_ub);
}

int main(){
	int Nscm[D] = {120,120,120};
	float lb[D] = {-2.5f,-2.5f,-2.5f};
	float ub[D] = {2.5f,2.5f,2.5f};

	pscm_bd(Nscm,lb,ub); //pure scm for global analysis on boundary finding

	return 0;
}